#!/usr/bin/env python3
"""phase18_localization.py -- couche de localisation subcellulaire prédite (P3.5).

Les prédicteurs de référence (SignalP6 sous licence, DeepTMHMM en cloud biolib, tmhmm.py qui ne compile
pas sous Python 3.14) n'étant pas exploitables ici, on livre une couche fully-local, sans compilation ni
licence, fondée sur des méthodes ÉTABLIES et citables, clairement étiquetées « prédiction » :

  - Hélices transmembranaires : hydropathie de Kyte-Doolittle (fenêtre 19, seuil 1.6) — Kyte & Doolittle 1982.
    Un unique segment hydrophobe à l'extrême N-terminal (centre < 30 aa) est signalé comme signal/ancre
    plutôt que TM interne.
  - Lipoprotéines : motif lipobox de SPaseII dans les ~40 premiers résidus (Cys lipidé en +1),
    consensus (myco)bactérien [LVIFMW][ASTVILGCM][GAST]C — Sutcliffe & Harrington 2004.

C'est une couche indicative de basse résolution (hydropathie sur-prédit ; ne remplace pas SignalP6/DeepTMHMM),
utile pour flaguer les hypothétiques membranaires / sécrétés / lipoprotéines.

Sortie : résultats/phase18_localization/localization.json + fusion `localization` dans les fiches.
Stdlib pure. Run: python analyses/phase18_localization.py [--write]
"""
from __future__ import annotations
import argparse, glob, json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "résultats" / "phase18_localization"

KD = {  # Kyte & Doolittle 1982 hydropathy
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
WIN = 19
TM_THRESH = 1.6
LIPOBOX = re.compile(r"[LVIFMW][ASTVILGCM][GAST]C")


def tm_segments(seq: str):
    """Segments où la moyenne d'hydropathie (fenêtre WIN) dépasse le seuil TM. Renvoie liste de centres."""
    vals = [KD.get(c, 0.0) for c in seq]
    n = len(vals)
    if n < WIN:
        return []
    half = WIN // 2
    smooth = []
    for i in range(half, n - half):
        smooth.append((i, sum(vals[i - half:i + half + 1]) / WIN))
    segs, in_seg, start = [], False, 0
    for i, v in smooth:
        if v > TM_THRESH and not in_seg:
            in_seg, start = True, i
        elif v <= TM_THRESH and in_seg:
            in_seg = False
            segs.append((start, i - 1))
    if in_seg:
        segs.append((start, smooth[-1][0]))
    # centres, en fusionnant les segments proches (<5 aa d'écart)
    merged = []
    for s, e in segs:
        if merged and s - merged[-1][1] < 5:
            merged[-1] = (merged[-1][0], e)
        else:
            merged.append((s, e))
    return [ (s + e) // 2 for s, e in merged ]


def lipobox(seq: str):
    """Cherche un lipobox dans les 40 premiers résidus ; renvoie la position (1-based) du Cys lipidé, ou None."""
    for m in LIPOBOX.finditer(seq[:40]):
        cys = m.end()  # position du C (1-based = index de fin du match)
        if 12 <= cys <= 40:  # après un n/h-region plausible
            return cys
    return None


def classify(seq: str):
    centers = tm_segments(seq)
    lb = lipobox(seq)
    n_tm = len(centers)
    signal_anchor = (n_tm == 1 and centers[0] < 30)
    if lb:
        loc = "predicted lipoprotein (lipobox)"
    elif signal_anchor:
        loc = "predicted signal peptide / membrane anchor (N-terminal)"
    elif n_tm >= 1:
        loc = f"predicted membrane protein ({n_tm} TM helix{'es' if n_tm > 1 else ''})"
    else:
        loc = "no TM/signal predicted (likely cytoplasmic)"
    return {"tm_helices": n_tm, "tm_centers": centers[:12], "lipoprotein": bool(lb),
            "lipobox_cys": lb, "signal_anchor": signal_anchor, "prediction": loc}


def main(argv=None):
    ap = argparse.ArgumentParser(); ap.add_argument("--write", action="store_true")
    ap.add_argument("--force", action="store_true",
                    help="écraser localization.json même s'il porte déjà la référence DeepTMHMM")
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    # Garde-fou : ce script est le SUBSTITUT Kyte-Doolittle, SUPERSÉDÉ par phase18b_deeptmhmm.py
    # (référence). Les deux écrivent le MÊME localization.json. Le relancer après phase18b
    # régresserait silencieusement la couche live (et donc un phase4-full) du DeepTMHMM vers le KD.
    lj = OUT / "localization.json"
    if lj.exists() and not a.force:
        try:
            existing = json.loads(lj.read_text())
            if any((v or {}).get("method") == "DeepTMHMM" for v in existing.values()):
                print("REFUS : localization.json porte déjà la référence DeepTMHMM (phase18b).")
                print("        Relancer ce substitut KD la régresserait. Utiliser --force pour outrepasser.")
                return
        except Exception:
            pass
    rec = {}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f)); seq = d.get("protein_mtbc0") or ""
        if seq:
            rec[d["rv"]] = classify(seq)
    (OUT / "localization.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))
    from collections import Counter
    c = Counter(v["prediction"].split(" (")[0] for v in rec.values())
    print(f"localization.json : {len(rec)} gènes")
    for k, n in c.most_common():
        print(f"  {k}: {n}")
    lip = sum(1 for v in rec.values() if v["lipoprotein"])
    tm = sum(1 for v in rec.values() if v["tm_helices"] >= 1)
    print(f"lipoprotéines prédites : {lip} | avec >=1 TM : {tm}")
    # intersection hypothétiques
    hyp_lip = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        if "hypothetical" in (d.get("product_h37rv") or "").lower() and rec.get(d["rv"], {}).get("lipoprotein"):
            hyp_lip += 1
    print(f"hypothétiques prédits lipoprotéines : {hyp_lip}")

    if a.write:
        n = 0
        for f in glob.glob(str(GENES / "*.json")):
            d = json.load(open(f))
            if d["rv"] in rec:
                d["localization"] = rec[d["rv"]]
                Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False)); n += 1
        print(f"\n(--write) champ localization écrit dans {n} fiches")


if __name__ == "__main__":
    raise SystemExit(main())
