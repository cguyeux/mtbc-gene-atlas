#!/usr/bin/env python3
"""phase18b_deeptmhmm.py -- upgrade de la couche localisation vers DeepTMHMM (référence) (P3.7).

phase18 livrait un substitut basse résolution (Kyte-Doolittle). DeepTMHMM (Hallgren et al. 2022) est le
prédicteur de RÉFÉRENCE (peptide signal + hélices TM + topologie), exécutable en cloud biolib ANONYME
(sans licence ni compte). Ce script parse sa sortie `.3line` et REMPLACE la couche `localization` par la
prédiction DeepTMHMM, en CONSERVANT le lipobox mycobactérien de phase18 (DeepTMHMM ne distingue pas les
lipoprotéines ; une lipoprotéine = type SP + lipobox).

Format .3line : 3 lignes/protéine — `>Rv | TYPE`, séquence, topologie (I=inside, M=TM, O=outside, S=signal).
TYPE ∈ {TM, SP, SP+TM, GLOB, BETA}.

Entrée : résultats/phase18_localization/predicted_topologies.3line + localization.json (pour le lipobox).
Sortie : réécrit localization.json (fusion) + fusionne `localization` dans les fiches.
Run: python analyses/phase18b_deeptmhmm.py [chemin_vers_3line] [--write]
"""
from __future__ import annotations
import argparse, glob, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
LOCDIR = ROOT / "résultats" / "phase18_localization"


def parse_3line(path: Path):
    out = {}
    lines = path.read_text().splitlines()
    i = 0
    while i + 2 < len(lines) + 1:
        if i >= len(lines) or not lines[i].startswith(">"):
            i += 1; continue
        header = lines[i]
        topo = lines[i + 2] if i + 2 < len(lines) else ""
        rv = header[1:].split("|")[0].strip()
        typ = header.split("|")[1].strip() if "|" in header else ""
        # compter les segments M (hélices TM)
        n_tm, in_m = 0, False
        for c in topo:
            if c == "M" and not in_m:
                n_tm += 1; in_m = True
            elif c != "M":
                in_m = False
        signal = "S" in topo or typ in ("SP", "SP+TM")
        out[rv] = {"type": typ, "tm_helices": n_tm, "signal_peptide": bool(signal)}
        i += 3
    return out


def summary(typ, n_tm, signal, lipo):
    if lipo:
        return "predicted lipoprotein (lipobox + signal peptide)"
    if typ == "BETA":
        return "predicted beta-barrel outer-membrane protein"
    if typ == "SP+TM":
        return f"predicted membrane protein with signal peptide ({n_tm} TM helix{'es' if n_tm > 1 else ''})"
    if typ == "SP":
        return "predicted secreted protein (signal peptide)"
    if typ == "TM":
        return f"predicted membrane protein ({n_tm} TM helix{'es' if n_tm > 1 else ''})"
    return "predicted cytoplasmic / globular"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("three_line", nargs="?", default=str(LOCDIR / "predicted_topologies.3line"))
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args(argv)

    dt = parse_3line(Path(a.three_line))
    print(f"DeepTMHMM : {len(dt)} protéines parsées")
    old = json.loads((LOCDIR / "localization.json").read_text()) if (LOCDIR / "localization.json").exists() else {}

    rec = {}
    for rv, d in dt.items():
        lipo = bool((old.get(rv) or {}).get("lipoprotein"))
        cys = (old.get(rv) or {}).get("lipobox_cys")
        rec[rv] = {"method": "DeepTMHMM", "type": d["type"], "tm_helices": d["tm_helices"],
                   "signal_peptide": d["signal_peptide"], "lipoprotein": lipo, "lipobox_cys": cys,
                   "prediction": summary(d["type"], d["tm_helices"], d["signal_peptide"], lipo)}
    (LOCDIR / "localization.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))

    from collections import Counter
    print("Types DeepTMHMM :", dict(Counter(v["type"] for v in rec.values())))
    print("lipoprotéines (lipobox) :", sum(1 for v in rec.values() if v["lipoprotein"]),
          "| signal peptide :", sum(1 for v in rec.values() if v["signal_peptide"]),
          "| >=1 TM :", sum(1 for v in rec.values() if v["tm_helices"] >= 1))

    if a.write:
        n = 0
        for f in glob.glob(str(GENES / "*.json")):
            d = json.load(open(f))
            if d["rv"] in rec:
                d["localization"] = rec[d["rv"]]
                Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False)); n += 1
        print(f"\n(--write) localization (DeepTMHMM) fusionné dans {n} fiches")


if __name__ == "__main__":
    raise SystemExit(main())
