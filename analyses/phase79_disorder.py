#!/usr/bin/env python3
"""phase79_disorder.py -- P16.13 : couche `disorder` (désordre intrinsèque / IDR).

Ajoute une PROPRIÉTÉ biophysique orthogonale : la fraction de désordre intrinsèque par protéine et ses régions
désordonnées (IDR). Deux signaux :
  - SÉQUENCE (primaire) : metapredict v3 (Holehouse lab) — score de désordre par résidu + frontières d'IDR.
  - STRUCTURE (corroborant) : mean pLDDT AlphaFold déjà en base (un pLDDT bas corrèle avec le désordre, Akdel 2022)
    — mais on ne l'a qu'en MOYENNE, donc secondaire ; metapredict tranche par résidu.

APPORT HONNÊTETÉ (la vraie valeur pour les dark) : un gène `dark` très DÉSORDONNÉ est dark **PAR CONSTRUCTION**
— une protéine sans repli globulaire n'a rien à matcher pour Foldseek/homologie, donc son obscurité est ATTENDUE,
pas un échec du pipeline. On sépare ainsi « dark car orphelin-IDR » de « dark car homologue pas encore trouvé ».
Flag `idr_orphan` posé sur les dark à forte fraction de désordre.

GARDE-FOU (piste P16.13) : la couche annote une PROPRIÉTÉ, pas une fonction. AUCUN verdict n'est modifié. Aucun
appel à la séparation de phase (LLPS/condensat) : le désordre ne prouve PAS la LLPS (il faudrait une évidence
indépendante — liaison ADN, littérature LLPS expérimentale). On ne l'écrit donc pas.

Écrit EN PLACE. Run: /tmp/venv_disorder/bin/python analyses/phase79_disorder.py   (metapredict dans le venv)
"""
from __future__ import annotations
import glob, json, sys
from pathlib import Path

try:
    import metapredict as mp
except ImportError:
    sys.exit("metapredict absent — lancer avec le venv : /tmp/venv_disorder/bin/python analyses/phase79_disorder.py")

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
THRESH = 0.5          # score metapredict > 0.5 = résidu désordonné
HIGH = 0.40           # fraction de désordre au-dessus de laquelle on parle de protéine "très désordonnée"
PARTIAL = 0.15
SOURCE = "metapredict v3 (Emenecker/Holehouse) per-residue disorder + AlphaFold mean pLDDT (annotation_mtbc P16.13)"


def classify(frac: float) -> str:
    if frac >= HIGH:
        return "highly disordered"
    if frac >= PARTIAL:
        return "partially disordered"
    return "ordered"


def main() -> None:
    files = sorted(glob.glob(str(GENES / "*.json")))
    n = 0
    cls_count = {"ordered": 0, "partially disordered": 0, "highly disordered": 0}
    idr_orphans = []
    dark_total = 0
    for f in files:
        d = json.loads(Path(f).read_text())
        seq = d.get("protein_mtbc0") or ""
        if not seq or len(seq) < 20:
            continue
        scores = mp.predict_disorder(seq)
        frac = round(sum(1 for s in scores if s > THRESH) / len(scores), 3)
        try:
            doms = mp.predict_disorder_domains(seq).disordered_domain_boundaries
        except Exception:
            doms = []
        regions = [[int(a), int(b)] for a, b in doms]
        longest = max((b - a for a, b in regions), default=0)
        total_idr = sum(b - a for a, b in regions)
        mean_plddt = (d.get("plddt_af") or {}).get("mean_plddt")
        cls = classify(frac)
        cls_count[cls] += 1
        verdict = d.get("verdict")
        if verdict == "dark":
            dark_total += 1
        # idr_orphan = CONVERGENCE des DEUX signaux (anti-sur-appel) : dark ET séquence très désordonnée
        # (metapredict >= HIGH) ET structure NON confiante (mean pLDDT < 70). Si AlphaFold replie la protéine
        # avec confiance (pLDDT >= 70), metapredict sur-appelle -> on NE flague PAS (le désordre n'explique pas
        # l'obscurité). Mesuré : ce gate retire 21 faux positifs (100 -> 79), dont Rv2541 (vulnérable), Rv0810c.
        seq_dis = frac >= HIGH or longest > 0.5 * len(seq)
        struct_dis = mean_plddt is None or mean_plddt < 70
        idr_orphan = bool(verdict == "dark" and seq_dis and struct_dis)
        if idr_orphan:
            idr_orphans.append((d["rv"], d.get("gene") or "-", frac, len(seq)))

        if idr_orphan:
            interp = ("both sequence (metapredict) and structure (low AlphaFold pLDDT) indicate little stable "
                      "globular structure: this gene is dark BY CONSTRUCTION — there is no confident fold for "
                      "structure/homology search to match, so its obscurity is expected, not a pipeline failure")
        elif seq_dis and not struct_dis:
            interp = ("sequence-based disorder is high but AlphaFold folds it confidently (pLDDT>=70): treat the "
                      "disorder call with caution (possible metapredict over-call)")
        elif cls == "ordered":
            interp = "predominantly ordered / globular"
        else:
            interp = (f"carries a substantial disordered region ({total_idr}/{len(seq)} residues); "
                      "disorder is a property, not a function")

        d["disorder"] = {
            "frac_disordered": frac,
            "classification": cls,
            "n_idr": len(regions),
            "longest_idr": longest,
            "total_idr_residues": total_idr,
            "idr_regions": regions[:8],
            "mean_plddt": mean_plddt,
            "idr_orphan": idr_orphan,
            "interpretation": interp,
            "caveat": "A property (biophysics), not a function. No LLPS/condensate claim is made from disorder alone. "
                      "Verdict unchanged.",
            "source": SOURCE,
        }
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
        if n % 500 == 0:
            print(f"  {n} fiches…")

    print(f"\ncouche `disorder` écrite sur {n} fiches.")
    print(f"  classification : {cls_count}")
    print(f"\n>>> DARK expliqués par le désordre (idr_orphan) : {len(idr_orphans)} / {dark_total} dark")
    print("    = leur obscurité est ATTENDUE (pas de repli à matcher), pas un échec du pipeline :")
    for rv, nm, frac, L in sorted(idr_orphans, key=lambda x: -x[2])[:20]:
        print(f"      {rv:9} {nm:9} {int(frac*100):3}% désordre, {L} aa")


if __name__ == "__main__":
    main()
