#!/usr/bin/env python3
"""phase58_ntm_conservation.py -- P10.2 : intègre le profil de conservation NTM à la couche `outgroup`.

Lit résultats/phase57_ntm/ntm_presence_per_gene.tsv (produit par phase57) et ajoute `outgroup.ntm` à chaque
fiche (à côté de `outgroup.canettii`). Répond « cœur ancien vs innovation spécifique-MTBC » par présence/absence
robuste (tblastn protéine-vs-génome sur ~53 NTM assemblés). Écrit EN PLACE.
Run: python analyses/phase58_ntm_conservation.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
TSV = ROOT / "résultats" / "phase57_ntm" / "ntm_presence_per_gene.tsv"

INTERP = {
    "genus-core": "conserved across the genus (present in {n}/{tot} non-MTBC Mycobacterium genomes, incl. distant "
                  "relatives) — an ancient core gene predating the genus radiation",
    "restricted": "present in a subset of the genus ({n}/{tot} NTM; in {mtbap} of the 4 closest MTBAP relatives) — "
                  "partial/intermediate conservation",
    "MTBC-near-specific": "absent from the closest MTBAP relatives and nearly all NTM ({n}/{tot}) — a strong "
                          "MTBC-restricted candidate (possible host-adaptation innovation, confirm by synteny)",
    "MTBC-specific": "absent from ALL {tot} non-MTBC Mycobacterium genomes tested (incl. the closest MTBAP relatives) "
                     "— a candidate MTBC-specific genetic innovation / host-adaptation factor (confirm by synteny; "
                     "rule out a short/low-complexity ORF and extreme divergence)",
}

# [2026-07-30] Interprétation SUBSTITUÉE quand phase57 a flaggé des hits sous-seuil : le gène est présent mais
# divergent chez au moins un NTM, donc PAS une innovation MTBC. Sans cela, la phrase « candidate MTBC-specific
# genetic innovation / host-adaptation factor » ci-dessus part telle quelle dans la fiche et fonde des projets
# entiers sur un artefact de seuil (vécu : Rv2438A, homologue M. decipiens 79 % id / 32 % cov manqué par qcov>=50).
INTERP_SUBTHRESHOLD = (
    "NOT MTBC-specific despite passing the strict presence filter: no hit at pident>=30 & qcov>=50 across the "
    "{tot} non-MTBC genomes, BUT sub-threshold tblastn hit(s) exist ({sub}) — the gene is PRESENT BUT DIVERGENT "
    "in at least one non-MTBC Mycobacterium, not a genus-level innovation. Do not frame as MTBC-specific nor as "
    "a host-adaptation factor. Human non-homology, if needed, must be established directly (BLASTp vs human "
    "proteome), never inferred from this field."
)


def main():
    if not TSV.exists():
        raise SystemExit(f"{TSV} introuvable — lancer phase57_ntm_orthology.py d'abord.")
    tbl = {}
    for l in open(TSV):
        if l.startswith("rv"):
            continue
        p = l.rstrip("\n").split("\t")
        tbl[p[0]] = {"n_present": int(p[1]), "n_total": int(p[2]), "frac": float(p[3]),
                     "in_mtbap": int(p[4]), "in_mkc": int(p[5]), "mean_pident": float(p[6]),
                     "classification": p[7],
                     "species_present": [s for s in (p[8].split(",") if len(p) > 8 and p[8] else [])],
                     # colonnes de sensibilité ajoutées le 2026-07-30 ; absentes des TSV antérieurs (rétro-compat).
                     "sensitivity_flag": p[9] if len(p) > 9 else "not-tested",
                     "subthreshold_best": p[10] if len(p) > 10 else "-"}
    n = n_sub = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        rec = tbl.get(d["rv"])
        if not rec:
            continue
        rec = dict(rec)
        flagged = rec["sensitivity_flag"].startswith("SUBTHRESHOLD_HITS")
        if flagged:
            rec["interpretation"] = INTERP_SUBTHRESHOLD.format(tot=rec["n_total"], sub=rec["subthreshold_best"])
            # la classification brute reste lisible, mais ne doit plus être consommée telle quelle en aval.
            rec["classification_superseded"] = True
            rec["corrected_classification"] = "present-but-divergent-in-NTM"
            n_sub += 1
        else:
            rec["interpretation"] = INTERP[rec["classification"]].format(
                n=rec["n_present"], tot=rec["n_total"], mtbap=rec["in_mtbap"])
        rec["source"] = ("presence/absence by tblastn vs ~53 assembled non-MTBC Mycobacterium genomes "
                         "(annotation_mtbc P10.2), STRICT thresholds pident>=30 & qcov>=50, with intrinsic "
                         "sub-threshold sensitivity test (P10.2 revision 2026-07-30)")
        og = d.get("outgroup") or {}
        og["ntm"] = rec
        d["outgroup"] = og
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
    from collections import Counter
    c = Counter(r["classification"] for r in tbl.values())
    print(f"couche `outgroup.ntm` écrite sur {n} fiches. Classification : {dict(c)}")
    if n_sub:
        print(f"  /!\\ {n_sub} fiche(s) « MTBC-specific » RÉTROGRADÉE(S) en 'present-but-divergent-in-NTM' "
              f"(hits sous-seuil détectés par phase57).")
    if any(r["sensitivity_flag"] == "not-tested" for r in tbl.values()):
        print("  /!\\ TSV sans colonnes de sensibilité (antérieur au 2026-07-30) : relancer phase57 pour que "
              "les « MTBC-specific » soient testés contre la divergence.")


if __name__ == "__main__":
    main()
