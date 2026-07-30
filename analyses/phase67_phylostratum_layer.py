#!/usr/bin/env python3
"""phase67_phylostratum_layer.py -- P10.3 : écrit la couche `outgroup.phylostratum` (âge du gène).

Combine la profondeur CROSS-GENRE (phase66) avec la présence DANS le genre déjà écrite en `outgroup.ntm`
(phase58) pour poser, sur chaque fiche, un phylostrate ordinal = rang taxonomique le plus profond où le
gène est encore détecté :

    0 MTBC-specific  <  1 Mycobacterium  <  2 Mycobacteriaceae  <  3 Corynebacteriales
      <  4 Actinomycetia  <  5 Bacteria (universel)

Écrit EN PLACE dans site/content/genes/*.json, fusionné dans le dict `outgroup` (à côté de canettii + ntm).

GARDE-FOU null-first (KB [2026-07-11], Moyers & Zhang) : un phylostrate SUPERFICIEL peut être un artefact
d'échec de détection d'homologie (les gènes courts détectent moins loin). Le script (a) porte un flag
`short_gene_caveat` (len<100 aa) sur chaque fiche superficielle, (b) imprime la distribution de longueur
PAR strate pour que la lecture vérifie que « jeune » n'est pas juste « court ». On ne revendique JAMAIS une
nouveauté biologique sur la seule base d'un phylostrate superficiel : c'est un axe descriptif, pas une preuve.

Run: python analyses/phase67_phylostratum_layer.py    (après phase66)
"""
from __future__ import annotations
import json, glob, statistics
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
TSV = ROOT / "résultats" / "phase66_phylostratum" / "crossgenus_presence_per_gene.tsv"

STRATUM_LABEL = {0: "MTBC-specific", 1: "Mycobacterium", 2: "Mycobacteriaceae",
                 3: "Corynebacteriales", 4: "Actinomycetia", 5: "Bacteria"}
# rang cross-genre (phase66) -> strate ordinale
XG_TO_STRATUM = {1: 2, 2: 3, 3: 4, 4: 5}

INTERP = {
    5: "detected down to outside the phylum (Proteobacteria/Firmicutes controls) — a universally conserved, "
       "ancient bacterial gene",
    4: "detected across the class Actinomycetia (beyond Corynebacteriales) but not outside the phylum — an "
       "Actinobacteria-level ancient gene",
    3: "detected across the order Corynebacteriales (Corynebacterium/Nocardia/Rhodococcus/…) but not in more "
       "distant Actinomycetia — a Corynebacteriales-level gene",
    2: "detected in the sister family Mycobacteriaceae (M. abscessus) but not in the broader Corynebacteriales — "
       "a Mycobacteriaceae-restricted gene (single-genome rung: interpret with care)",
    1: "present across the genus Mycobacterium (NTM) but not detected in any non-Mycobacterium genome — a "
       "Mycobacterium-genus gene",
    0: "absent from the whole genus and from all outgroups tested — a candidate MTBC-specific innovation "
       "(confirm by synteny; rule out detection failure for short/divergent ORFs)",
}


def genus_present(d: dict) -> bool:
    """le gène est-il détecté dans >=1 génome NTM (couche outgroup.ntm de phase58) ?"""
    ntm = ((d.get("outgroup") or {}).get("ntm")) or {}
    return bool(ntm.get("n_present"))


def main():
    if not TSV.exists():
        raise SystemExit(f"{TSV} introuvable — lancer phase66_phylostratum_scan.py d'abord.")
    xg = {}
    for l in open(TSV):
        if l.startswith("rv\t"):
            continue
        p = l.rstrip("\n").split("\t")
        xg[p[0]] = {"len_aa": int(p[1]), "n_hits": int(p[2]), "deepest_rank": int(p[3]),
                    "deepest_label": p[4], "mean_pident": float(p[9]),
                    "genomes_hit": [g for g in (p[10].split(",") if len(p) > 10 and p[10] else [])]}

    n = 0
    strat_counts = Counter()
    len_by_stratum = defaultdict(list)
    short_by_stratum = Counter()
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        rv = d["rv"]
        rec = xg.get(rv)
        if rec is None:
            continue
        deepest = rec["deepest_rank"]
        if deepest >= 1:
            stratum = XG_TO_STRATUM[deepest]
        else:                              # aucun hit cross-genre : le genre tranche
            stratum = 1 if genus_present(d) else 0
        length = rec["len_aa"]
        short = length < 100
        layer = {
            "stratum": stratum,
            "stratum_label": STRATUM_LABEL[stratum],
            "deepest_outgroup_clade": rec["deepest_label"] if deepest else ("genus Mycobacterium" if stratum == 1 else "none"),
            "n_crossgenus_genomes": rec["n_hits"],
            "crossgenus_genomes": rec["genomes_hit"],
            "mean_pident": rec["mean_pident"],
            "length_aa": length,
            "short_gene_caveat": bool(short and stratum <= 1),   # risque d'échec de détection surtout aux strates jeunes
            "interpretation": INTERP[stratum],
            "source": "cross-genus phylostratigraphy: tblastn vs 13 non-Mycobacterium reference genomes "
                      "(Mycobacteriaceae>Corynebacteriales>Actinomycetia>outside-phylum) + intra-genus NTM "
                      "presence (annotation_mtbc P10.3)",
        }
        og = d.get("outgroup") or {}
        og["phylostratum"] = layer
        d["outgroup"] = og
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
        strat_counts[stratum] += 1
        len_by_stratum[stratum].append(length)
        if short:
            short_by_stratum[stratum] += 1

    print(f"couche `outgroup.phylostratum` écrite sur {n} fiches.\n")
    print("Distribution par strate (0=MTBC-spécifique … 5=Bacteria) :")
    for s in range(6):
        lens = len_by_stratum.get(s, [])
        med = int(statistics.median(lens)) if lens else 0
        pct_short = (100 * short_by_stratum.get(s, 0) / len(lens)) if lens else 0
        print(f"  {s} {STRATUM_LABEL[s]:<16} n={strat_counts.get(s,0):<5} médiane {med:>4} aa  "
              f"| <100aa : {pct_short:.0f}%")
    print("\nCONTRÔLE MODÈLE NUL (Moyers & Zhang) : si les strates JEUNES (0-1) sont dominées par des gènes")
    print("courts alors que les strates ANCIENNES (4-5) sont longues, une part du signal = échec de détection,")
    print("pas de la vraie jeunesse. Lire le gradient de médiane ci-dessus AVANT toute interprétation biologique.")


if __name__ == "__main__":
    main()
