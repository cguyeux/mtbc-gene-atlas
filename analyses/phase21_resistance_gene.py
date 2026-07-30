#!/usr/bin/env python3
"""phase21_resistance_gene.py -- croisement RÉSISTANCE WHO par gène (P5.9).

Champ de parité Mycobrowser à gain facile : la résistance aux antituberculeux est DÉJÀ
embarquée dans le conteneur (catalogue consolidé WHO 2e éd. 2023 + tb-profiler, servi par
le testeur /resistance). On la SURFACE sur la fiche gène : « ce gène porte N variants
catalogués R-associés pour tel médicament ».

Garde-fou anti-sur-appel (crucial) : le catalogue consolidé mélange des grades WHO 1 à 5
ET des entrées empiriques CRyPTIC (467 « gènes », majoritairement du bruit type GWAS). On
ne retient QUE le tier R-associé WHO = grades « 1) Assoc w R » et « 2) Assoc w R - Interim »
(+ variantes de libellé). Résultat : ~26 gènes canoniques (katG/INH, rpoB/RIF, pncA/PZA,
gyrA/FQ, ethA, gid, embB...), dont 24 mappent sur une fiche protéique (rrs/rrl = ARNr, hors
protéome). Aucun « Uncertain / Not assoc / empirical » n'est surfacé (trompeur sur une fiche).

Mapping gène->Rv : le champ `gene` du catalogue est tantôt un symbole (katG), tantôt un
locus tag (Rv2752c) ; on résout par locus tag direct puis par nom de gène de la fiche.

Sortie : résultats/phase21_resistance/resistance_by_gene.json (keyed Rv). Fusion directe
dans les fiches (idempotent) + enregistrement dans phase4.
Run: python analyses/phase21_resistance_gene.py
"""
from __future__ import annotations
import json, csv, glob
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent
CATALOGUE = ROOT / "site" / "content" / "resistance" / "catalogue_consolide.tsv"
OUT = ROOT / "résultats" / "phase21_resistance"
GENES = ROOT / "site" / "content" / "genes"

SOURCE = "Catalogue consolidé WHO 2e éd. (2023) + tb-profiler (tier R-associé, grades 1–2)"
REFS = [
    {"authors": "World Health Organization", "year": 2023,
     "title": "Catalogue of mutations in Mycobacterium tuberculosis complex and their association with drug resistance, 2nd edition",
     "journal": "WHO", "doi": "10.2665/9789240082410"},
]


def is_R_associated(grade: str | None) -> bool:
    g = (grade or "").strip().lower()
    return g.startswith("1)") or g.startswith("2)") or g.startswith("assoc w r")


def main():
    print("== phase21 : croisement résistance WHO par gène (P5.9) ==")
    # gène (identifiant catalogue) -> drug -> nb de variants R-associés
    gene_drug = defaultdict(lambda: defaultdict(int))
    with open(CATALOGUE, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not is_R_associated(row.get("source_grade")):
                continue
            g = (row.get("gene") or "").strip()
            d = (row.get("drug") or "").strip()
            if g and d:
                gene_drug[g][d] += 1
    print(f"Gènes R-associés (grade 1/2) dans le catalogue : {len(gene_drug)}")

    # index de mapping vers les fiches
    rv_set = set()
    name2rv = {}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        rv_set.add(d["rv"])
        gn = (d.get("gene") or "").strip().lower()
        if gn:
            name2rv.setdefault(gn, d["rv"])

    def to_rv(g):
        if g in rv_set:
            return g
        return name2rv.get(g.lower())

    rec = {}
    unmapped = []
    for g, drugs in gene_drug.items():
        rv = to_rv(g)
        if not rv:
            unmapped.append(g)
            continue
        drug_list = sorted(({"drug": d, "n_variants": n} for d, n in drugs.items()),
                           key=lambda x: -x["n_variants"])
        rec[rv] = {
            "is_resistance_gene": True,
            "catalogue_gene": g,
            "drugs": drug_list,
            "total_variants": sum(drugs.values()),
            "source": SOURCE,
            "refs": REFS,
        }
    print(f"Mappés sur une fiche : {len(rec)} | non mappés (ARNr/non-CDS) : {sorted(unmapped)}")

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "resistance_by_gene.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))
    print(f"Écrit {OUT/'resistance_by_gene.json'} ({len(rec)} gènes de résistance)")

    # ── fusion directe dans les fiches (idempotent) ──
    n_written = n_res = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        d["resistance"] = rec.get(d["rv"], {})
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n_written += 1
        if d["rv"] in rec:
            n_res += 1
    print(f"Fusionné 'resistance' dans {n_written} fiches ({n_res} gènes de résistance surfacés)")
    for rv, r in sorted(rec.items(), key=lambda kv: -kv[1]["total_variants"]):
        drugs = ", ".join(f"{x['drug']}({x['n_variants']})" for x in r["drugs"])
        print(f"  {rv} [{r['catalogue_gene']}] : {drugs}")


if __name__ == "__main__":
    main()
