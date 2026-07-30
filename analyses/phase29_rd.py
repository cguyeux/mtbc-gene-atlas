#!/usr/bin/env python3
"""phase29_rd.py -- régions de différence (RD) / délétions par gène (P5.10).

Dernier champ de parité Mycobrowser (contexte évolutif) : ce gène tombe-t-il dans une RD
connue, et dans quelles lignées est-il DÉLÉTÉ ? Une RD est une grande délétion, souvent
synapomorphie de lignée (RD9 des lignées animales, RD1 de BCG, TbD1 des lignées modernes...).
Pour un « hypothetical » délété dans une lignée, l'absence dispensable est un signal ; pour
tout gène, la carte de délétion est un contexte évolutif/diagnostique précieux.

Sources (natives du projet, analyse RD consolidée sur ~145k souches, bien plus riche que la
table RD historique) :
  - `investigate_phylo/resources/regions/rd.bed` : 223 RD canoniques, coordonnées H37Rv (BED).
  - `global_supplementary/RD_canonical_frequencies_2026-05-24.tsv` : fréquence d'absence par
    lignée (f_L1..f_L8, f_Bovis, f_Orygis, f_Caprae, animales...).
  - `global_supplementary/RD_canonical_branch_mapping_2026-05-24.tsv` : branche d'origine.
Gènes en coordonnées H37Rv via le GFF (NC_000962.3.gff3). Overlap = gène délété (tout ou partie)
dans les lignées où la RD est absente.

Sortie : résultats/phase29_rd/rd.json (keyed Rv). Fusion dans les fiches (idempotent) + phase4.
Run: python analyses/phase29_rd.py
"""
from __future__ import annotations
import json, csv, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GFF = ROOT.parent / "investigate_phylo" / "resources" / "NC_000962.3.gff3"
BED = ROOT.parent / "investigate_phylo" / "resources" / "regions" / "rd.bed"
FREQ = ROOT.parent / "global_supplementary" / "RD_canonical_frequencies_2026-05-24.tsv"
BRANCH = ROOT.parent / "global_supplementary" / "RD_canonical_branch_mapping_2026-05-24.tsv"
OUT = ROOT / "résultats" / "phase29_rd"
GENES = ROOT / "site" / "content" / "genes"

F_ABSENT = 0.5    # fréquence d'absence pour dire « délété dans la lignée »
MIN_COVERED = 5   # nb minimal de souches couvertes pour retenir une lignée

SOURCE = "Consolidated MTBC Regions-of-Difference analysis (145k strains; H37Rv coordinates)"
REFS = [
    {"authors": "Brosch R, Gordon SV, Marmiesse M, et al.", "year": 2002,
     "title": "A new evolutionary scenario for the Mycobacterium tuberculosis complex",
     "journal": "PNAS", "doi": "10.1073/pnas.052548299"},
    {"authors": "Gagneux S, Small PM", "year": 2007,
     "title": "Global phylogeography of Mycobacterium tuberculosis and implications for tuberculosis product development",
     "journal": "The Lancet Infectious Diseases", "doi": "10.1016/S1473-3099(07)70108-1"},
]


def load_beds():
    rds = []
    for line in BED.read_text().splitlines():
        p = line.split("\t")
        if len(p) < 4:
            continue
        rds.append({"name": p[3], "start": int(p[1]), "end": int(p[2])})
    return rds


def load_freqs():
    """rd_name -> [(lineage, f_absence)] pour f>=F_ABSENT et couverture suffisante."""
    out = {}
    with open(FREQ, encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        lin_cols = [c[2:] for c in (r.fieldnames or []) if c.startswith("f_")]
        for row in r:
            deleted = []
            for lin in lin_cols:
                try:
                    f = float(row.get(f"f_{lin}") or 0)
                    cov = float(row.get(f"covered_{lin}") or 0)
                except ValueError:
                    continue
                if f >= F_ABSENT and cov >= MIN_COVERED:
                    deleted.append({"lineage": lin, "freq": round(f, 2)})
            deleted.sort(key=lambda x: -x["freq"])
            out[row["rd_name"]] = deleted
    return out


def load_branch():
    out = {}
    if BRANCH.exists():
        for row in csv.DictReader(open(BRANCH, encoding="utf-8"), delimiter="\t"):
            out[row["rd_name"]] = (row.get("branch_of_origin") or "").strip()
    return out


def load_genes():
    g = []
    for line in GFF.read_text().splitlines():
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 9 or p[2] != "gene":
            continue
        attr = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
        g.append({"locus": attr.get("locus_tag", ""), "start": int(p[3]), "end": int(p[4])})
    return g


def main():
    print("== phase29 : régions de différence / délétions par gène (P5.10) ==")
    rds = load_beds()
    freqs = load_freqs()
    branch = load_branch()
    genes = load_genes()
    print(f"{len(rds)} RD (bed) | {len(freqs)} RD avec fréquences | {len(genes)} gènes")

    # overlap gène x RD (balayage simple ; 223 RD × 3978 gènes = ok)
    rec = {}
    for g in genes:
        hits = []
        for rd in rds:
            if g["start"] <= rd["end"] and g["end"] >= rd["start"]:
                deleted = freqs.get(rd["name"], [])
                # ne retenir la RD que si au moins une lignée est réellement délétée
                if not deleted:
                    continue
                # portion du gène couverte par la RD
                ov = min(g["end"], rd["end"]) - max(g["start"], rd["start"]) + 1
                frac = round(ov / (g["end"] - g["start"] + 1), 2)
                hits.append({
                    "rd": rd["name"],
                    "branch_of_origin": branch.get(rd["name"], ""),
                    "overlap_frac": frac,
                    "deleted_in": deleted[:12],
                })
        if hits:
            hits.sort(key=lambda h: -h["overlap_frac"])
            rec[g["locus"]] = {"rds": hits, "source": SOURCE, "refs": REFS}

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "rd.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))
    print(f"Écrit {OUT/'rd.json'} ({len(rec)} gènes dans >=1 RD délétée)")

    # fusion dans les fiches
    n_written = n_hit = hyp_hit = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        d["rd"] = rec.get(d["rv"], {})
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n_written += 1
        if d["rv"] in rec:
            n_hit += 1
            if "hypothetical" in (d.get("product_h37rv") or "").lower():
                hyp_hit += 1
    print(f"Fusionné 'rd' dans {n_written} fiches ({n_hit} dans une RD délétée, dont {hyp_hit} hypothétiques)")

    # bilan : lignées les plus concernées
    from collections import Counter
    linc = Counter()
    for r in rec.values():
        seen = set()
        for h in r["rds"]:
            for dl in h["deleted_in"]:
                seen.add(dl["lineage"])
        for lin in seen:
            linc[lin] += 1
    print("Gènes délétés par lignée (top) :", linc.most_common(10))


if __name__ == "__main__":
    main()
