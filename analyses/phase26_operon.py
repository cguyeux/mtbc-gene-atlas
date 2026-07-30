#!/usr/bin/env python3
"""phase26_operon.py -- contexte génomique + opéron prédit (P5.7).

Champ de parité Mycobrowser : le contexte génomique explicite d'un gène (voisins immédiats,
brin, distance intergénique) et son appartenance à une unité de transcription (opéron). Utile
pour interpréter un « hypothetical » : un gène co-transcrit avec un opéron caractérisé hérite
d'une hypothèse de fonction (complément concret du canal neighborhood de STRING).

Source : GFF H37Rv `investigate_phylo/resources/NC_000962.3.gff3` (3978 features `gene`,
tous locus tags). 100 % local, reproductible.

Méthode opéron : prédiction par distance intergénique co-directionnelle (Salgado et al. 2000) :
gènes consécutifs sur le MÊME brin dont tous les gaps internes <= 50 bp = même unité de
transcription prédite. Distribution H37Rv : 28 % des paires same-strand chevauchent, médiane
29 bp -> seuil 50 bp = 59 % des paires (fourchette « co-transcrit probable »). PRÉDICTION,
étiquetée comme telle.

Sortie : résultats/phase26_operon/operon.json (keyed Rv). Fusion dans les fiches (idempotent)
+ enregistrement phase4. Stdlib.
Run: python analyses/phase26_operon.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GFF = ROOT.parent / "investigate_phylo" / "resources" / "NC_000962.3.gff3"
OUT = ROOT / "résultats" / "phase26_operon"
GENES = ROOT / "site" / "content" / "genes"
GAP_MAX = 50

SOURCE = "Genomic context from H37Rv annotation; operon predicted by co-directional intergenic distance (<=50 bp)"
REFS = [
    {"authors": "Salgado H, Moreno-Hagelsieb G, Smith TF, Collado-Vides J", "year": 2000,
     "title": "Operons in Escherichia coli: genomic analyses and predictions",
     "journal": "PNAS", "doi": "10.1073/pnas.030539397"},
]


def load_genes():
    g = []
    for line in GFF.read_text().splitlines():
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 9 or p[2] != "gene":
            continue
        attr = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
        g.append({"start": int(p[3]), "stop": int(p[4]), "strand": p[6],
                  "locus": attr.get("locus_tag", ""), "gene": attr.get("gene", "")})
    g.sort(key=lambda x: x["start"])
    return g


def main():
    print("== phase26 : contexte génomique + opéron prédit (P5.7) ==")
    genes = load_genes()
    print(f"{len(genes)} features gene chargées")

    # opérons = runs maximaux même-brin, gaps internes <= GAP_MAX
    operon_of = {}   # locus -> [membres locus dans l'ordre génomique]
    run = [genes[0]]
    def flush(run):
        loci = [x["locus"] for x in run]
        for x in run:
            operon_of[x["locus"]] = loci
    for i in range(1, len(genes)):
        prev, cur = genes[i - 1], genes[i]
        gap = cur["start"] - prev["stop"] - 1
        if cur["strand"] == prev["strand"] and gap <= GAP_MAX:
            run.append(cur)
        else:
            flush(run); run = [cur]
    flush(run)

    idx = {g["locus"]: i for i, g in enumerate(genes)}

    def ctx(nb, gap):
        return {"locus": nb["locus"], "gene": nb["gene"] or None, "strand": nb["strand"], "gap": gap}

    rec = {}
    for g in genes:
        i = idx[g["locus"]]
        left = genes[i - 1] if i > 0 else None
        right = genes[i + 1] if i < len(genes) - 1 else None
        left_ctx = ctx(left, g["start"] - left["stop"] - 1) if left else None
        right_ctx = ctx(right, right["start"] - g["stop"] - 1) if right else None
        members = [m for m in operon_of.get(g["locus"], [g["locus"]])]
        member_objs = [{"locus": m, "gene": genes[idx[m]]["gene"] or None} for m in members]
        rec[g["locus"]] = {
            "strand": g["strand"],
            "left": left_ctx,
            "right": right_ctx,
            "operon": member_objs,
            "operon_size": len(member_objs),
            "source": SOURCE, "refs": REFS,
        }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "operon.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))
    print(f"Écrit {OUT/'operon.json'} ({len(rec)} loci)")

    # fusion dans les fiches
    n_written = n_hit = n_multi = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        r = rec.get(d["rv"])
        d["genomic_context"] = r or {}
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n_written += 1
        if r:
            n_hit += 1
            if r["operon_size"] >= 2:
                n_multi += 1
    print(f"Fusionné 'genomic_context' dans {n_written} fiches ({n_hit} avec contexte, "
          f"{n_multi} dans un opéron multi-gènes)")

    # bilan opérons
    from collections import Counter
    sizes = Counter(len(v) for v in {tuple(o): o for o in operon_of.values()}.values())
    n_operons = len({tuple(o) for o in operon_of.values()})
    print(f"Unités de transcription prédites : {n_operons} ; tailles {dict(sorted(sizes.items()))}")


if __name__ == "__main__":
    main()
