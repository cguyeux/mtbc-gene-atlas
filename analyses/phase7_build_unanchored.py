#!/usr/bin/env python3
"""phase7_build_unanchored.py -- Catalogue the H37Rv genes with no MTBC0 1:1 anchor.

~920 H37Rv CDS have no 1:1 RefSeq-NP anchor to an MTBC0 gene (PE/PPE/PE-PGRS,
paralogues that share an NP, IS/transposases, and CDS whose PGAP MTBC0 entry
lacks an `inference` NP). They were therefore absent from gene_xref.tsv.

This builds gene_unanchored.tsv in the SAME column layout as gene_xref.tsv, so
the rest of the pipeline (phase2b Pfam, phase5 autocurate, phase4 export) can
consume it via the GENE_XREF env var. Caveat: for these genes the protein is the
H37Rv one (no ancestral MTBC0 sequence available), and coordinates are H37Rv.
mtbc0 is set to '-' to flag the absence of an ancestral anchor.

Output: data/gene_unanchored.tsv  +  data/gene_xref_all.tsv (concat of both)
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parent.parent
H37RV_CDS = ROOT.parent / "investigate_phylo" / "resources" / "NC_000962.3_CDS.fasta"
XREF = ROOT / "data" / "gene_xref.tsv"
OUT = ROOT / "data" / "gene_unanchored.tsv"
OUT_ALL = ROOT / "data" / "gene_xref_all.tsv"

FIELD = re.compile(r"\[(\w+)=([^\]]*)\]")
LOC = re.compile(r"(\d+)\.\.(\d+)")
COLS = ["rv", "mtbc0", "np_id", "gene", "product_h37rv", "product_mtbc0_pgap",
        "is_hypothetical_h37rv", "len_aa", "strand", "start_mtbc0", "end_mtbc0",
        "protein_mtbc0"]


def translate(dna: Seq) -> str:
    n = len(dna) - (len(dna) % 3)
    aa = str(dna[:n].translate(table=11, to_stop=True))
    return ("M" + aa[1:]) if aa and aa[0] != "M" else aa


def main() -> None:
    existing = {r["rv"] for r in csv.DictReader(open(XREF), delimiter="\t")}
    rows = []
    for rec in SeqIO.parse(str(H37RV_CDS), "fasta"):
        f = dict(FIELD.findall(rec.description))
        rv = f.get("locus_tag", "")
        if not rv or rv in existing:
            continue
        loc = f.get("location", "")
        strand = "-" if "complement" in loc else "+"
        m = LOC.search(loc)
        start, end = (m.group(1), m.group(2)) if m else ("", "")
        prot = translate(rec.seq)
        product = f.get("protein", "")
        rows.append({
            "rv": rv, "mtbc0": "-", "np_id": f.get("protein_id", ""),
            "gene": f.get("gene", ""), "product_h37rv": product,
            "product_mtbc0_pgap": "",
            "is_hypothetical_h37rv": str(product.strip().lower() == "hypothetical protein"),
            "len_aa": str(len(prot)), "strand": strand,
            "start_mtbc0": start, "end_mtbc0": end, "protein_mtbc0": prot,
        })
    rows.sort(key=lambda r: r["rv"])
    with open(OUT, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in COLS) + "\n")
    # combined catalogue (anchored + unanchored)
    with open(OUT_ALL, "w") as out:
        out.write(Path(XREF).read_text())
        with open(OUT) as un:
            next(un)  # skip header
            out.write(un.read())
    print(f"Unanchored genes: {len(rows)} -> {OUT.name}")
    print(f"Combined catalogue -> {OUT_ALL.name}")


if __name__ == "__main__":
    main()
