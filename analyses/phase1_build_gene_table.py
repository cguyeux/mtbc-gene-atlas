#!/usr/bin/env python3
"""phase1_build_gene_table.py -- Reference gene table for the Mycobrowser successor.

Builds the cross-reference between the MTBC ancestral genome MTBC0
(Harrison et al. 2024) and H37Rv (NC_000962.3), and extracts the
*ancestral* protein sequence of every gene from MTBC0.

The project starts from MTBC0: the protein we annotate is the ancestral
MTBC0 CDS, not the H37Rv one. H37Rv stays as the historical anchor (its
locus tag Rv#### and its legacy product, including the 'hypothetical
protein' label we aim to resolve).

Join key: the MTBC0 PGAP annotation carries, for each CDS,
`inference=COORDINATES: similar to AA sequence:RefSeq:NP_######.1`, which
is exactly the H37Rv RefSeq protein_id found in the H37Rv CDS catalogue
header. We join on that NP_ id -- deterministic, offline, no scraping.

Inputs (all local):
  data/MTBC0/MTBC0_v1.1.fasta              ancestral genome
  data/MTBC0/MTBC0v1.1_PGAP_annot.gff      MTBC0 PGAP gene annotation
  ../investigate_phylo/resources/NC_000962.3_CDS.fasta   H37Rv CDS catalogue

Output:
  data/gene_xref.tsv      one row per joined gene
  résultats/phase1_summary.txt
"""
from __future__ import annotations

import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parent.parent
MTBC0_FASTA = ROOT / "data" / "MTBC0" / "MTBC0_v1.1.fasta"
MTBC0_GFF = ROOT / "data" / "MTBC0" / "MTBC0v1.1_PGAP_annot.gff"
H37RV_CDS = ROOT.parent / "investigate_phylo" / "resources" / "NC_000962.3_CDS.fasta"
OUT_TSV = ROOT / "data" / "gene_xref.tsv"
OUT_SUMMARY = ROOT / "résultats" / "phase1_summary.txt"

# H37Rv CDS header, e.g.:
#  >lcl|NC_000962.3_cds_NP_214518.1_4 [locus_tag=Rv0004] [protein=hypothetical protein] [protein_id=NP_214518.1] ...
H37RV_FIELD = re.compile(r"\[(\w+)=([^\]]*)\]")
NP_FROM_INFERENCE = re.compile(r"RefSeq:(NP_\d+\.\d+)")


def parse_h37rv() -> dict[str, dict]:
    """NP protein_id -> {rv, gene, product}."""
    by_np: dict[str, dict] = {}
    for rec in SeqIO.parse(str(H37RV_CDS), "fasta"):
        fields = dict(H37RV_FIELD.findall(rec.description))
        np_id = fields.get("protein_id")
        if not np_id:
            continue
        by_np[np_id] = {
            "rv": fields.get("locus_tag", ""),
            "gene": fields.get("gene", ""),
            "product_h37rv": fields.get("protein", ""),
        }
    return by_np


def parse_gff_attrs(col9: str) -> dict[str, str]:
    out = {}
    for kv in col9.strip().rstrip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def parse_mtbc0() -> list[dict]:
    """One dict per MTBC0 CDS, with coords, strand and the H37Rv NP id."""
    cds = []
    with open(MTBC0_GFF) as fh:
        for line in fh:
            if line.startswith("#") or "\tCDS\t" not in line:
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            attrs = parse_gff_attrs(c[8])
            np_match = NP_FROM_INFERENCE.search(attrs.get("inference", ""))
            cds.append({
                "mtbc0": attrs.get("locus_tag", ""),
                "gene_pgap": attrs.get("gene", ""),
                "product_mtbc0_pgap": attrs.get("product", ""),
                "np_id": np_match.group(1) if np_match else "",
                "start": int(c[3]),
                "end": int(c[4]),
                "strand": c[6],
            })
    return cds


def translate_cds(genome: Seq, start: int, end: int, strand: str) -> str:
    sub = genome[start - 1:end]
    if strand == "-":
        sub = sub.reverse_complement()
    # Genetic table 11; force alternative initiators (GTG/TTG/CTG) to M.
    aa = str(sub.translate(table=11, to_stop=True))
    if aa and aa[0] != "M":
        aa = "M" + aa[1:]
    return aa


def main() -> None:
    by_np = parse_h37rv()
    mtbc0_cds = parse_mtbc0()
    genome = next(SeqIO.parse(str(MTBC0_FASTA), "fasta")).seq

    rows = []
    for cds in mtbc0_cds:
        h = by_np.get(cds["np_id"])
        if not h:
            continue  # no 1:1 H37Rv RefSeq anchor -> skip for the cross-ref table
        prot = translate_cds(genome, cds["start"], cds["end"], cds["strand"])
        rows.append({
            "rv": h["rv"],
            "mtbc0": cds["mtbc0"],
            "np_id": cds["np_id"],
            "gene": h["gene"] or cds["gene_pgap"],
            "product_h37rv": h["product_h37rv"],
            "product_mtbc0_pgap": cds["product_mtbc0_pgap"],
            "is_hypothetical_h37rv": str(h["product_h37rv"].strip().lower() == "hypothetical protein"),
            "len_aa": str(len(prot)),
            "strand": cds["strand"],
            "start_mtbc0": str(cds["start"]),
            "end_mtbc0": str(cds["end"]),
            "protein_mtbc0": prot,
        })

    rows.sort(key=lambda r: r["rv"])
    cols = ["rv", "mtbc0", "np_id", "gene", "product_h37rv", "product_mtbc0_pgap",
            "is_hypothetical_h37rv", "len_aa", "strand", "start_mtbc0", "end_mtbc0",
            "protein_mtbc0"]
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_TSV, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in cols) + "\n")

    n_hyp = sum(r["is_hypothetical_h37rv"] == "True" for r in rows)
    # How many H37Rv 'hypothetical' get a *non*-hypothetical PGAP product on MTBC0?
    rescued = [r for r in rows
               if r["is_hypothetical_h37rv"] == "True"
               and r["product_mtbc0_pgap"].strip().lower() not in ("hypothetical protein", "")]
    lines = [
        f"Joined genes (MTBC0 <-> H37Rv via RefSeq NP id): {len(rows)}",
        f"H37Rv 'hypothetical protein': {n_hyp}",
        f"  of which PGAP/MTBC0 already gives a non-hypothetical product: {len(rescued)}",
        "",
        "First 5 H37Rv-hypothetical genes (pilot batch candidates):",
    ]
    pilot = [r for r in rows if r["is_hypothetical_h37rv"] == "True"][:5]
    for r in pilot:
        lines.append(f"  {r['rv']:9s} {r['mtbc0']:14s} len={r['len_aa']:>4} aa  "
                     f"PGAP='{r['product_mtbc0_pgap']}'")
    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    OUT_SUMMARY.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote {OUT_TSV}")


if __name__ == "__main__":
    main()
