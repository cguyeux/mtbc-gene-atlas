#!/usr/bin/env python3
"""phase2f_conservation.py -- Intra-MTBC selection & pseudogene signal per gene.

This is the layer no frozen database (Mycobrowser, UniProt) can have: it reads
the SPDI variants of the whole local strain collection (bdd/actuelle, ~146k
strains) and derives, per gene:
  - pN/pS  : nonsynonymous vs synonymous polymorphism, normalised by the number
             of possible sites -> a purifying-selection signal. A 'hypothetical'
             under strong purifying selection (low pN/pS) is very likely a real,
             functional gene; one drifting near pN/pS ~ 1 may be a pseudogene.
  - disruption load : how many strains carry a premature stop (nonsense SNP) or a
             frameshifting indel inside the CDS -> a pseudogenisation signal
             (lineage-specific gene loss shows up as a disrupting allele shared by
             many strains).
  - SNP density.

Two passes (the heavy one is cached):
  Pass 1  walk every spdi.txt, accumulate allele counts per position. Writes
          résultats/phase2f_conservation/{snp_counts.tsv, indel_counts.tsv, meta.json}.
          Reusable for any downstream variant analysis.
  Pass 2  map each variable position to its codon (auto-detecting the SPDI
          coordinate convention against the genome), classify syn/missense/
          nonsense/frameshift, aggregate per gene -> conservation.json (keyed by
          rv, consumed by phase4).

Inputs:  ../../bdd/actuelle/*/*/NC_000962.3/spdi.txt ; H37Rv genome + CDS fasta.
Usage:   python phase2f_conservation.py                 # both passes
         CONS_FORCE=1 python phase2f_conservation.py     # recompute pass 1
"""
from __future__ import annotations

import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent  # .../mtbc/annotation_mtbc
MTBC = ROOT.parent                             # .../mtbc
BDD = Path(os.environ.get("TBANNOTATOR_BDD", MTBC / "bdd")) / "actuelle"
RES = ROOT.parent / "investigate_phylo" / "resources"
GENOME_FA = RES / "NC_000962.3.fasta"
CDS_FA = RES / "NC_000962.3_CDS.fasta"
OUTDIR = ROOT / "résultats" / "phase2f_conservation"
CHROM = "NC_000962.3"

# Disruption-frequency threshold to flag a pseudogene candidate (fraction of strains
# sharing the most common disrupting allele). Configurable.
PSEUDO_FRAC = float(os.environ.get("CONS_PSEUDO_FRAC", "0.01"))
# Minimum allele FREQUENCY for a SNP site to enter the pN/pS estimate. At ~146k
# strains a fixed count is meaningless (count>=3 is a 0.002% frequency, admitting
# nearly every rare deleterious variant and biasing pN/pS toward the neutral
# mutational spectrum). A frequency floor keeps only variants that have actually
# RISEN despite purifying selection -> these reveal the real constraint. The
# effective threshold is max(MIN_COUNT, MIN_FRAC * n_strains).
MIN_COUNT = int(os.environ.get("CONS_MIN_COUNT", "3"))
MIN_FRAC = float(os.environ.get("CONS_MIN_FRAC", "0.001"))

CODON_TABLE = {  # standard genetic code; '*' = stop
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
COMP = str.maketrans("ACGT", "TGCA")


# ---------------------------------------------------------------- pass 1
def pass1_count() -> dict:
    snp = defaultdict(int)        # (pos:int, alt:str) -> count   (pos as in SPDI)
    indel = defaultdict(int)      # (pos:int, ref:str, alt:str) -> count
    n_strains = 0
    OUTDIR.mkdir(parents=True, exist_ok=True)
    limit = int(os.environ.get("CONS_LIMIT", "0"))  # 0 = all
    paths = BDD.glob(f"*/*/{CHROM}/spdi.txt")
    for i, p in enumerate(paths):
        if limit and i >= limit:
            break
        try:
            data = p.read_text()
        except OSError:
            continue
        n_strains += 1
        for line in data.split("\n"):
            if not line:
                continue
            # CHROM:pos:ref:alt
            parts = line.split(":")
            if len(parts) != 4:
                continue
            pos, ref, alt = parts[1], parts[2], parts[3]
            if len(ref) == 1 and len(alt) == 1:
                snp[(int(pos), alt)] += 1
            else:
                indel[(int(pos), ref, alt)] += 1
        if i % 20000 == 0 and i:
            print(f"  ...{i} strains, {len(snp)} SNP keys", flush=True)

    with open(OUTDIR / "snp_counts.tsv", "w") as fh:
        fh.write("pos\talt\tcount\n")
        for (pos, alt), c in snp.items():
            fh.write(f"{pos}\t{alt}\t{c}\n")
    with open(OUTDIR / "indel_counts.tsv", "w") as fh:
        fh.write("pos\tref\talt\tcount\n")
        for (pos, ref, alt), c in indel.items():
            fh.write(f"{pos}\t{ref}\t{alt}\t{c}\n")
    meta = {"n_strains": n_strains, "n_snp_keys": len(snp), "n_indel_keys": len(indel)}
    (OUTDIR / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[pass1] {n_strains} strains, {len(snp)} SNP keys, {len(indel)} indel keys")
    return meta


# ---------------------------------------------------------------- helpers
def load_genome() -> str:
    seq = []
    for line in GENOME_FA.read_text().split("\n"):
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq)


_LOC = re.compile(r"\[location=(complement\()?(\d+)\.\.(\d+)")
_TAG = re.compile(r"\[locus_tag=([^\]]+)\]")


def load_cds() -> list[dict]:
    """Parse CDS fasta -> list of {rv, start, end (1-based incl), strand, seq}."""
    genes = []
    cur_hdr, cur_seq = None, []
    def flush():
        if cur_hdr is None:
            return
        m = _LOC.search(cur_hdr); t = _TAG.search(cur_hdr)
        if not m or not t:
            return
        strand = "-" if m.group(1) else "+"
        genes.append({"rv": t.group(1), "start": int(m.group(2)),
                      "end": int(m.group(3)), "strand": strand,
                      "seq": "".join(cur_seq)})
    for line in CDS_FA.read_text().split("\n"):
        if line.startswith(">"):
            flush(); cur_hdr, cur_seq = line, []
        else:
            cur_seq.append(line.strip())
    flush()
    return genes


def count_sites(seq: str) -> tuple[float, float]:
    """Nei-Gojobori possible synonymous / nonsynonymous sites for a CDS."""
    syn = nsyn = 0.0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon)
        if aa is None:
            continue
        for j in range(3):
            s = 0
            for b in "ACGT":
                if b == codon[j]:
                    continue
                mut = codon[:j] + b + codon[j+1:]
                if CODON_TABLE.get(mut) == aa:
                    s += 1
            syn += s / 3.0
            nsyn += (3 - s) / 3.0
    return syn, nsyn


# ---------------------------------------------------------------- pass 2
def detect_offset(genome: str, snp_rows: list[tuple[int, str, int]]) -> int:
    """Return the offset (0 or -1 added to SPDI pos) giving 1-based genome match.

    We test whether genome[pos-1] (SPDI 1-based) or genome[pos] (SPDI 0-based)
    equals the SPDI ref on a sample. The SPDI line does not store ref for SNPs in
    our compact table, so we recover ref from the genome and just require internal
    consistency: under the right convention, the *alt* differs from the genome
    base; under the wrong one, alt would frequently equal the genome base.
    """
    glen = len(genome)
    # Use high-count SNPs (lineage markers) as the probe set.
    probes = sorted(snp_rows, key=lambda r: -r[2])[:5000]
    best, best_conv = -1, 0
    for conv, shift in (("1-based", 1), ("0-based", 0)):
        ok = 0
        for pos, alt, _ in probes:
            gi = pos - shift
            if 0 <= gi < glen and genome[gi] != alt:  # alt should differ from ref
                ok += 1
        if ok > best:
            best, best_conv = ok, shift
    print(f"[detect] coordinate shift = {best_conv} "
          f"({best}/{len(probes)} probes have alt != genome base)")
    return best_conv


def pass2_analyse() -> dict:
    genome = load_genome()
    glen = len(genome)
    genes = load_cds()
    meta = json.loads((OUTDIR / "meta.json").read_text())
    n_strains = meta["n_strains"]

    # pos (1-based genome) -> gene index
    pos2gene = [-1] * (glen + 2)
    for gi, g in enumerate(genes):
        for p in range(g["start"], g["end"] + 1):
            if 0 < p <= glen:
                pos2gene[p] = gi

    snp_rows = []
    for line in (OUTDIR / "snp_counts.tsv").read_text().split("\n")[1:]:
        if not line:
            continue
        pos, alt, c = line.split("\t")
        snp_rows.append((int(pos), alt, int(c)))
    shift = detect_offset(genome, snp_rows)
    thr = max(MIN_COUNT, int(MIN_FRAC * n_strains))
    print(f"[pass2] allele-frequency floor = {thr} strains "
          f"(max of count {MIN_COUNT}, {MIN_FRAC:.1%} of {n_strains})")

    agg = defaultdict(lambda: {"syn": 0, "missense": 0, "nonsense": 0,
                               "syn_occ": 0, "mis_occ": 0, "non_occ": 0,
                               "fs": 0, "max_disrupt": 0, "n_snp_sites": 0,
                               "disrupt_sites": 0})

    for pos, alt, c in snp_rows:
        if c < thr:               # keep only variants above the frequency floor
            continue
        gpos = pos - shift + 1   # -> 1-based genome coordinate
        if not (0 < gpos <= glen):
            continue
        gi = pos2gene[gpos]
        if gi < 0:
            continue
        g = genes[gi]
        # offset within the (oriented) CDS
        if g["strand"] == "+":
            off = gpos - g["start"]
            cds_alt = alt
        else:
            off = g["end"] - gpos
            cds_alt = alt.translate(COMP)
        seq = g["seq"]
        if not (0 <= off < len(seq)):
            continue
        ci = (off // 3) * 3
        codon = seq[ci:ci+3]
        if len(codon) != 3:
            continue
        j = off % 3
        ref_aa = CODON_TABLE.get(codon)
        mut_codon = codon[:j] + cds_alt + codon[j+1:]
        alt_aa = CODON_TABLE.get(mut_codon)
        if ref_aa is None or alt_aa is None:
            continue
        a = agg[g["rv"]]
        a["n_snp_sites"] += 1
        if alt_aa == ref_aa:
            a["syn"] += 1; a["syn_occ"] += c
        elif alt_aa == "*":
            a["nonsense"] += 1; a["non_occ"] += c
            a["disrupt_sites"] += 1
            a["max_disrupt"] = max(a["max_disrupt"], c)
        else:
            a["missense"] += 1; a["mis_occ"] += c

    # frameshifting indels
    for line in (OUTDIR / "indel_counts.tsv").read_text().split("\n")[1:]:
        if not line:
            continue
        pos, ref, alt, c = line.split("\t")
        pos, c = int(pos), int(c)
        if c < thr:
            continue
        gpos = pos - shift + 1
        if not (0 < gpos <= glen):
            continue
        gi = pos2gene[gpos]
        if gi < 0:
            continue
        if (len(alt) - len(ref)) % 3 != 0:
            a = agg[genes[gi]["rv"]]
            a["fs"] += 1
            a["disrupt_sites"] += 1
            a["max_disrupt"] = max(a["max_disrupt"], c)

    # possible-site normalisation + verdicts
    sites = {g["rv"]: count_sites(g["seq"]) for g in genes}
    out = {}
    for rv, a in agg.items():
        syn_s, nsyn_s = sites.get(rv, (0.0, 0.0))
        pS = a["syn"] / syn_s if syn_s else 0.0
        pN = (a["missense"] + a["nonsense"]) / nsyn_s if nsyn_s else 0.0
        pnps = (pN / pS) if pS > 0 else None
        disrupt_frac = a["max_disrupt"] / n_strains if n_strains else 0.0
        if pnps is None:
            sel = "n/a"
        elif pnps < 0.2:
            sel = "strong purifying"
        elif pnps < 0.5:
            sel = "purifying"
        elif pnps <= 1.2:
            sel = "relaxed/neutral"
        else:
            sel = "diversifying/relaxed"
        # Disruption interpretation. A disrupting allele in >50% of strains means
        # H37Rv (the reference) is the odd one out: it is H37Rv-specific, NOT a
        # sublineage loss (reference bias). A true clade pseudogene sits at an
        # intermediate frequency. Many distinct disrupting alleles = convergent
        # loss-of-function (typical of drug-resistance genes).
        disrupt_mode = ""
        if a["disrupt_sites"]:
            if disrupt_frac >= 0.5:
                disrupt_mode = "reference-fixed"
            elif a["disrupt_sites"] >= 4:
                disrupt_mode = "convergent"
            else:
                disrupt_mode = "clonal"
        out[rv] = {
            "n_strains": n_strains,
            "snp_sites": a["n_snp_sites"],
            "syn": a["syn"], "missense": a["missense"], "nonsense": a["nonsense"],
            "frameshift": a["fs"],
            "pN_pS": round(pnps, 3) if pnps is not None else None,
            "selection": sel,
            "disrupt_sites": a["disrupt_sites"],
            "disrupt_mode": disrupt_mode,
            "disrupt_strains": a["max_disrupt"],
            "disrupt_frac": round(disrupt_frac, 4),
            "pseudogene_flag": bool(PSEUDO_FRAC <= disrupt_frac < 0.5 and a["disrupt_sites"] < 4),
        }
    (OUTDIR / "conservation.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))

    n_pseudo = sum(1 for v in out.values() if v["pseudogene_flag"])
    n_strong = sum(1 for v in out.values() if v["selection"] == "strong purifying")
    print(f"[pass2] {len(out)} genes scored | {n_strong} strong-purifying | "
          f"{n_pseudo} pseudogene candidates (disrupt >= {PSEUDO_FRAC:.0%})")
    print(f"Wrote {OUTDIR / 'conservation.json'}")
    return out


def main() -> None:
    have = (OUTDIR / "snp_counts.tsv").exists() and (OUTDIR / "meta.json").exists()
    if not have or os.environ.get("CONS_FORCE"):
        pass1_count()
    else:
        print("[pass1] cached (snp_counts.tsv present); set CONS_FORCE=1 to recompute")
    pass2_analyse()


if __name__ == "__main__":
    main()
