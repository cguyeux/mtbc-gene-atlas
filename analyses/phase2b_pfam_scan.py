#!/usr/bin/env python3
"""phase2b_pfam_scan.py -- Pfam domain assignment for MTBC0 ancestral proteins.

Runs `hmmscan --cut_ga` (HMMER) against the local Pfam-A database on the
ancestral MTBC0 protein of each gene in a batch. This is the *traceable*
domain layer the project requires: a Pfam accession (PFxxxxx) with an
E-value, as opposed to the exploratory ESM SAE labels. Even genes left
'dark' after PGAP + literature get a confirmed (or newly discovered)
domain assignment here, with a citable accession.

This is the right tool precisely for the 'still unknown' genes: a DUF
(Domain of Unknown Function) is itself a Pfam family, so hmmscan confirms
it with an accession; and a more sensitive scan sometimes reveals an
additional domain PGAP did not list.

Inputs:  data/gene_xref.tsv, a batch file (list of Rv)
Output:  résultats/phase2b_pfam/<rv>.domtblout (raw) and pfam.json (parsed)

Usage:  python phase2b_pfam_scan.py data/pilot_batch_01.txt
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
XREF = Path(os.environ.get("GENE_XREF", ROOT / "data" / "gene_xref.tsv"))
OUTDIR = ROOT / "résultats" / "phase2b_pfam"

# Tool / DB locations are configurable via env vars (PFAM_DB, HMMSCAN_BIN);
# the defaults point at the known local copies (cf. ~/.claude/knowledge/tuberculosis.md).
# Stabilise by exporting, e.g.:
#   export PFAM_DB=/path/to/Pfam-A.hmm
#   export HMMSCAN_BIN=/path/to/hmmscan
_DEFAULT_HMMSCAN = ROOT.parent / "L8" / "eggnog-mapper-2.1.12" / "eggnogmapper" / "bin" / "hmmscan"
_DEFAULT_PFAM = ROOT.parent / "projets_abandonnes" / "Mycobacterium_sp_novel" / "data" / "pfam" / "Pfam-A.hmm"
HMMSCAN = Path(os.environ.get("HMMSCAN_BIN", _DEFAULT_HMMSCAN))
PFAM_DB = Path(os.environ.get("PFAM_DB", _DEFAULT_PFAM))


def load_xref() -> dict[str, dict]:
    with open(XREF) as fh:
        return {r["rv"]: r for r in csv.DictReader(fh, delimiter="\t")}


def parse_domtblout(path: Path) -> list[dict]:
    """Parse HMMER domtblout. One row per significant domain hit."""
    hits = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        # domtblout columns (space-separated, description may contain spaces -> 22 fixed + tail)
        hits.append({
            "pfam_name": f[0],
            "pfam_acc": f[1],
            "i_evalue": float(f[12]),
            "ali_from": int(f[17]),
            "ali_to": int(f[18]),
            "description": " ".join(f[22:]),
        })
    # keep best i-Evalue per Pfam accession
    best: dict[str, dict] = {}
    for h in hits:
        k = h["pfam_acc"]
        if k not in best or h["i_evalue"] < best[k]["i_evalue"]:
            best[k] = h
    return sorted(best.values(), key=lambda h: h["ali_from"])


def main(batch_file: str) -> None:
    if not HMMSCAN.exists() or not PFAM_DB.exists():
        sys.exit(f"Missing hmmscan ({HMMSCAN.exists()}) or Pfam-A ({PFAM_DB.exists()})")
    xref = load_xref()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rvs = [l.strip() for l in Path(batch_file).read_text().splitlines() if l.strip()]

    # Merge with any prior results so pfam.json accumulates across batches.
    pj = OUTDIR / "pfam.json"
    result = json.loads(pj.read_text()) if pj.exists() else {}
    for rv in rvs:
        row = xref.get(rv)
        if not row:
            print(f"[skip] {rv}"); continue
        faa = OUTDIR / f"{rv}.faa"
        faa.write_text(f">{rv}|{row['mtbc0']}\n{row['protein_mtbc0']}\n")
        domtbl = OUTDIR / f"{rv}.domtblout"
        cmd = [str(HMMSCAN), "--cut_ga", "--domtblout", str(domtbl), str(PFAM_DB), str(faa)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            print(f"[fail] {rv}: {proc.stderr.strip()[:200]}"); continue
        hits = parse_domtblout(domtbl)
        result[rv] = hits
        dom = ", ".join(f"{h['pfam_name']}({h['pfam_acc']},E={h['i_evalue']:.1e})" for h in hits) or "NONE"
        print(f"{rv:8s} {row['mtbc0']:14s} -> {dom}")

    (OUTDIR / "pfam.json").write_text(json.dumps(result, indent=2))
    print(f"\nWrote {OUTDIR / 'pfam.json'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/pilot_batch_01.txt")
