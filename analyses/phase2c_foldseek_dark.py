#!/usr/bin/env python3
"""phase2c_foldseek_dark.py -- Structural homology for genes Pfam can't resolve.

For genes left 'dark' after PGAP + Pfam + literature (DUF / conserved
hypotheticals), sequence-level profiles fail by definition. The next lever is
STRUCTURE: fetch the ESMFold model (via ESM Atlas, fold-on-demand) and run a
Foldseek search against a structural database. A confident structural hit to a
characterised fold often names a function that no sequence method could.

Targets: by default the genes whose curation verdict is 'dark'; or pass a batch
file (list of Rv) as argv[1].

Tools / DB are env-configurable:
  FOLDSEEK_BIN  (default tools/foldseek/bin/foldseek)
  FOLDSEEK_DB   (default tools/foldseek_db/pdb)   -- build with:
                tools/foldseek/bin/foldseek databases PDB tools/foldseek_db/pdb tmp

Output: résultats/phase2c_foldseek/<rv>.pdb (ESMFold model)
        résultats/phase2c_foldseek/foldseek.json  {rv: [hits]}
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT.parent / ".claude" / "skills" / "esm-atlas-cli" / "src"))

ESM_DIR = ROOT / "résultats" / "phase2_esm"
XREF = Path(os.environ.get("GENE_XREF", ROOT / "data" / "gene_xref.tsv"))
CURATION = ROOT / "data" / "pilot_curation.json"
OUTDIR = ROOT / "résultats" / "phase2c_foldseek"
FOLDSEEK = Path(os.environ.get("FOLDSEEK_BIN", ROOT / "tools" / "foldseek" / "bin" / "foldseek"))
FOLDSEEK_DB = Path(os.environ.get("FOLDSEEK_DB", ROOT / "tools" / "foldseek_db" / "pdb"))
TOPN = 8


def targets(argv: list[str]) -> list[str]:
    if len(argv) > 1:
        return [l.strip() for l in Path(argv[1]).read_text().splitlines() if l.strip()]
    cur = json.loads(CURATION.read_text())
    return [rv for rv, c in cur.items() if not rv.startswith("_") and c.get("verdict") == "dark"]


def load_xref() -> dict[str, dict]:
    with open(XREF) as fh:
        return {r["rv"]: r for r in csv.DictReader(fh, delimiter="\t")}


def esmfold_api(seq: str, timeout: int = 200) -> str | None:
    """Fold a sequence with the PUBLIC ESMFold API (server-side, no GPU needed).
    Limited to ~400 aa. Used when the ESM Atlas has no pre-folded structure."""
    if not seq or len(seq) > 400:
        return None
    try:
        req = urllib.request.Request("https://api.esmatlas.com/foldSequence/v1/pdb/",
                                     data=seq.encode(), method="POST")
        with urllib.request.urlopen(req, timeout=timeout) as r:
            pdb = r.read().decode()
        return pdb if pdb.lstrip().startswith(("HEADER", "ATOM", "REMARK")) else None
    except Exception:
        return None


def fetch_pdb(rv: str, xref: dict) -> str | None:
    """ESMFold model for a gene: (1) the ESM Atlas pre-folded structure (fast,
    cached) if available, else (2) the public ESMFold API (folds server-side,
    <=400 aa). Local ESMFold is not used (no GPU)."""
    seq = (xref.get(rv) or {}).get("protein_mtbc0", "")
    try:
        from esm_atlas_cli import EsmAtlasClient
        c = EsmAtlasClient()
        ef = ESM_DIR / f"{rv}.json"
        if ef.exists():
            h = json.loads(ef.read_text())["condensed"]["protein_hash"]
            pdb = c.lookup(h, fold_on_miss=True).get("pdb")
            if pdb:
                return pdb
        elif seq:
            pdb = c.lookup_sequence(seq, fold_on_miss=True).get("pdb")
            if pdb:
                return pdb
    except Exception:
        pass
    return esmfold_api(seq)


def run_foldseek(pdb_path: Path, rv: str) -> list[dict]:
    out_m8 = OUTDIR / f"{rv}.m8"
    tmp = OUTDIR / "tmp"
    # Permissive search: mycobacterial DUF folds rarely have a *significant* PDB
    # match. We report the best structural neighbours ranked by Foldseek
    # homology probability, and flag whether each clears a significance bar
    # (E < 0.01). Anything else is a suggestive fold hint, not an assignment.
    cmd = [str(FOLDSEEK), "easy-search", str(pdb_path), str(FOLDSEEK_DB), str(out_m8), str(tmp),
           "--format-output", "target,evalue,prob,alntmscore,theader", "-e", "10", "--max-seqs", "50"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"  [{rv}] foldseek failed: {proc.stderr.strip()[:200]}")
        return []
    hits = []
    for line in out_m8.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        ev = float(f[1])
        hits.append({"target": f[0], "evalue": ev, "prob": float(f[2]),
                     "tmscore": float(f[3]), "description": f[4],
                     "significant": ev < 0.01})
    hits.sort(key=lambda h: h["prob"], reverse=True)
    return hits[:TOPN]


def main(argv: list[str]) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    xref = load_xref()
    rvs = targets(argv)
    db_ready = FOLDSEEK_DB.with_suffix(".dbtype").exists() or Path(str(FOLDSEEK_DB) + ".dbtype").exists()
    if not db_ready:
        print(f"NOTE: Foldseek DB not found at {FOLDSEEK_DB} -- will fetch structures only.\n"
              f"Build it with: {FOLDSEEK} databases PDB {FOLDSEEK_DB} tmp")
    # Merge with any prior results so foldseek.json accumulates across batches.
    fj = OUTDIR / "foldseek.json"
    result = json.loads(fj.read_text()) if fj.exists() else {}
    for i, rv in enumerate(rvs, 1):
        pdb = fetch_pdb(rv, xref)
        if not pdb:
            print(f"{rv}: no structure"); continue
        pdb_path = OUTDIR / f"{rv}.pdb"
        pdb_path.write_text(pdb)
        if db_ready:
            hits = run_foldseek(pdb_path, rv)
            result[rv] = hits
            top = hits[0] if hits else None
            print(f"[{i}/{len(rvs)}] {rv}: {len(hits)} hits"
                  + (f" | {top['description'][:50]} (E={top['evalue']:.0e}, TM={top['tmscore']:.2f})" if top else ""))
        else:
            print(f"{rv}: structure saved ({len(pdb)} chars), foldseek skipped (no DB)")
        if i % 20 == 0:   # checkpoint so a long mass run survives interruption
            fj.write_text(json.dumps(result, indent=2))
    if result:
        fj.write_text(json.dumps(result, indent=2))
        print(f"\nWrote {fj}")


if __name__ == "__main__":
    main(sys.argv)
