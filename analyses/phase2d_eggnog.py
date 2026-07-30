#!/usr/bin/env python3
"""phase2d_eggnog.py -- Controlled-vocabulary annotation of the whole proteome.

Runs eggNOG-mapper (diamond mode) on every gene of the catalogue and parses
the orthology transfer into a compact per-gene record: orthologous group,
COG functional category, EC number, KEGG KO / pathway / module, Gene Ontology
terms, CAZy family, eggNOG description and preferred name.

Unlike the Pfam (phase2b) and Foldseek (phase2c) layers -- which target the
residual 'dark' genes -- this layer lifts *all* 3906 genes at once with a
standard, citable, machine-readable vocabulary (EC/GO/KEGG/COG/CAZy). It is
the single highest-ROI enrichment for a Mycobrowser successor, and makes the
atlas interoperable with KEGG, QuickGO and the COG database.

Database note: eggNOG-mapper 2.1.12 hardcodes the dead host eggnogdb.embl.de;
the live mirror is eggnog5.embl.de (patched in download_eggnog_data.py on
2026-06-02). Set EGGNOG_DATA_DIR to the directory holding eggnog.db +
eggnog_proteins.dmnd (default = the install's data/ dir).

Inputs:  data/gene_xref_all.tsv (rv -> protein), the eggNOG 5.0.2 DB.
Output:  résultats/phase2d_eggnog/{proteome.faa, mtbc.emapper.annotations,
         eggnog.json}.  eggnog.json is keyed by rv and consumed by phase4.

Usage:   python phase2d_eggnog.py            # whole proteome (gene_xref_all)
         GENE_XREF=data/gene_xref.tsv python phase2d_eggnog.py
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
XREF = Path(os.environ.get("GENE_XREF", ROOT / "data" / "gene_xref_all.tsv"))
OUTDIR = ROOT / "résultats" / "phase2d_eggnog"

_EGGNOG_HOME = ROOT.parent / "L8" / "eggnog-mapper-2.1.12"
EMAPPER = Path(os.environ.get("EMAPPER_BIN", _EGGNOG_HOME / "emapper.py"))
DATA_DIR = Path(os.environ.get("EGGNOG_DATA_DIR", _EGGNOG_HOME / "data"))

# Columns of an eggNOG-mapper 2.1.12 .emapper.annotations file (after #query).
COLS = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
    "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
]


def build_proteome(rows: dict[str, dict]) -> Path:
    """One FASTA, one entry per Rv (dedup paralogues, skip empty)."""
    OUTDIR.mkdir(parents=True, exist_ok=True)
    faa = OUTDIR / "proteome.faa"
    n = 0
    with open(faa, "w") as fh:
        for rv, r in rows.items():
            seq = (r.get("protein_mtbc0") or "").strip().rstrip("*")
            if not seq or set(seq) <= {"-"}:
                continue
            fh.write(f">{rv}\n{seq}\n")
            n += 1
    print(f"[proteome] {n} proteins -> {faa}")
    return faa


def run_emapper(faa: Path) -> Path:
    ann = OUTDIR / "mtbc.emapper.annotations"
    if ann.exists() and os.environ.get("EGGNOG_REUSE"):
        print(f"[emapper] reuse {ann}")
        return ann
    cpu = str(os.cpu_count() or 4)
    cmd = [
        sys.executable, str(EMAPPER),
        "-i", str(faa), "--itype", "proteins",
        "-m", "diamond", "--data_dir", str(DATA_DIR),
        "-o", "mtbc", "--output_dir", str(OUTDIR),
        "--cpu", cpu, "--override", "--no_file_comments",
    ]
    print("[emapper]", " ".join(cmd))
    proc = subprocess.run(cmd, text=True)
    if proc.returncode != 0:
        sys.exit(f"emapper failed (exit {proc.returncode})")
    return ann


def _clean(v: str) -> str:
    return "" if v in ("-", "", None) else v


def _split(v: str) -> list[str]:
    v = _clean(v)
    return [x for x in v.split(",") if x] if v else []


def _narrow_og(field: str) -> str:
    """eggNOG_OGs = 'OG@tax|name,OG@tax|name,...' -> first (narrowest) OG id."""
    ogs = _split(field)
    return ogs[0].split("@")[0] if ogs else ""


def parse_annotations(ann: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for line in ann.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < len(COLS):
            f += [""] * (len(COLS) - len(f))
        d = dict(zip(COLS, f))
        rv = d["query"]
        kos = [k.replace("ko:", "") for k in _split(d["KEGG_ko"])]
        rec = {
            "og": _narrow_og(d["eggNOG_OGs"]),
            "cog_cat": _clean(d["COG_category"]),
            "description": _clean(d["Description"]),
            "preferred_name": _clean(d["Preferred_name"]),
            "ec": _split(d["EC"]),
            "kegg_ko": kos,
            "kegg_pathway": [p for p in _split(d["KEGG_Pathway"]) if p.startswith("map")],
            "kegg_module": _split(d["KEGG_Module"]),
            "cazy": _split(d["CAZy"]),
            "go": _split(d["GOs"]),
            "pfams": _split(d["PFAMs"]),
            "seed_ortholog": _clean(d["seed_ortholog"]),
            "evalue": _clean(d["evalue"]),
        }
        # Drop completely empty hits (eggNOG sometimes emits a query with no transfer).
        if any(rec[k] for k in ("og", "description", "ec", "kegg_ko", "go", "cog_cat")):
            out[rv] = rec
    return out


def main() -> None:
    if not EMAPPER.exists():
        sys.exit(f"emapper.py not found at {EMAPPER}")
    if not (DATA_DIR / "eggnog.db").exists():
        sys.exit(f"eggNOG DB missing: {DATA_DIR / 'eggnog.db'} "
                 f"(run download_eggnog_data.py first)")
    with open(XREF) as fh:
        rows: dict[str, dict] = {}
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.setdefault(r["rv"], r)  # first occurrence per Rv (H37Rv-centric)
    faa = build_proteome(rows)
    ann = run_emapper(faa)
    recs = parse_annotations(ann)
    (OUTDIR / "eggnog.json").write_text(json.dumps(recs, indent=2, ensure_ascii=False))

    n_ec = sum(1 for v in recs.values() if v["ec"])
    n_ko = sum(1 for v in recs.values() if v["kegg_ko"])
    n_cog = sum(1 for v in recs.values() if v["cog_cat"])
    n_cazy = sum(1 for v in recs.values() if v["cazy"])
    print(f"\n[eggnog] {len(recs)} genes annotated "
          f"({n_cog} COG cat, {n_ec} EC, {n_ko} KEGG KO, {n_cazy} CAZy)")
    print(f"Wrote {OUTDIR / 'eggnog.json'}")


if __name__ == "__main__":
    main()
