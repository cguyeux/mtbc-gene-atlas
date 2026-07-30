# MTBC Gene Atlas — web front-end

Successor to Mycobrowser: a continuously updated functional re-annotation of the
*Mycobacterium tuberculosis* complex gene set, anchored on the **MTBC0 ancestral
genome** (Harrison et al. 2024). Built on the same FastAPI + PostgreSQL + Docker
pattern as the sibling project `atlas_mtbc` (per-lineage), but centred on **genes**.

## What it serves

- A browsable catalogue of all genes joined between MTBC0 and H37Rv (`gene_xref.tsv`).
- Per-gene fiches with: legacy → PGAP → revised annotation, Pfam domains
  (`hmmscan --cut_ga`), optional Foldseek structural hits, ESM Atlas signal
  (exploratory), evidence, a graded **verdict** + confidence, and a **Sources**
  section citing the provenance of every field.
- Full-text-ish search, verdict / hypothetical / enriched filters, a JSON API
  (`/api/genes`), and a per-gene feedback form.

## Data contract

The site is decoupled from the pipeline by files under `content/`:

| File | Produced by | Role |
|---|---|---|
| `content/gene_xref.tsv` | `analyses/phase1_build_gene_table.py` | full catalogue (one row per joined gene) |
| `content/genes/<rv>.json` | `analyses/phase4_export_json.py` | enriched per-gene records |

Refresh after a new batch:

```bash
cd ..                       # annotation_mtbc/
python analyses/phase4_export_json.py data/pilot_batch_01.txt
cp data/gene_xref.tsv site/content/gene_xref.tsv
```

## Run locally (SQLite, no Docker)

```bash
python -m venv .venv && .venv/bin/pip install -e ".[dev]"
.venv/bin/python -m backend.ingest      # optional: pre-populate; startup also auto-ingests
./bin/dev                               # http://127.0.0.1:9624
```

## Run dockerised (PostgreSQL)

```bash
cp ../data/gene_xref.tsv content/gene_xref.tsv
cd deploy && docker compose up --build  # http://localhost:9624
```

Port `9624` = `NC_000962.3` + 1, chosen to coexist with `atlas_mtbc` on `9623`.

## Configuration (env vars, prefix `ANNOT_`)

- `ANNOT_DATABASE_URL` (default SQLite under `data/`)
- `ANNOT_AUTO_INGEST` (default `true`): populate the DB on startup if empty
- `ANNOT_GENE_XREF`, `ANNOT_GENES_CONTENT`: override data locations
