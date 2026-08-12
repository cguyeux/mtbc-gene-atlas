# MTBC Gene Atlas

A continuously updated, structure- and population-scale functional re-annotation
of the *Mycobacterium tuberculosis* complex (MTBC) gene set, anchored on the MTBC0
ancestral genome (Harrison et al. 2024). This repository holds the companion code,
the annotation pipeline and the derived per-gene records for the web resource served
at https://mtbc.gclab.fr (also reachable via the persistent, domain-independent
identifier https://w3id.org/mtbc-atlas/).

The project takes over from Mycobrowser (EPFL/UNIBE), the reference *Mycobacterium*
annotation database, which is no longer maintained. Its central aim is to verify,
reconstruct and complete the functional annotation of the MTBC proteome, with priority
on the large stock of genes labelled "hypothetical protein". Every proposed function
carries its supporting evidence, a dated provenance and a confidence level; a graded
verdict marks each gene as requalified, family-assigned or still dark.

## What this repository contains

- `analyses/` : the phased annotation pipeline (`phaseN_*.py`). Each phase adds one
  evidence layer per gene (orthology, structure, conservation, interaction network,
  essentiality, localisation, and so on) and writes back into the per-gene records.
- `site/` : the web application that serves the enriched fiches. FastAPI, SQLAlchemy
  and Jinja2, containerised for deployment (`site/deploy/`). See `site/README.md` for
  the data contract and configuration.
- `site/content/` : the derived data served by the site, namely the joined gene
  catalogue (`gene_xref.tsv`) and one enriched JSON record per gene
  (`content/genes/<locus>.json`). This is the output of the pipeline.

Heavy input datasets consumed by the pipeline (for example the MtbTnDB Tn-seq
compendium, the TFOE expression tables, the MTBC0 genome and the Mycobrowser
reference proteomes) are not tracked here. They are either fetched from their original
public sources by the pipeline scripts or archived alongside the manuscript on Zenodo
(DOI 10.5281/zenodo.20815247).

## Evidence layers

Each gene aggregates, where available: orthology (eggNOG COG/EC/KO/GO/CAZy), UniProt
curation, intra-MTBC conservation (pN/pS over more than 145,000 strains), predicted
structure (ESMFold and AlphaFold) with Foldseek structural leads and pLDDT gating,
the STRING interaction network, M-CSA catalytic sites, Tn-seq essentiality, CRISPRi
vulnerability, PaxDb proteomics, DeepTMHMM localisation, operon and regulon context,
iModulon expression, PTM, Regions of Difference, mutant phenotypes, Pfam domains, and
a field-by-field comparison against the legacy Mycobrowser record.

## Quick start (local, SQLite, no Docker)

```
cd site
python -m venv .venv && .venv/bin/pip install -e ".[dev]"
./bin/dev            # http://127.0.0.1:9624
```

The database is not shipped: on first start the site auto-ingests the records under
`site/content/` into a local SQLite database (`ANNOT_AUTO_INGEST=true` by default),
so the resource is browsable right after a clone. See `site/README.md` for the
dockerised (PostgreSQL) setup and the environment variables.

## Licensing

The pipeline and application code are released under the MIT licence (`LICENSE`).
The annotation records and derived data are released under CC-BY 4.0
(`DATA_LICENSE.md`). Redistributed third-party annotations retain the licences of
their sources: UniProt and STRING under CC-BY 4.0, Pfam, eggNOG, and the other
databases cited on each fiche.

## Citation

Guyeux C. From hypothetical to functional: a continuously updated, structure- and
population-scale re-annotation of the *Mycobacterium tuberculosis* complex gene set,
anchored on the MTBC0 ancestral genome. FEMTO-ST Institute, Besançon, France.
Archived on Zenodo, DOI 10.5281/zenodo.20815247.
