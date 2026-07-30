#!/usr/bin/env python3
"""phase2_esm_annotate.py -- ESM Atlas enrichment of MTBC0 ancestral proteins.

For each gene of a batch (list of Rv locus tags, one per line), takes the
*ancestral MTBC0* protein from data/gene_xref.tsv and queries the ESM Atlas
(EvolutionaryScale x BioHub) via the esm-atlas-cli client:

  * lookup_sequence  -> top SAE features (key: 'sae_features')
  * cluster          -> cluster size + % Pfam-characterised
  * similarity_search-> nearest proteins in ESM space (putative homologs)

Writes one JSON per gene under résultats/phase2_esm/<rv>.json (raw payloads
kept verbatim for transparency) plus a condensed view.

ESM SAE labels are EXPLORATORY (cf. mtbc-gene-function caveat): they orient a
hypothesis, they do not prove a domain. Corroborate with PGAP/orthology/
literature before writing a function.

Usage:
  python phase2_esm_annotate.py data/pilot_batch_01.txt
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ESM_SRC = ROOT.parent / ".claude" / "skills" / "esm-atlas-cli" / "src"
sys.path.insert(0, str(ESM_SRC))

XREF = ROOT / "data" / "gene_xref.tsv"
OUTDIR = ROOT / "résultats" / "phase2_esm"


def load_xref() -> dict[str, dict]:
    with open(XREF) as fh:
        return {r["rv"]: r for r in csv.DictReader(fh, delimiter="\t")}


def first(d: dict, *keys, default=None):
    for k in keys:
        if k in d and d[k] not in (None, ""):
            return d[k]
    return default


def condense(rv: str, row: dict, lk: dict | None, cl: dict | None, sim: dict | None) -> dict:
    feats = []
    if lk:
        for f in (lk.get("sae_features") or lk.get("topk_features") or [])[:12]:
            feats.append({
                "index": first(f, "index", "feature_index"),
                "activation": first(f, "activation", "value"),
                "label": first(f, "label", "description", default=""),
            })
    homologs, consensus = [], []
    if sim:
        for h in (sim.get("similar_proteins") or [])[:10]:
            homologs.append({
                "accession": first(h, "protein_accession", "accession"),
                "similarity": first(h, "similarity_score", "similarity", "score"),
                "len_aa": first(h, "sequence_length", "length"),
            })
        for f in (sim.get("top_features_across_results") or [])[:12]:
            consensus.append({
                "index": first(f, "index", "feature_index"),
                "label": first(f, "label", "description", default=""),
            })
    return {
        "rv": rv,
        "mtbc0": row["mtbc0"],
        "gene": row["gene"],
        "product_h37rv": row["product_h37rv"],
        "product_mtbc0_pgap": row["product_mtbc0_pgap"],
        "len_aa": int(row["len_aa"]),
        "protein_hash": first(lk or {}, "protein_hash", "hash"),
        "mean_plddt": first(lk or {}, "mean_plddt"),
        "cluster_pct_characterized": first(cl or {}, "pct_characterized",
                                           "cluster_pct_characterized"),
        "cluster_size": first(cl or {}, "size", "cluster_size"),
        "top_sae_features": feats,
        "consensus_sae_features": consensus,
        "esm_homologs": homologs,
    }


def main(batch_file: str) -> None:
    from esm_atlas_cli import EsmAtlasClient, EsmAtlasUnavailable

    xref = load_xref()
    client = EsmAtlasClient()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    rvs = [l.strip() for l in Path(batch_file).read_text().splitlines() if l.strip()]
    for rv in rvs:
        row = xref.get(rv)
        if not row:
            print(f"[skip] {rv}: not in gene_xref.tsv"); continue
        prot = row["protein_mtbc0"]
        lk = cl = sim = None
        status = "ok"
        try:
            lk = client.lookup_sequence(prot, topk_features=12)
            rep = first(lk, "cluster_rep_protein_hash", "cluster_rep_hash")
            h = first(lk, "hash", "protein_hash")
            if rep or h:
                try:                      # cluster is optional: 404 when the
                    cl = client.cluster(rep or h)  # protein has no referenced cluster
                except Exception:
                    pass
            try:
                sim = client.similarity_search(prot, topk_results=10)
            except Exception:
                pass
        except EsmAtlasUnavailable as e:
            status = f"esm_unavailable: {e}"

        cond = condense(rv, row, lk, cl, sim)
        cond["status"] = status
        payload = {"condensed": cond, "raw": {"lookup": lk, "cluster": cl, "similarity": sim}}
        (OUTDIR / f"{rv}.json").write_text(json.dumps(payload, indent=2, default=str))
        nf = len(cond["top_sae_features"])
        nh = len(cond["esm_homologs"])
        print(f"[{status:>4}] {rv:8s} {row['mtbc0']:14s} {nf} SAE feats, {nh} homologs, "
              f"cluster%char={cond['cluster_pct_characterized']}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/pilot_batch_01.txt")
