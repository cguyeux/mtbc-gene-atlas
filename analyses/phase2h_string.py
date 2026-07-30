#!/usr/bin/env python3
"""phase2h_string.py -- Functional interaction network (STRING v12) layer.

Adds a *guilt-by-association* layer on top of the orthology (eggNOG, phase2d),
curation (UniProt, phase2e), domain (Pfam, phase2b) and structure (Foldseek,
phase2c) layers. Where those say "this gene IS in COG X / has this fold",
STRING says "this gene WORKS WITH that characterised operon" -- a signal none
of the other layers produces, and the one that lets us propose a citable
function for a hypothetical from its genomic-context neighbours.

Why STRING here, and the caveats that shape the code:
  * STRING v12.0 (active since 2023-07, von Mering/Jensen, EMBL/SIB), CC-BY 4.0,
    so it can be redistributed on the site with attribution.
  * STRING's own docs say the REST API is for *occasional* access; for "all
    interactions of a species" use the bulk files. So we download the per-
    organism flat file ONCE (taxon 83332 = H37Rv) and parse locally -- never
    hammer interaction_partners 3906x at 1/s.
  * Every edge is decomposed into 7 channels: neighborhood, fusion, cooccurence
    (= phylogenetic profiling), coexpression, experimental, database,
    textmining. The first three are prokaryote-strong, computational, and
    INDEPENDENT of eggNOG orthology / ESMFold structure.
  * The textmining channel is empty for true hypotheticals and merely echoes
    UniProt for known genes -> we recompute a `combined_no_tm` score (STRING's
    own noisy-OR with prior 0.041, minus the textmining channel) so a link's
    strength from data alone is visible. `context_driven` flags edges carried
    by the genomic-context channels.
  * Consequence: coverage of hypotheticals is sparser than of known genes --
    that is a quality filter, not a defect.

IDs: STRING uses `83332.Rv####` for H37Rv, so the rv key is the suffix after
the taxon prefix. The MTBC0 bridge is already handled by gene_xref.

Inputs:  data/gene_xref_all.tsv (rv, is_hypothetical_h37rv, product, gene),
         the STRING bulk files (downloaded on first run, cached).
Output:  résultats/phase2h_string/{83332.protein.links.detailed.v12.0.txt.gz,
         83332.protein.info.v12.0.txt.gz, string.json, meta.json}.
         string.json is keyed by rv and consumed by phase4.

Usage:   python phase2h_string.py                  # download (if needed) + build
         STRING_REUSE=1 python phase2h_string.py    # never re-download
         GENE_XREF=data/gene_xref.tsv python phase2h_string.py
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
XREF = Path(os.environ.get("GENE_XREF", ROOT / "data" / "gene_xref_all.tsv"))
OUTDIR = Path(os.environ.get("STRING_DATA_DIR", ROOT / "résultats" / "phase2h_string"))

TAXON = "83332"  # M. tuberculosis H37Rv
VERSION = "v12.0"
BASE = "https://stringdb-downloads.org/download"
LINKS_URL = f"{BASE}/protein.links.detailed.{VERSION}/{TAXON}.protein.links.detailed.{VERSION}.txt.gz"
INFO_URL = f"{BASE}/protein.info.{VERSION}/{TAXON}.protein.info.{VERSION}.txt.gz"
LINKS_GZ = OUTDIR / f"{TAXON}.protein.links.detailed.{VERSION}.txt.gz"
INFO_GZ = OUTDIR / f"{TAXON}.protein.info.{VERSION}.txt.gz"

# STRING per-organism detailed file channels, in column order after the two ids.
CHANNELS = ["neighborhood", "fusion", "cooccurence", "coexpression",
            "experimental", "database", "textmining"]
CONTEXT_CHANNELS = ["neighborhood", "fusion", "cooccurence"]
PRIOR = 0.041  # STRING combination prior

# Thresholds (STRING scores are 0-1000).
MEDIUM = 400       # keep a partner edge at/above medium confidence
HIGH = 700         # high confidence
CONTEXT_MIN = 400  # a context channel at/above this => context_driven edge
TOP_N = 15         # partners stored per gene


def download(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[cache] {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
        return
    if os.environ.get("STRING_REUSE"):
        sys.exit(f"STRING_REUSE set but {dest} is missing")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    print(f"[download] {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "annotation_mtbc/phase2h (research)"})
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as fh:
        while chunk := r.read(1 << 20):
            fh.write(chunk)
    tmp.rename(dest)
    print(f"[download] -> {dest} ({dest.stat().st_size/1e6:.1f} MB)")


def combine(scores: list[float]) -> int:
    """STRING noisy-OR combination of channel scores (each 0..1000) -> 0..1000.

    Reproduces STRING's published formula: prior-correct each channel, combine
    by noisy-OR, then re-add the prior. Used to recompute `combined_no_tm`.
    """
    prod = 1.0
    for s in scores:
        sp = (s / 1000.0 - PRIOR) / (1.0 - PRIOR)
        if sp > 0:
            prod *= (1.0 - sp)
    tot_nopr = 1.0 - prod
    tot = tot_nopr + PRIOR * (1.0 - tot_nopr)
    return round(tot * 1000)


def load_xref() -> dict[str, dict]:
    with open(XREF) as fh:
        rows: dict[str, dict] = {}
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.setdefault(r["rv"], r)
    return rows


def parse_info() -> dict[str, dict]:
    """STRING protein.info -> {rv: {gene, annotation}} (rv = id minus taxon)."""
    info: dict[str, dict] = {}
    with gzip.open(INFO_GZ, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            sid, name, _size, annot = f[0], f[1], f[2], f[3]
            rv = sid.split(".", 1)[1] if "." in sid else sid
            info[rv] = {"gene": name, "annotation": annot}
    return info


def parse_links() -> tuple[dict[str, list], int, float]:
    """Stream the links file -> {rv: [partner_edge, ...]} above MEDIUM.

    Each edge keeps per-channel scores, the raw combined, and a recomputed
    combined excluding the textmining channel. Also returns a sanity check of
    the recomputed combined against the file's combined column.
    """
    partners: dict[str, list] = {}
    n_edges = 0
    err_sum = 0.0
    err_n = 0
    with gzip.open(LINKS_GZ, "rt") as fh:
        header = fh.readline().split()
        # Expected: protein1 protein2 <7 channels> combined_score
        ci = {name: header.index(name) for name in CHANNELS if name in header}
        comb_idx = header.index("combined_score")
        for line in fh:
            f = line.split()
            if len(f) < comb_idx + 1:
                continue
            combined = int(f[comb_idx])
            if combined < MEDIUM:
                continue  # half the file is below medium; skip early
            p1 = f[0].split(".", 1)[1] if "." in f[0] else f[0]
            p2 = f[1].split(".", 1)[1] if "." in f[1] else f[1]
            ch = {name: int(f[ci[name]]) for name in CHANNELS if name in ci}
            # Recompute combined w/o textmining (data-only strength).
            no_tm = combine([v for name, v in ch.items() if name != "textmining"])
            # Sanity: full recompute vs file combined.
            full = combine(list(ch.values()))
            err_sum += abs(full - combined); err_n += 1
            context = max(ch.get(c, 0) for c in CONTEXT_CHANNELS)
            edge = {
                "p2": p2,
                "channels": ch,
                "combined": combined,
                "combined_no_tm": no_tm,
                "context": context,
                "context_driven": context >= CONTEXT_MIN,
            }
            partners.setdefault(p1, []).append(edge)
            n_edges += 1
    mae = err_sum / err_n if err_n else float("nan")
    return partners, n_edges, mae


def build_record(rv: str, edges: list, xref: dict, info: dict) -> dict:
    """Per-gene STRING record: top partners + a guilt-by-association anchor."""
    # Sort by data-only strength first, then raw combined.
    edges = sorted(edges, key=lambda e: (e["combined_no_tm"], e["combined"]), reverse=True)
    out_partners = []
    for e in edges[:TOP_N]:
        p2 = e["p2"]
        prow = xref.get(p2, {})
        pinfo = info.get(p2, {})
        gene = (prow.get("gene") or pinfo.get("gene") or "").strip()
        product = (prow.get("product_h37rv") or pinfo.get("annotation") or "").strip()
        hypo = str(prow.get("is_hypothetical_h37rv", "")).strip().lower() == "true"
        out_partners.append({
            "rv": p2,
            "gene": gene if gene and gene != p2 else "",
            "product": product[:140],
            "hypothetical": hypo,
            "combined": e["combined"],
            "combined_no_tm": e["combined_no_tm"],
            "context_driven": e["context_driven"],
            "channels": e["channels"],
        })

    # Guilt-by-association anchor: best context-driven, data-strong partner that
    # is itself a CHARACTERISED gene (not hypothetical). This is the citable seed
    # of a function hypothesis for a hypothetical query gene.
    anchor = None
    for p in out_partners:
        if p["hypothetical"] or not p["context_driven"]:
            continue
        if p["combined_no_tm"] >= MEDIUM:
            anchor = {
                "rv": p["rv"],
                "gene": p["gene"],
                "product": p["product"],
                "combined_no_tm": p["combined_no_tm"],
                "confidence": "high" if p["combined_no_tm"] >= HIGH else "medium",
                "channels": p["channels"],
            }
            break

    return {
        "string_id": f"{TAXON}.{rv}",
        "n_partners": len(edges),
        "n_context_driven": sum(1 for e in edges if e["context_driven"]),
        "partners": out_partners,
        "anchor": anchor,
    }


def main() -> None:
    download(LINKS_URL, LINKS_GZ)
    download(INFO_URL, INFO_GZ)
    xref = load_xref()
    info = parse_info()
    partners, n_edges, mae = parse_links()
    print(f"[links] {n_edges} edges >= {MEDIUM} over {len(partners)} genes "
          f"(combined recompute MAE = {mae:.1f})")
    if mae > 15:
        print("[warn] combined recompute drifts from file column; "
              "channel order or formula may differ -- inspect before trusting no_tm")

    records: dict[str, dict] = {}
    for rv, edges in partners.items():
        if rv not in xref:
            continue  # STRING id not in our catalogue (rare)
        records[rv] = build_record(rv, edges, xref, info)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "string.json").write_text(json.dumps(records, indent=2, ensure_ascii=False))

    # Stats + demo on hypotheticals that gain a citable anchor.
    n_anchor = sum(1 for r in records.values() if r["anchor"])
    hypo_rvs = {rv for rv, r in xref.items()
                if str(r.get("is_hypothetical_h37rv", "")).strip().lower() == "true"}
    hypo_anchored = [(rv, records[rv]) for rv in hypo_rvs
                     if rv in records and records[rv]["anchor"]]
    hypo_anchored.sort(key=lambda kv: kv[1]["anchor"]["combined_no_tm"], reverse=True)

    print(f"\n[string] {len(records)} genes with >=1 partner "
          f"({n_anchor} have a non-textmining context anchor)")
    print(f"[string] {len(hypo_anchored)}/{len(hypo_rvs)} hypotheticals get a "
          f"guilt-by-association anchor toward a characterised gene")
    print("\n--- demo: top hypotheticals with a citable function hypothesis ---")
    for rv, r in hypo_anchored[:15]:
        a = r["anchor"]
        prod_self = xref[rv].get("product_h37rv", "")
        tag = a["gene"] or a["rv"]
        print(f"  {rv} ({prod_self[:32]:32s}) -> {tag} "
              f"[{a['confidence']}, no_tm={a['combined_no_tm']}] {a['product'][:46]}")
    print(f"\nWrote {OUTDIR / 'string.json'}")

    (OUTDIR / "meta.json").write_text(json.dumps({
        "source": "STRING", "version": VERSION, "taxon": TAXON,
        "license": "CC-BY 4.0", "url": "https://string-db.org",
        "thresholds": {"medium": MEDIUM, "high": HIGH,
                       "context_min": CONTEXT_MIN, "top_n": TOP_N},
        "channels": CHANNELS, "context_channels": CONTEXT_CHANNELS,
        "combined_recompute_mae": round(mae, 2),
        "n_genes": len(records), "n_edges_kept": n_edges,
        "n_anchor": n_anchor,
        "n_hypotheticals": len(hypo_rvs),
        "n_hypotheticals_anchored": len(hypo_anchored),
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
