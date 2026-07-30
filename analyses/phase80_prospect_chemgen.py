#!/usr/bin/env python3
"""phase80 — PROSPECT chemical-genetic / druggability layer (P16.10).

Source: Bond AN et al. "Reference-based chemical-genetic interaction profiling to
elucidate small molecule mechanism of action." Nat Commun 2025;16:9673.
doi:10.1038/s41467-025-64662-x (PMID 41184262). Supplementary Data downloaded
(springer static-content) into résultats/phase80_prospect/S{3..8}.xlsx.

What we ingest, per gene (HONEST scope):
  - Hypomorph-panel membership: the PROSPECT collection is a set of TetON
    transcriptional-knockdown strains of *essential* genes, used as a sensitised
    background for chemical-genetic target deconvolution. A gene being in the
    panel = it is an essential/vulnerable target for which a validated knockdown
    tool strain exists. Table S4 ("strains and baseline doublings").
  - Baseline knockdown fitness (median doublings) + whether the strain was used
    in the PCL (phenotypic-cluster) target-ID analysis.
  - Druggability cross-ref: if the gene name is an annotated mechanism-of-action
    (MOA) target of reference compounds (Table S3), we record how many reference
    compounds phenocopy its inhibition = it is a chemically-validated drug target.

Scope reality (measured): 465 panel genes all map to the atlas; only 7 are
verdict=dark and 0 are H37Rv-hypothetical. So this is a *druggability / tool-strain*
layer, not a dark-matter resolver. It is written for all 465 panel genes; the MOA
druggability field is added wherever the gene name matches a reference MOA target.

Writes layer `prospect` in place into site/content/genes/Rv*.json. Read-only on
verdicts (never flips a gene). Idempotent.
"""
from __future__ import annotations
import json, re, glob, os, warnings
warnings.filterwarnings("ignore")
from openpyxl import load_workbook

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SUP = os.path.join(ROOT, "résultats", "phase80_prospect")
CG = os.path.join(ROOT, "site", "content", "genes")

CITATION = ("Bond AN et al., Nat Commun 2025;16:9673 (doi:10.1038/s41467-025-64662-x); "
            "PROSPECT chemical-genetic platform")

RV = re.compile(r"Rv\d{4}[A-Za-z]?")


def num(x):
    try:
        return round(float(x), 3)
    except (TypeError, ValueError):
        return None


def s(x):
    """Safely stringify an openpyxl cell value -> stripped str or None."""
    if x is None:
        return None
    v = str(x).strip()
    return v or None


def load_sheet(fn, sheet):
    wb = load_workbook(os.path.join(SUP, fn), read_only=True, data_only=True)
    ws = wb[sheet]
    rows = list(ws.iter_rows(values_only=True))
    hdr = [str(c).strip() if c is not None else "" for c in rows[0]]
    return hdr, rows[1:]


def main():
    # --- atlas gene_name -> rv (for MOA druggability match) ---
    name2rv, rv_set = {}, set()
    for f in glob.glob(os.path.join(CG, "Rv*.json")):
        d = json.load(open(f))
        rv = d.get("rv")
        rv_set.add(rv)
        gn = (d.get("gene") or "").strip()
        if gn:
            name2rv.setdefault(gn.lower(), rv)

    # --- S3: reference-compound MOA targets -> count per target token ---
    hdr, rows = load_sheet("S3.xlsx", "reference set")
    mi = hdr.index("Annotated MOA")
    moa_count: dict[str, int] = {}
    for r in rows:
        v = r[mi]
        if not v:
            continue
        for tok in str(v).split("|"):
            tok = tok.strip()
            if tok:
                moa_count[tok] = moa_count.get(tok, 0) + 1
    # map MOA target tokens that are actual atlas gene names
    moa_by_rv: dict[str, tuple[str, int]] = {}
    for tok, n in moa_count.items():
        rv = name2rv.get(tok.lower())
        if rv:
            moa_by_rv[rv] = (tok, n)

    # --- S4: hypomorph panel per gene ---
    hdr, rows = load_sheet("S4.xlsx", "strains and baseline doublings")
    ix = {h: i for i, h in enumerate(hdr)}
    per_rv: dict[str, dict] = {}
    for r in rows:
        locus = r[ix["H37Rv locus"]]
        if not locus:
            continue
        m = RV.search(str(locus))
        if not m:
            continue
        rv = m.group(0)
        cat = s(r[ix["Strain category"]]) or ""
        used = (s(r[ix["Strain used in PCL analysis"]]) or "").lower() == "yes"
        entry = {
            "strain": s(r[ix["Strain name"]]),
            "gene": s(r[ix["Gene"]]),
            "teton_promoter": s(r[ix["TetON promoter"]]),
            "baseline_doublings_median": num(r[ix["Median number of doublings"]]),
            "n_screens": int(ns) if (ns := s(r[ix["Number of screens included in strain pool"]])) and ns.isdigit() else None,
            "strain_category": cat,
            "used_in_pcl": used,
        }
        # keep the strain with most screens if a gene has several
        prev = per_rv.get(rv)
        if prev is None or (entry["n_screens"] or 0) > (prev["n_screens"] or 0):
            per_rv[rv] = entry

    # --- assemble + attach druggability, write in place ---
    written = 0
    dark_hits, moa_hits = [], []
    for rv, entry in per_rv.items():
        fn = os.path.join(CG, f"{rv}.json")
        if not os.path.exists(fn):
            continue
        d = json.load(open(fn))
        layer = dict(entry)
        layer["in_hypomorph_panel"] = True
        if rv in moa_by_rv:
            tok, n = moa_by_rv[rv]
            layer["drug_target_moa"] = tok
            layer["reference_compounds_phenocopying"] = n
            moa_hits.append((rv, tok, n))
        layer["source"] = CITATION
        d["prospect"] = layer
        json.dump(d, open(fn, "w"), ensure_ascii=False, indent=2)
        written += 1
        if d.get("verdict") == "dark":
            dark_hits.append(rv)

    print(f"phase80: couche prospect écrite sur {written} gènes (panel hypomorphe PROSPECT)")
    print(f"  cibles druggables (MOA nommée + composés de référence): {len(moa_hits)}")
    for rv, tok, n in sorted(moa_hits, key=lambda x: -x[2])[:12]:
        print(f"    {rv}  {tok:12s}  {n} composés de référence")
    print(f"  gènes dark dans le panel: {len(dark_hits)} -> {sorted(dark_hits)}")


if __name__ == "__main__":
    main()
