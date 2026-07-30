#!/usr/bin/env python3
"""phase83 — proteogenomics microprotein catalogue (P16.5, light entry).

Source: de Souza EV, Dalberto PF, ..., Saghatelian A, Bizarro CV.
"Large-scale proteogenomics characterization of microproteins in Mycobacterium
tuberculosis." Sci Rep 2024;14:31186 (doi:10.1038/s41598-024-82465-w). OPEN ACCESS.
Supplementary .xlsx downloaded via the springer static-content OA path into
résultats/phase83_microproteins/S{1..10}.xlsx.

Scope (perimeter EXPANSION, separate track): the H37Rv annotation excluded small
ORFs (smORFs). This ingests the MS-PROVEN microproteome only: Tier T1+T2 (72
microproteins, all "Spectrum Forest = hc" = high-confidence peptide-spectrum
evidence). T3+T4 (1350, "lc") are prediction-only and NOT ingested. Each entry is
tagged:
  - tORF = novel intergenic microprotein (true new gene, does NOT overlap any
    annotated CDS) -> 26 of the 72.
  - gORF = ORF overlapping / internal to an annotated gene (alternative reading
    frame / dual coding) -> 46 of the 72; cross-linked to its host locus.

These are a SEPARATE track: "microprotein, existence proven (MS), function unknown".
They are NOT counted in the 3906 canonical CDS nor in the 220 dark stock.

Merges per microprotein: coordinates/strand, start codon, Shine-Dalgarno, folding
free energy, sequence/length (S8 Supplementary Data 1); predicted localization
(S8 Domains); essentiality (S10); cross-species conservation (S7); host/surrounding
genes for gORFs (S1).

Output: résultats/phase83_microproteins/microprotein_catalog.{tsv,json} (read-only
on the atlas; no fiche is modified here — integration into the site is a separate
step once the display architecture is agreed).
"""
from __future__ import annotations
import json, os, re, warnings
warnings.filterwarnings("ignore")
from openpyxl import load_workbook

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
D = os.path.join(ROOT, "résultats", "phase83_microproteins")
GFF = os.path.join(ROOT, "..", "investigate_phylo", "resources", "NC_000962.3.gff3")
CONTENT = os.path.join(ROOT, "site", "content", "microproteins.json")
CITATION = ("de Souza EV et al., Sci Rep 2024;14:31186 "
            "(doi:10.1038/s41598-024-82465-w); proteogenomics microprotein compendium")


def load_h37rv_cds():
    """Return [(start, end, strand, locus_tag)] for all H37Rv CDS (ground-truth for overlap)."""
    cds = []
    for line in open(GFF):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "CDS":
            continue
        m = re.search(r"locus_tag=([^;]+)", p[8])
        cds.append((int(p[3]), int(p[4]), p[6], m.group(1) if m else "?"))
    return cds


def classify_overlap(a, b, strand, cds):
    """Classify an ORF [a,b]/strand against the reference CDS set by geometry (NOT the source label)."""
    hits = [c for c in cds if not (c[1] < a or c[0] > b)]
    if not hits:
        return "intergenic", None
    same = [c for c in hits if c[2] == strand]
    if same:
        return "sense_overlap", same[0][3]
    return "antisense", hits[0][3]


def sheet(fn, name):
    ws = load_workbook(os.path.join(D, fn), data_only=True)[name]
    rows = list(ws.iter_rows(values_only=True))
    hdr = [str(c).strip() if c is not None else "" for c in rows[0]]
    return hdr, rows[1:]


def norm_id(x):
    """Normalise smORF ids: strip stray brackets/quotes/tabs, collapse gORF__ -> gORF_."""
    s = re.sub(r"[\"\]\[\t].*$", "", str(x)).strip()
    return s.replace("gORF__", "gORF_")


def parse_coords(s):
    m = re.match(r"\s*(\d+)\s*-\s*(\d+)\s*", str(s))
    if not m:
        return None, None, None
    a, b = int(m.group(1)), int(m.group(2))
    if a <= b:
        return a, b, "+"
    return b, a, "-"          # reverse strand encoded as start>end in the source


def main():
    # --- master table S8 Supplementary Data 1 ---
    h, rows = sheet("S8.xlsx", "Supplementary Data 1")
    ix = {k: h.index(k) for k in h}
    cat = {}
    for r in rows:
        tier = str(r[ix["Tier"]])
        if tier not in ("T1", "T2"):
            continue
        sid = norm_id(r[ix["smORFs"]])
        start, end, strand = parse_coords(r[ix["Genome Coordinates"]])
        seq = str(r[ix["Microprotein Sequence"]] or "").strip()
        cat[sid] = {
            "id": sid,
            "kind": "novel_intergenic" if sid.startswith("tORF") else "overlapping_gene",
            "tier": tier,
            "ms_evidence": "high-confidence peptide-spectrum match (Spectrum Forest hc)",
            "detection_method": "MS peptide-spectrum (compendium)",
            "start": start, "end": end, "strand": strand,
            "length_aa": len(seq),
            "start_codon": str(r[ix["Start Codon"]] or "").strip() or None,
            "shine_dalgarno": str(r[ix["Shine Dalgarno"]] or "").strip() or None,
            "folding_free_energy": r[ix["Free Energy"]],
            "sequence": seq,
            "source": CITATION,
        }

    # --- localization (S8 Domains) ---
    h, rows = sheet("S8.xlsx", "Domains")
    li = {k: h.index(k) for k in h}
    for r in rows:
        sid = norm_id(r[0])
        if sid in cat:
            loc = str(r[li.get("Localization", 1)] or "").strip()
            rel = str(r[li.get("Reliability", 2)] or "").strip()
            if loc:
                cat[sid]["localization"] = loc
                cat[sid]["localization_reliability"] = rel or None

    # --- essentiality (S10); TnSeq codes: ES=essential, GD=growth-defect, GA=growth-advantage ---
    ESS_DECODE = {"ES": "essential", "GD": "growth-defect when disrupted",
                  "GA": "growth-advantage when disrupted"}
    h, rows = sheet("S10.xlsx", "Supplementary Table 1 fixed")
    ei = {k: h.index(k) for k in h}
    ecol = "essentiality" if "essentiality" in ei else h[-1]
    for r in rows:
        sid = norm_id(r[0])
        if sid in cat:
            code = str(r[ei[ecol]] or "").strip()
            if code:
                cat[sid]["essentiality_code"] = code
                cat[sid]["essentiality"] = ESS_DECODE.get(code, code)

    # --- host / surrounding genes for gORFs (S1) ---
    h, rows = sheet("S1.xlsx", "Supplementary Data 2 fixed")
    si = {k: h.index(k) for k in h}
    for r in rows:
        sid = norm_id(r[0])
        if sid in cat:
            host = str(r[si.get("Overlapping gene name", 2)] or "").strip()
            hostid = str(r[si.get("Overlapping gene id", 1)] or "").strip()
            if host and host.lower() not in ("none", "na", ""):
                cat[sid]["host_gene"] = host
                m = re.search(r"Rv\w+", hostid)
                cat[sid]["host_rv"] = m.group(0) if m else None
            sur = str(r[si.get("Surrounding gene names", 4)] or "").strip()
            if sur and sur.lower() not in ("none", "na", ""):
                cat[sid]["surrounding_genes"] = sur

    # --- conservation (S7): mark taxa where the microprotein has a hit ---
    wb = load_workbook(os.path.join(D, "S7.xlsx"), data_only=True)
    for sh in wb.sheetnames:
        taxon = sh.split("_")[0]
        for row in wb[sh].iter_rows(values_only=True):
            if not row or not row[0]:
                continue
            sid = norm_id(row[0])
            if sid in cat:
                cat[sid].setdefault("conserved_in", set()).add(taxon)
    for v in cat.values():
        if "conserved_in" in v:
            v["conserved_in"] = sorted(v["conserved_in"])

    # --- overlap classification by GEOMETRY vs H37Rv CDS (ground truth, not the source label) ---
    cds = load_h37rv_cds()
    for it in cat.values():
        if it["start"] is None:
            it["overlap_class"], it["overlap_locus"] = "unknown", None
            continue
        oc, loc = classify_overlap(it["start"], it["end"], it["strand"], cds)
        it["overlap_class"], it["overlap_locus"] = oc, loc

    # --- write ---
    items = sorted(cat.values(), key=lambda x: (x["overlap_class"] != "intergenic", x["id"]))
    with open(os.path.join(D, "microprotein_catalog.json"), "w") as f:
        json.dump(items, f, ensure_ascii=False, indent=2)
    cols = ["id", "kind", "tier", "start", "end", "strand", "length_aa", "start_codon",
            "localization", "essentiality", "host_gene", "host_rv"]
    with open(os.path.join(D, "microprotein_catalog.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for it in items:
            f.write("\t".join(str(it.get(c, "")) for c in cols) + "\n")
    # baked into the site image (separate microproteome track)
    os.makedirs(os.path.dirname(CONTENT), exist_ok=True)
    with open(CONTENT, "w") as f:
        json.dump(items, f, ensure_ascii=False, indent=2)

    from collections import Counter
    oc = Counter(it["overlap_class"] for it in items)
    n_ess = sum(1 for it in items if it.get("essentiality_code") == "ES")
    n_loc = sum(1 for it in items if it.get("localization"))
    n_cons = sum(1 for it in items if it.get("conserved_in"))
    print(f"phase83: catalogue microprotéines MS-prouvées (T1+T2) = {len(items)}")
    print(f"  classification par géométrie H37Rv: {dict(oc)}")
    print(f"  avec localisation prédite: {n_loc} | avec essentialité: {sum(1 for it in items if it.get('essentiality'))} "
          f"(dont essentielles ES: {n_ess}) | avec conservation: {n_cons}")
    print(f"  écrit: {D}/microprotein_catalog.json/.tsv + {CONTENT}")
    print("\n  === microprotéines essentielles (existence MS + TnSeq ES) ===")
    for it in items:
        if it.get("essentiality_code") == "ES":
            print(f"    {it['id']} ({it['kind']}, {it['length_aa']} aa) loc={it.get('localization')} "
                  f"host={it.get('host_gene')} ess={it.get('essentiality')}")


if __name__ == "__main__":
    main()
