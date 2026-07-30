#!/usr/bin/env python3
"""phase86b — add Mycobrowser-only genes to the canonical set (P16.5d, foundation stage).

The atlas was built on RefSeq NC_000962.3 (3906 CDS). Mycobrowser (the resource this
atlas succeeds) annotates 4031 CDS: 133 non-pseudogene protein-coding genes are absent
from the atlas. 55 are mobile-element/transposase fragments and 10 are split-gene halves
(RefSeq excludes these by design); the remaining 68 "serious" genes (28 named/functional,
35 conserved-hypothetical, 5 PE/PPE; 47 A-suffixed small genes) are a genuine coverage gap.

CG chose FULL enrichment. That is a multi-stage project (the pipeline is 90 phases, MTBC0-
anchored via RefSeq NP_ ids these genes lack; and the heavy layers need hmmscan/PFAM_DB,
eggNOG 39GB, ESMFold, Foldseek, DeepTMHMM/biolib, intra-MTBC conservation over 145k strains
— not all active in this shell). This script does the FOUNDATION stage, feasible now and
required by every later layer:
  - translate the protein from the H37Rv genome at the Mycobrowser coordinates
  - build a canonical fiche with the rich Mycobrowser-NATIVE annotation (product, function,
    functional category, Pfam, EC, GO, comments, UniProt AC) + Biopython ProtParam
  - register the gene in the catalogue (gene_xref.tsv) and as an enriched fiche
  - verdict from the Mycobrowser product (functional -> family_assigned, hypothetical -> dark),
    auto=true, needs_review=true, and a curation_note flagging provenance + pending heavy layers

Heavy layers (structure, orthology by eggNOG, intra-MTBC conservation, STRING, Foldseek,
localization, essentiality, proteomics, ...) remain a tracked follow-up (P16.5d-cont).

Idempotent. Read the H37Rv genome from investigate_phylo resources.
"""
from __future__ import annotations
import csv, json, os, re, warnings
warnings.filterwarnings("ignore")
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MB = os.path.join(ROOT, "data", "mycobrowser", "H37Rv.txt")
GENOME = os.path.join(ROOT, "..", "investigate_phylo", "resources", "NC_000962.3.fasta")
XREF = os.path.join(ROOT, "site", "content", "gene_xref.tsv")
CG = os.path.join(ROOT, "site", "content", "genes")
MB_REF = {"authors": "Kapopoulou A, Lew JM, Cole ST", "year": 2011,
          "title": "The MycoBrowser portal", "journal": "Tuberculosis 91:8-13"}


def load_genome():
    seq = []
    for line in open(GENOME):
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq)


def translate(genome, start, stop, strand):
    sub = genome[start - 1:stop]                 # 1-based inclusive
    s = Seq(sub)
    if strand == "-":
        s = s.reverse_complement()
    aa = str(s.translate(table=11))
    return aa.rstrip("*")


def protparam(aa):
    clean = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", aa)
    if len(clean) < 5:
        return None
    p = ProteinAnalysis(clean)
    return {
        "length_aa": len(aa),
        "mw_kda": round(p.molecular_weight() / 1000, 1),
        "pi": round(p.isoelectric_point(), 2),
        "gravy": round(p.gravy(), 3),
        "instability_index": round(p.instability_index(), 1),
        "instability_class": "stable" if p.instability_index() < 40 else "unstable",
        "aromaticity": round(p.aromaticity(), 3),
        "source": "Biopython ProtParam (ExPASy method)",
        "refs": [{"authors": "Gasteiger E, Hoogland C, Gattiker A, et al.", "year": 2005,
                  "title": "Protein Identification and Analysis Tools on the ExPASy Server"}],
    }


def parse_pfam(cell):
    """Mycobrowser PFAM cell -> [{pfam_name, pfam_acc, source}]."""
    out = []
    for m in re.finditer(r"(PF\d{5})(?:\.\d+)?\s*([^;,|]*)", str(cell or "")):
        out.append({"source": "mycobrowser", "pfam_acc": m.group(1),
                    "pfam_name": (m.group(2) or "").strip() or m.group(1)})
    return out


def is_mobile_or_fragment(product):
    p = (product or "").lower()
    if re.search(r"transposase|integrase|insertion sequence|\bIS\d|mobile|recombinase", p):
        return True
    if "fragment" in p:
        return True
    return False


def main():
    genome = load_genome()
    import glob as _glob
    # Build the set of atlas genes, EXCLUDING this phase's own prior additions so the script
    # is re-runnable (overwrites its P16.5d fiches instead of skipping them as "already present").
    atlas_up = set()
    for f in _glob.glob(os.path.join(CG, "Rv*.json")):
        try:
            note = json.load(open(f)).get("curation_note", "")
        except Exception:
            note = ""
        if not str(note).startswith("Added P16.5d"):
            atlas_up.add(os.path.basename(f)[:-5].upper())
    rows = list(csv.DictReader(open(MB), delimiter="\t"))

    selected = []
    for r in rows:
        if r["Feature"] != "CDS":
            continue
        if r["Locus"].upper() in atlas_up:
            continue
        if str(r.get("Is_Pseudogene", "")).strip().lower() in ("yes", "true", "1"):
            continue
        if is_mobile_or_fragment(r["Product"]):
            continue                              # 55 mobile + 10 fragment excluded by design
        selected.append(r)

    # --- build fiches ---
    xref_new = []
    written = 0
    for r in selected:
        rv = r["Locus"]
        try:
            start, stop = int(r["Start"]), int(r["Stop"])
        except ValueError:
            continue
        strand = r["Strand"] if r["Strand"] in ("+", "-") else "+"
        aa = translate(genome, start, stop, strand)
        product = (r["Product"] or "").strip() or "hypothetical protein"
        func = (r["Function"] or "").strip()
        is_hyp = ("hypothetical" in product.lower()) or ("unknown" in product.lower())
        # A Mycobrowser "Function" like "Unknown" / "Function unknown[, but ...]" is NOT a real
        # function; treat it as absent (anti-oversell: a conserved hypothetical must stay dark).
        real_func = bool(func) and not func.lower().startswith(
            ("unknown", "function unknown", "conserved hypothetical", "hypothetical"))
        ec = (r["Enzyme Classification"] or "").strip()
        verdict = "family_assigned" if (real_func or not is_hyp) else "dark"

        fiche = {
            "schema_version": 1, "rv": rv, "mtbc0": None, "np_id": None,
            "gene": (r["Name"] or "").strip() or None,
            "len_aa": len(aa), "strand": strand,
            "start_mtbc0": start, "end_mtbc0": stop,   # H37Rv coords (no MTBC0 anchor for these)
            "product_h37rv": product, "product_mtbc0_pgap": None,
            "function_revised": (func if real_func else (None if is_hyp else product)),
            "verdict": verdict, "confidence": "low",
            "auto": True, "needs_review": True,
            "protein_mtbc0": aa or None,
            "is_hypothetical": is_hyp,
            "curation_note": ("Added P16.5d: Mycobrowser-annotated gene absent from RefSeq "
                              "NC_000962.3 (the atlas' original 3906-CDS reference). Foundation "
                              "record from the Mycobrowser release + translated H37Rv sequence; "
                              "heavy enrichment layers (structure, orthology, conservation, "
                              "interaction) pending."),
            "mycobrowser": {
                "legacy_product": product, "legacy_function": func or None,
                "legacy_ec": ec or None,
                "functional_category": (r["Functional_Category"] or "").strip() or None,
                "comments": (r["Comments"] or "").strip() or None,
                "uniprot_ac": (r["UniProt_AC"] or "").strip() or None,
                "go": (r["Gene Ontology"] or "").strip() or None,
                "source": "Mycobrowser release 5 (H37Rv)", "refs": [MB_REF],
            },
            "funccat": ({"category": (r["Functional_Category"] or "").strip(),
                         "source": "TubercuList functional category (via Mycobrowser release 5)"}
                        if (r["Functional_Category"] or "").strip() else None),
            "protparam": protparam(aa),
        }
        dom = parse_pfam(r["PFAM"])
        if dom:
            fiche["domains"] = dom
        fiche = {k: v for k, v in fiche.items() if v is not None}
        json.dump(fiche, open(os.path.join(CG, f"{rv}.json"), "w"), ensure_ascii=False, indent=2)
        written += 1
        xref_new.append({
            "rv": rv, "mtbc0": "", "np_id": "", "gene": fiche.get("gene", ""),
            "product_h37rv": product, "product_mtbc0_pgap": "",
            "is_hypothetical_h37rv": str(is_hyp), "len_aa": str(len(aa)),
            "strand": strand, "start_mtbc0": str(start), "end_mtbc0": str(stop),
            "protein_mtbc0": aa,
        })

    # --- append to catalogue (idempotent: drop any pre-existing rows for these rv) ---
    with open(XREF) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        existing = [ln for ln in fh if ln.strip()]
    new_rv = {x["rv"] for x in xref_new}
    kept = [ln for ln in existing if ln.split("\t", 1)[0] not in new_rv]
    with open(XREF, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.writelines(kept)
        for x in xref_new:
            fh.write("\t".join(x.get(c, "") for c in header) + "\n")

    from collections import Counter
    verd = Counter(json.load(open(os.path.join(CG, f"{x['rv']}.json")))["verdict"] for x in xref_new)
    print(f"phase86b: {written} gènes Mycobrowser ajoutés au set canonique (fondation)")
    print(f"  verdicts: {dict(verd)}")
    print(f"  catalogue gene_xref.tsv -> {len(kept) + len(xref_new)} lignes")
    with_dom = sum(1 for x in xref_new if json.load(open(os.path.join(CG, f"{x['rv']}.json"))).get("domains"))
    print(f"  avec domaines Pfam (Mycobrowser): {with_dom} | avec ProtParam: {written}")


if __name__ == "__main__":
    main()
