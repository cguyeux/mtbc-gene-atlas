#!/usr/bin/env python3
"""phase2e_uniprot.py -- Curated cross-reference layer from UniProt/SwissProt.

eggNOG-mapper (phase2d) transfers a *computed* controlled vocabulary by
orthology. For a model organism like H37Rv, UniProt additionally carries a
*manually curated* layer: reviewed (SwissProt) function text, reviewed EC
numbers, GO, KEGG cross-references and an explicit protein-existence evidence
level. This phase harvests that layer for the whole proteome so each fiche can
show curated evidence next to the computed one (and so a curator can see when
the two disagree).

The UniProt reference proteome of M. tuberculosis H37Rv is UP000001584; the
Rv locus tag is the "ordered locus name" (gene_oln). When several UniProt
entries map to one Rv, the reviewed (SwissProt) entry wins.

Inputs:  UniProt REST (rest.uniprot.org), no local DB needed.
Output:  résultats/phase2e_uniprot/{uniprot_raw.tsv, uniprot.json}.
         uniprot.json is keyed by rv and consumed by phase4.

Usage:   python phase2e_uniprot.py            # whole H37Rv proteome
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "résultats" / "phase2e_uniprot"
PROTEOME = "UP000001584"  # M. tuberculosis H37Rv reference proteome

FIELDS = ["accession", "gene_oln", "protein_name", "ec", "go_id",
          "cc_function", "xref_kegg", "protein_existence", "reviewed"]
STREAM = (
    "https://rest.uniprot.org/uniprotkb/stream?"
    f"query=proteome:{PROTEOME}&format=tsv&fields={','.join(FIELDS)}"
)

_RV = re.compile(r"\bRv\d{4}[A-Za-z]?\b")
_ECO = re.compile(r"\s*\{ECO:[^}]*\}")           # evidence tags
_FUNC = re.compile(r"^FUNCTION:\s*")
_PUBMED = re.compile(r"\s*\(PubMed:[^)]*\)")
_GO = re.compile(r"GO:\d{7}")


def fetch(raw: Path) -> Path:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    if raw.exists() and raw.stat().st_size > 0:
        print(f"[uniprot] reuse {raw} ({raw.stat().st_size} bytes)")
        return raw
    print(f"[uniprot] streaming proteome {PROTEOME} ...")
    proc = subprocess.run(
        ["curl", "-sS", "--fail", "--max-time", "300", "-o", str(raw), STREAM],
        text=True, capture_output=True,
    )
    if proc.returncode != 0 or not raw.exists() or raw.stat().st_size == 0:
        sys.exit(f"UniProt fetch failed: {proc.stderr.strip()[:300]}")
    return raw


def clean_function(text: str) -> str:
    text = _FUNC.sub("", text)
    text = _ECO.sub("", text)
    text = _PUBMED.sub("", text)
    # UniProt joins multiple FUNCTION comments with '; '. Keep it readable.
    return text.strip().strip(".").strip()[:600]


def split_list(v: str) -> list[str]:
    return [x.strip() for x in v.split(";") if x.strip()] if v else []


def parse(raw: Path) -> dict[str, dict]:
    lines = raw.read_text().splitlines()
    if not lines:
        return {}
    header = lines[0].split("\t")
    idx = {name: i for i, name in enumerate(header)}
    out: dict[str, dict] = {}
    for line in lines[1:]:
        f = line.split("\t")
        if len(f) < len(header):
            f += [""] * (len(header) - len(f))

        def col(key: str) -> str:
            # UniProt TSV header labels differ from the field codes; match loosely.
            for h, i in idx.items():
                hl = h.lower()
                if key == "oln" and "locus" in hl:
                    return f[i]
                if key == "acc" and h == "Entry":
                    return f[i]
                if key == "name" and "protein names" in hl:
                    return f[i]
                if key == "ec" and hl == "ec number":
                    return f[i]
                if key == "go" and hl.startswith("gene ontology"):
                    return f[i]
                if key == "func" and "function" in hl:
                    return f[i]
                if key == "kegg" and "kegg" in hl:
                    return f[i]
                if key == "pe" and "existence" in hl:
                    return f[i]
                if key == "rev" and hl == "reviewed":
                    return f[i]
            return ""

        oln = col("oln")
        rvs = _RV.findall(oln)
        if not rvs:
            continue
        reviewed = col("rev").strip().lower() == "reviewed"
        rec = {
            "acc": col("acc"),
            "reviewed": reviewed,
            "protein_name": col("name").split(" (")[0].strip(),
            "ec": [e for e in split_list(col("ec")) if e],
            "function": clean_function(col("func")),
            "go": _GO.findall(col("go")),
            "kegg": [k for k in split_list(col("kegg"))],
            "evidence": col("pe").strip(),
        }
        for rv in rvs:
            prev = out.get(rv)
            # Reviewed entry wins; otherwise keep the first / the richer one.
            if prev is None or (reviewed and not prev["reviewed"]):
                out[rv] = rec
    return out


def main() -> None:
    raw = fetch(OUTDIR / "uniprot_raw.tsv")
    recs = parse(raw)
    (OUTDIR / "uniprot.json").write_text(json.dumps(recs, indent=2, ensure_ascii=False))
    n_rev = sum(1 for v in recs.values() if v["reviewed"])
    n_func = sum(1 for v in recs.values() if v["function"])
    n_ec = sum(1 for v in recs.values() if v["ec"])
    print(f"\n[uniprot] {len(recs)} Rv mapped "
          f"({n_rev} reviewed/SwissProt, {n_func} with curated function, {n_ec} with EC)")
    print(f"Wrote {OUTDIR / 'uniprot.json'}")


if __name__ == "__main__":
    main()
