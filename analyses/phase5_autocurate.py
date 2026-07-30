#!/usr/bin/env python3
"""phase5_autocurate.py -- Rule-based semi-automatic curation.

Derives a verdict + function + evidence for each gene from the deterministic
signals already computed (MTBC0 PGAP product, Pfam domains via hmmscan,
Foldseek structural hits), WITHOUT a per-gene manual pass. Output goes to a
SEPARATE file (data/autocuration.json) so it never overwrites the hand-curated
data/pilot_curation.json (which always takes priority downstream).

Each entry is flagged:
  auto: true
  needs_review: true   when a human should look (requalified claims, or a
                       verdict resting only on a structural Foldseek hint,
                       or a PGAP-hypothetical that a domain/structure rescued)

Verdict logic (graded, conservative):
  - requalified     PGAP gives a specific named function (an enzyme/role,
                    not a '... family'/'... domain-containing' label).
  - family_assigned PGAP gives a family/domain, OR a non-DUF Pfam domain is
                    present, OR a significant Foldseek hit to a characterised
                    fold (structure-based).
  - dark            PGAP 'hypothetical' (or DUF-only) with no functional Pfam
                    and no significant structural hit.

Usage:
  python phase5_autocurate.py data/some_batch.txt      # specific genes
  python phase5_autocurate.py --all-hypothetical       # every hypothetical not already hand-curated
"""
from __future__ import annotations

import csv
import json
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
XREF = Path(os.environ.get("GENE_XREF", ROOT / "data" / "gene_xref.tsv"))
PFAM_FILE = ROOT / "résultats" / "phase2b_pfam" / "pfam.json"
FOLDSEEK_FILE = ROOT / "résultats" / "phase2c_foldseek" / "foldseek.json"
MANUAL = ROOT / "data" / "pilot_curation.json"
OUT = ROOT / "data" / "autocuration.json"

GENERIC = ("family protein", "domain-containing protein", "domain-containing",
           " family", "-like protein", "-containing protein", "binding protein",
           "secretion-associated", "transcriptional regulator", "regulatory protein",
           " subunit", "domain nuclease", "alpha/beta hydrolase", "transporter permease",
           "abc transporter", "primosomal protein")
# 'putative/probable/...' = an uncertain call -> family-level, not a confident requalification.
UNCERTAIN = ("putative", "probable", "possible", "uncharacterized")
# Generic localisation/descriptor products that are family-level ONLY when they are the
# WHOLE product (exact match) -- so 'membrane protein' is family, but
# 'membrane protein insertase YidC' stays a named requalification.
EXACT_GENERIC = {"membrane protein", "integral membrane protein", "conserved membrane protein",
                 "possible membrane protein", "secreted protein", "exported protein",
                 "lipoprotein", "conserved protein", "putative protein", "cell wall protein"}

# --- paralogue-of-fold guard ---------------------------------------------------
# A Foldseek hit to a NAMED enzyme does not mean THIS gene has that function: the
# real enzyme may be another, already-named locus (e.g. Rv3268 hits 'acetyl-CoA
# synthetase' but acsA is Rv3667). We index named functions/gene names of the whole
# genome and refuse to promote a dark gene to a function already owned elsewhere.
ENZ_RE = re.compile(r"([a-z][a-z0-9'/-]{2,})\s+([a-z]{4,}ase)\b")
KNOWN_ENZ: dict[str, set] = {}   # "homoserine kinase" -> {rv,...}
KNOWN_GENE: dict[str, str] = {}  # "pkni" -> rv


def _functions(text: str) -> set[str]:
    return {f"{m.group(1)} {m.group(2)}" for m in ENZ_RE.finditer(text.lower())}


def build_known(xref_all: dict) -> None:
    for rv, r in xref_all.items():
        prod = (r.get("product_mtbc0_pgap") or r.get("product_h37rv") or "").lower()
        if prod and "hypothetical" not in prod:
            for f in _functions(prod):
                KNOWN_ENZ.setdefault(f, set()).add(rv)
        gn = (r.get("gene") or "").strip().lower()
        if len(gn) >= 3:
            KNOWN_GENE[gn] = rv


def paralogue_of(desc: str, rv: str) -> str | None:
    """If the structural hit's function is already a different named locus, return it."""
    for f in _functions(desc):
        owners = KNOWN_ENZ.get(f, set()) - {rv}
        if owners:
            return f"{f} = {sorted(owners)[0]}"
    dl = desc.lower()
    for gn, owner in KNOWN_GENE.items():
        if owner != rv and re.search(r"\b" + re.escape(gn) + r"\b", dl):
            return f"{gn} = {owner}"
    return None


def load_xref() -> dict[str, dict]:
    with open(XREF) as fh:
        return {r["rv"]: r for r in csv.DictReader(fh, delimiter="\t")}


def is_duf(name: str) -> bool:
    return name.upper().startswith("DUF") or name.upper().startswith("UPF")


def cap(s: str) -> str:
    s = (s or "").strip()
    return s[:1].upper() + s[1:] if s else s


def sig_fold_hits(hits: list) -> list:
    out = []
    for h in hits or []:
        # A dark gene is promoted to family only on a genuinely strong structural
        # match: either E-significant AND a real fold overlap (TM >= 0.4), or a
        # near-certain same-fold (prob >= 0.95 AND TM >= 0.8). E-significant with a
        # tiny TM (partial alignment) or high prob with low TM stays a text note,
        # not a verdict change -- this tracks manual judgment.
        tm = h.get("tmscore", 0)
        strong = (h.get("significant") and tm >= 0.4) or (h.get("prob", 0) >= 0.95 and tm >= 0.8)
        desc = (h.get("description") or "").lower()
        named = not any(w in desc for w in ("uncharacteriz", "hypothetical", "duf", "designed"))
        if strong and named:
            out.append(h)
    return out


def classify(row: dict, pfam: list, fold: list) -> dict:
    h37 = row.get("product_h37rv", "")
    pgap_raw = (row.get("product_mtbc0_pgap") or "").strip()
    unanchored = not pgap_raw                     # no MTBC0 anchor -> annotate on H37Rv product
    pgap = (pgap_raw or h37).strip()
    pgap_l = pgap.lower()
    legacy = h37
    is_hypo = pgap_l in ("", "hypothetical protein")
    is_duf_product = ("duf" in pgap_l and "domain-containing" in pgap_l)
    func_pfam = [d for d in pfam if not is_duf(d["pfam_name"])]
    duf_pfam = [d for d in pfam if is_duf(d["pfam_name"])]
    sig = sig_fold_hits(fold)

    needs_review = False
    paralogue = None
    if is_hypo or is_duf_product:
        if func_pfam:
            verdict, conf = "family_assigned", "medium"   # rescued by a real Pfam domain -- reliable
        elif sig:
            paralogue = paralogue_of(sig[0].get("description", ""), row.get("rv", ""))
            if paralogue:
                verdict, conf = "dark", "low"             # paralogue-of-fold trap: do NOT propagate
            else:
                verdict, conf = "family_assigned", "low"
                needs_review = True                       # structure-only promotion -- still to validate
        else:
            verdict, conf = "dark", "low"
    else:
        # PGAP carries an annotation. Generic descriptors and uncertain wording
        # ('putative/probable') are a family-level call, not a confident requalification.
        if (pgap_l in EXACT_GENERIC or any(m in pgap_l for m in GENERIC)
                or any(u in pgap_l for u in UNCERTAIN)):
            verdict = "family_assigned"
            conf = "medium" if func_pfam else "low"
        else:
            verdict = "requalified"          # a specific, confidently named PGAP function (reviewed: trusted)
            conf = "high" if func_pfam else "medium"

    # ---- function text ----
    pfam_str = ", ".join(f"{d['pfam_name']} ({d['pfam_acc']})" for d in pfam)
    if verdict == "dark" and paralogue:
        fn = (f"Conserved hypothetical; no reliable functional assignment. A Foldseek hit "
              f"resembles {sig[0]['description'][:55]}, but that function is encoded elsewhere "
              f"in the genome ({paralogue}) -- treated as a paralogue-of-fold artefact, not propagated.")
    elif verdict == "dark":
        dom = (f"DUF domain(s) {', '.join(d['pfam_name'] for d in duf_pfam)}"
               if duf_pfam else "no recognised domain")
        fn = f"Conserved hypothetical protein; {dom}. Function unknown."
        if fold:
            b = fold[0]
            fn += (f" Foldseek best (non-significant) hit: {b['description'][:70]} "
                   f"(prob {b.get('prob',0):.2f}, TM {b.get('tmscore',0):.2f}).")
    elif is_hypo and not func_pfam and sig:
        b = sig[0]
        fn = (f"No Pfam domain above threshold; Foldseek indicates a fold similar to "
              f"{b['description'][:80]} (prob {b.get('prob',0):.2f}, TM {b.get('tmscore',0):.2f}). "
              f"Structure-based, putative.")
    elif is_hypo and func_pfam:
        fn = f"Contains {pfam_str} domain(s); putative function inferred from the domain architecture."
    else:
        fn = cap(pgap) + "."
        if func_pfam:
            fn += f" Pfam: {pfam_str}."

    if unanchored:
        evidence = [f"Annotation from H37Rv (no MTBC0 1:1 anchor; H37Rv protein used): {legacy or 'hypothetical protein'}"]
    else:
        evidence = [f"Legacy H37Rv annotation: {legacy}",
                    f"MTBC0 PGAP product: {pgap or 'hypothetical protein'}"]
    if pfam:
        evidence.append("Pfam (hmmscan --cut_ga): "
                        + ", ".join(f"{d['pfam_name']} {d['pfam_acc']} (E={d['i_evalue']:.0e})" for d in pfam))
    if fold:
        b = fold[0]
        evidence.append(f"Foldseek best: {b['description'][:80]} "
                        f"(prob {b.get('prob',0):.2f}, E={b.get('evalue',0):.0e}, TM={b.get('tmscore',0):.2f})")
    evidence.append("(auto-curated by rules from PGAP + Pfam + Foldseek; not hand-reviewed)")

    return {"verdict": verdict, "confidence": conf, "function_revised": fn,
            "evidence": evidence, "references": [], "auto": True, "needs_review": needs_review}


def main(argv: list[str]) -> None:
    xref = load_xref()
    pfam_all = json.loads(PFAM_FILE.read_text()) if PFAM_FILE.exists() else {}
    fold_all = json.loads(FOLDSEEK_FILE.read_text()) if FOLDSEEK_FILE.exists() else {}
    manual = {k for k in json.loads(MANUAL.read_text()) if not k.startswith("_")}

    if argv and argv[0] == "--all-hypothetical":
        rvs = [rv for rv, r in xref.items()
               if r.get("is_hypothetical_h37rv") == "True" and rv not in manual]
    elif argv:
        rvs = [l.strip() for l in Path(argv[0]).read_text().splitlines() if l.strip()]
    else:
        sys.exit("usage: phase5_autocurate.py <batch.txt> | --all-hypothetical")

    # Build the paralogue-of-fold index from the FULL genome catalogue.
    all_path = ROOT / "data" / "gene_xref_all.tsv"
    src = all_path if all_path.exists() else XREF
    with open(src) as fh:
        build_known({r["rv"]: r for r in csv.DictReader(fh, delimiter="\t")})

    out = json.loads(OUT.read_text()) if OUT.exists() else {}
    counts = {"requalified": 0, "family_assigned": 0, "dark": 0}
    review = traps = 0
    for rv in rvs:
        row = xref.get(rv)
        if not row:
            continue
        entry = classify(row, pfam_all.get(rv, []), fold_all.get(rv, []))
        out[rv] = entry
        counts[entry["verdict"]] += 1
        review += int(entry["needs_review"])
        traps += "paralogue-of-fold" in entry["function_revised"]
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(f"Auto-curated {len(rvs)} genes -> {OUT.name}")
    print(f"  requalified={counts['requalified']}  family_assigned={counts['family_assigned']}  "
          f"dark={counts['dark']}  | needs_review={review}  paralogue_traps_caught={traps}")


if __name__ == "__main__":
    main(sys.argv[1:])
