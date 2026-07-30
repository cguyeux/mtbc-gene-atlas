"""Deterministic MTBC resistance profiler -- self-contained engine for the site.

Ports the reconciled catalogue-matching of the Resistance_antibio project
(`analyses/resistance_profile.py` + `variant_match.py`) into a dependency-light
module (pure stdlib, no pandas) so it can ship inside the gene-atlas container as
a second tool. Given a strain's variants (SPDI list, VCF text, or a single HGVS
mutation), it returns, per drug, a verdict "R" (>=1 catalogued R-associated
determinant carried) with the determinants hit.

The matching is RECONCILED: SNP exact + decomposition of observed MNVs into their
component SNPs + indel exact. This captures the representation blind spots the
upstream project corrected (notably the gyrA MNV for fluoroquinolones) that a raw
SPDI equality would miss.

Catalogue: content/resistance/catalogue_consolide.tsv (WHO 2nd ed. 2023 + tb-profiler
+ CRyPTIC empirical, consolidated). Only `call == "R-associated"` rows drive a verdict.

RESEARCH TOOL, NOT A CLINICAL DIAGNOSTIC. The catalogue only covers resistance
attributable to a known target mutation; absence of a determinant is not evidence
of susceptibility.
"""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

CHROM = "NC_000962.3"
CATALOGUE = Path(__file__).resolve().parent.parent / "content" / "resistance" / "catalogue_consolide.tsv"

# Display order: first line, fluoroquinolones, injectables, second line, new drugs.
DRUG_ORDER = ["isoniazid", "rifampicin", "ethambutol", "pyrazinamide", "streptomycin",
              "moxifloxacin", "levofloxacin", "amikacin", "kanamycin", "capreomycin",
              "ethionamide", "bedaquiline", "clofazimine", "linezolid", "delamanid"]

_AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "E": "Glu",
        "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
        "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
        "Y": "Tyr", "V": "Val", "*": "*"}
_ONE_LETTER = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY*])$")


# --- SPDI helpers (vendored from variant_match.py) -------------------------

def parse_spdi(spdi: str):
    _, p0, ref, alt = spdi.split(":")
    return int(p0), ref, alt


def is_snp(spdi: str) -> bool:
    _, ref, alt = parse_spdi(spdi)
    return len(ref) == 1 and len(alt) == 1


def decompose(spdi: str):
    """Component single-base SNPs of a substitution (SNP -> itself; MNV -> differing bases)."""
    p0, ref, alt = parse_spdi(spdi)
    if len(ref) != len(alt):
        return []
    return [f"{CHROM}:{p0 + i}:{ref[i]}:{alt[i]}"
            for i in range(len(ref)) if ref[i] != alt[i]]


# --- Catalogue (loaded once, cached) ---------------------------------------

@lru_cache(maxsize=1)
def _catalogue():
    det: set[str] = set()
    rows_by_spdi: dict[str, list[dict]] = defaultdict(list)
    drugs: set[str] = set()
    variant_index: dict[str, str] = {}     # canonical "gene_variant" -> spdi (any call)
    n_rows = 0
    with open(CATALOGUE, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            spdi = r.get("spdi", "")
            if not spdi:
                continue
            n_rows += 1
            if r.get("variant"):
                variant_index.setdefault(_canon(r["variant"]), spdi)
            if r.get("call") == "R-associated":
                det.add(spdi)
                drugs.add(r["drug"])
                rows_by_spdi[spdi].append({
                    "gene": r.get("gene", ""), "variant": r.get("variant", ""),
                    "drug": r["drug"], "grade": r.get("source_grade", ""),
                    "source": r.get("source", ""),
                })
    snp_det = {s for s in det if is_snp(s)}
    ordered = [d for d in DRUG_ORDER if d in drugs] + sorted(d for d in drugs if d not in DRUG_ORDER)
    return {"det": det, "snp_det": snp_det, "rows_by_spdi": dict(rows_by_spdi),
            "drugs": ordered, "variant_index": variant_index, "n_determinant_spdi": len(det)}


def _canon(variant: str) -> str:
    """Normalise a mutation label so user input and catalogue labels compare equal.

    'katG S315T' / 'katG_p.Ser315Thr' / 'katG_S315T' all map to 'katg_p.ser315thr'.
    """
    v = variant.strip().replace(" ", "_")
    parts = v.split("_", 1)
    gene = parts[0].lower()
    rest = parts[1] if len(parts) > 1 else ""
    rest = rest.replace("p.", "").replace("P.", "")
    m = _ONE_LETTER.match(rest)
    if m:                                  # 1-letter shorthand (S315T) -> 3-letter (Ser315Thr)
        a, num, b = m.groups()
        rest = f"{_AA3[a]}{num}{_AA3[b]}"
    rest = rest.lower()
    return f"{gene}_{rest}" if rest else gene


# --- Reconciled matching (ported from resistance_profile.py) ----------------

def matched_determinants(strain_spdis, det, snp_det):
    """{catalogued determinant SPDI hit -> 'exact'|'MNV'}, reconciled (SNP/MNV/indel)."""
    matched: dict[str, str] = {}
    for s in strain_spdis:
        try:
            _, ref, alt = parse_spdi(s)
        except Exception:
            continue
        if len(ref) == len(alt):                       # SNP or MNV substitution
            if s in det:
                matched[s] = "exact"
            multibase = len(ref) > 1
            for comp in decompose(s):
                if comp in snp_det:
                    matched.setdefault(comp, "MNV" if multibase else "exact")
        else:                                          # indel
            if s in det:
                matched[s] = "exact"
    return matched


def profile(strain_spdis):
    """set(SPDI) -> {drug -> {verdict, determinants:[{gene,variant,spdi,grade,source,match}]}}."""
    cat = _catalogue()
    matched = matched_determinants(strain_spdis, cat["det"], cat["snp_det"])
    by_drug: dict[str, dict[str, dict]] = defaultdict(dict)   # drug -> {spdi -> det}
    for spdi, mode in matched.items():
        for row in cat["rows_by_spdi"].get(spdi, []):
            by_drug[row["drug"]][spdi] = {**row, "spdi": spdi, "match": mode}
    out = {}
    for drug in cat["drugs"]:
        dets = sorted(by_drug.get(drug, {}).values(), key=lambda d: d["spdi"])
        out[drug] = {"verdict": "R" if dets else "S", "determinants": dets}
    return out


# --- Input parsing ----------------------------------------------------------

_SPDI_RE = re.compile(rf"{CHROM}:\d+:[ACGTNacgtn]*:[ACGTNacgtn]*")


def parse_spdi_text(text: str) -> set[str]:
    return set(_SPDI_RE.findall(text or ""))


def parse_vcf_text(text: str) -> set[str]:
    """VCF (H37Rv, CHROM=NC_000962.3) -> SPDI 0-based (POS-1, project convention)."""
    spdis = set()
    for line in (text or "").splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 5:
            f = line.split()
        if len(f) < 5:
            continue
        pos, ref, alt = f[1], f[3], f[4]
        for a in alt.split(","):
            if a in (".", "<NON_REF>", "*", ""):
                continue
            try:
                spdis.add(f"{CHROM}:{int(pos) - 1}:{ref}:{a}")
            except ValueError:
                continue
    return spdis


def resolve_mutation(text: str):
    """A single mutation -> (spdi, recognised:bool). Accepts a raw SPDI or an HGVS label."""
    text = (text or "").strip()
    if _SPDI_RE.fullmatch(text):
        return text, True
    spdi = _catalogue()["variant_index"].get(_canon(text))
    return (spdi, True) if spdi else (None, False)


def catalogue_stats() -> dict:
    cat = _catalogue()
    return {"n_determinants": cat["n_determinant_spdi"], "n_drugs": len(cat["drugs"]),
            "drugs": cat["drugs"]}
