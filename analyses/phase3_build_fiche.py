#!/usr/bin/env python3
"""phase3_build_fiche.py -- Assemble the enriched gene fiche (Mycobrowser successor).

Merges, per gene of a batch:
  * identity + ancestral MTBC0 protein  (data/gene_xref.tsv)
  * legacy H37Rv annotation             (the 'hypothetical protein' we resolve)
  * MTBC0 PGAP re-annotation             (first requalification layer)
  * Pfam domains, hmmscan --cut_ga       (résultats/phase2b_pfam/pfam.json) -- traceable
  * ESM Atlas SAE features / homologs    (résultats/phase2_esm/<rv>.json)   -- exploratory
  * curated verdict + literature         (data/pilot_curation.json)

Produces one Markdown fiche per gene under résultats/fiches/<rv>.md and a
batch INDEX.md. Every fiche carries an explicit verdict + confidence AND a
'Sources' section that always cites the provenance of each field (ancestral
sequence, product, domains, ESM signal, literature) -- so no fiche is ever
left without a citation, even when no primary functional study exists.

Usage:  python phase3_build_fiche.py data/pilot_batch_01.txt
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
XREF = ROOT / "data" / "gene_xref.tsv"
ESM_DIR = ROOT / "résultats" / "phase2_esm"
PFAM_FILE = ROOT / "résultats" / "phase2b_pfam" / "pfam.json"
CURATION = ROOT / "data" / "pilot_curation.json"
OUTDIR = ROOT / "résultats" / "fiches"

VERDICT_BADGE = {
    "requalified": "RESOLVED",
    "family_assigned": "FAMILY ASSIGNED",
    "dark": "STILL UNKNOWN",
}
# Standing data sources, cited on every fiche.
SRC_MTBC0 = ("Harrison LB et al. (2024), *An imputed ancestral reference genome for the "
             "MTBC...*, doi:[10.1101/2023.09.07.556366](https://doi.org/10.1101/2023.09.07.556366); "
             "repo github.com/lukebharrison/MTBC0")
SRC_ESM = ("ESM Atlas, EvolutionaryScale x BioHub Protein Atlas "
           "(ESMC sparse-autoencoder features) -- exploratory, not a validated domain call")


def load_xref() -> dict[str, dict]:
    with open(XREF) as fh:
        return {r["rv"]: r for r in csv.DictReader(fh, delimiter="\t")}


def wrap_seq(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def fiche_md(rv: str, row: dict, esm: dict, pfam: list, cur: dict) -> str:
    cond = esm.get("condensed", {})
    feats = cond.get("top_sae_features", [])
    homologs = cond.get("esm_homologs", [])
    sim_max = max((h["similarity"] for h in homologs if h.get("similarity")), default=None)
    badge = VERDICT_BADGE.get(cur.get("verdict", ""), "?")
    pfam_str = (", ".join(f"{h['pfam_name']} ({h['pfam_acc']})" for h in pfam)
                if pfam else "no domain above the Pfam gathering threshold")
    L = []
    L.append(f"# {rv}  ({row.get('gene') or '-'})")
    L.append("")
    L.append(f"**Verdict:** {badge} &nbsp;|&nbsp; **Confidence:** {cur.get('confidence','?')}")
    L.append("")
    L.append("## Identity")
    L.append("")
    L.append("| Field | Value |")
    L.append("|---|---|")
    L.append(f"| H37Rv locus tag | {rv} |")
    L.append(f"| MTBC0 locus tag | {row['mtbc0']} |")
    L.append(f"| RefSeq anchor | {row['np_id']} |")
    L.append(f"| Gene name | {row.get('gene') or '-'} |")
    L.append(f"| Protein length (MTBC0) | {row['len_aa']} aa |")
    L.append(f"| MTBC0 coordinates | {row['start_mtbc0']}..{row['end_mtbc0']} ({row['strand']}) |")
    L.append("")
    L.append("## Annotation: from legacy to revised")
    L.append("")
    L.append("| Source | Product / function |")
    L.append("|---|---|")
    L.append(f"| Legacy (H37Rv / Mycobrowser) | {row['product_h37rv']} |")
    L.append(f"| MTBC0 PGAP re-annotation | {row['product_mtbc0_pgap'] or '-'} |")
    L.append(f"| **Revised (this work)** | {cur.get('function_revised','-')} |")
    L.append("")
    L.append("## Domains (Pfam, hmmscan --cut_ga)")
    L.append("")
    if pfam:
        L.append("| Pfam | Accession | i-Evalue | Residues |")
        L.append("|---|---|---|---|")
        for h in pfam:
            L.append(f"| {h['pfam_name']} | {h['pfam_acc']} | {h['i_evalue']:.1e} "
                     f"| {h['ali_from']}-{h['ali_to']} |")
    else:
        L.append("No Pfam-A domain above the gathering threshold (genuinely uncharacterised "
                 "at the profile level; next step: HHpred profile-profile or Foldseek on the "
                 "ESMFold structure).")
    L.append("")
    L.append("## Evidence")
    L.append("")
    for e in cur.get("evidence", []):
        L.append(f"- {e}")
    L.append("")
    L.append("## ESM Atlas signal (exploratory)")
    L.append("")
    L.append(f"Ancestral MTBC0 protein hash `{cond.get('protein_hash','-')}`; "
             f"{len(homologs)} ESM-space neighbours"
             + (f" (max similarity {sim_max:.3f})" if sim_max else "") + ".")
    L.append("")
    if feats:
        L.append("Top SAE features (activation, label) -- orienting indices, not validated domains:")
        L.append("")
        for f in feats[:6]:
            act = f.get("activation")
            act_s = f"{act:.2f}" if isinstance(act, (int, float)) else str(act)
            L.append(f"- `{f.get('index')}` ({act_s}) {f.get('label') or '-'}")
        L.append("")
    # --- Sources: always present, one bullet per provenance channel ---
    L.append("## Sources")
    L.append("")
    L.append(f"- Ancestral sequence & coordinates: {SRC_MTBC0}")
    L.append(f"- Product annotation: NCBI PGAP on MTBC0; legacy from H37Rv NC_000962.3 "
             f"(RefSeq {row['np_id']})")
    L.append(f"- Domain assignment: Pfam-A via hmmscan --cut_ga -- {pfam_str}")
    L.append(f"- Sequence-level signal: {SRC_ESM}")
    refs = cur.get("references", [])
    if refs:
        for r in refs:
            doi = f" doi:[{r['doi']}](https://doi.org/{r['doi']})" if r.get("doi") else ""
            pmid = f" PMID:{r['pmid']}" if r.get("pmid") else ""
            L.append(f"- Primary literature: {r.get('authors','')} ({r.get('year','')}). "
                     f"*{r.get('title','')}* {r.get('journal','')}.{doi}{pmid}")
    else:
        L.append("- Primary literature: none located (tbmonitor 2021-2026 empty for this gene; "
                 "founder-era search via lit-review pending). Annotation rests on the "
                 "domain/homology sources above.")
    L.append("")
    L.append("## Ancestral MTBC0 protein sequence")
    L.append("")
    L.append("```")
    L.append(f">{row['mtbc0']}|{rv}|{row.get('gene') or ''}")
    L.append(wrap_seq(row["protein_mtbc0"]))
    L.append("```")
    L.append("")
    return "\n".join(L)


def main(batch_file: str) -> None:
    xref = load_xref()
    curation = json.loads(CURATION.read_text())
    pfam_all = json.loads(PFAM_FILE.read_text()) if PFAM_FILE.exists() else {}
    rvs = [l.strip() for l in Path(batch_file).read_text().splitlines() if l.strip()]
    OUTDIR.mkdir(parents=True, exist_ok=True)

    index = ["# Pilot batch 01 -- enriched gene fiches", "",
             "Source: MTBC0 ancestral genome (Harrison et al. 2024) + H37Rv legacy "
             "+ Pfam (hmmscan) + ESM Atlas + curated literature.", "",
             "| Gene | MTBC0 | Legacy | Revised | Pfam | Verdict | Conf. |",
             "|---|---|---|---|---|---|---|"]
    for rv in rvs:
        row = xref.get(rv)
        esm = json.loads((ESM_DIR / f"{rv}.json").read_text()) if (ESM_DIR / f"{rv}.json").exists() else {}
        pfam = pfam_all.get(rv, [])
        cur = curation.get(rv, {})
        (OUTDIR / f"{rv}.md").write_text(fiche_md(rv, row, esm, pfam, cur))
        revised = cur.get("function_revised", "-")
        revised_short = (revised[:48] + "...") if len(revised) > 51 else revised
        pfam_short = "/".join(h["pfam_name"] for h in pfam) or "-"
        index.append(f"| [{rv}]({rv}.md) | {row['mtbc0']} | {row['product_h37rv']} | "
                     f"{revised_short} | {pfam_short} | {VERDICT_BADGE.get(cur.get('verdict',''),'?')} | "
                     f"{cur.get('confidence','?')} |")
        print(f"  wrote fiches/{rv}.md  [{cur.get('verdict','?')}/{cur.get('confidence','?')}] "
              f"pfam={pfam_short}")
    (OUTDIR / "INDEX.md").write_text("\n".join(index) + "\n")
    print(f"\nWrote {OUTDIR / 'INDEX.md'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/pilot_batch_01.txt")
