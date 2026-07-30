#!/usr/bin/env python3
"""phase4_export_json.py -- Consolidate an enriched gene into one stable JSON.

Merges, per gene of a batch, the outputs of phase1/2/2b/2c + the curation
layer into a single self-contained record under site/content/genes/<rv>.json.
This is the artefact the web site ingests (the pipeline produces data; the
site renders it; they are decoupled by this JSON contract).

The site also ingests the full catalogue (data/gene_xref.tsv, 3022 genes)
directly; phase4 only emits the *enriched* genes (curation + ESM + Pfam +
optional Foldseek).

Usage:  python phase4_export_json.py data/pilot_batch_01.txt
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
XREF = Path(os.environ.get("GENE_XREF", ROOT / "data" / "gene_xref.tsv"))
ESM_DIR = ROOT / "résultats" / "phase2_esm"
PFAM_FILE = ROOT / "résultats" / "phase2b_pfam" / "pfam.json"
FOLDSEEK_FILE = ROOT / "résultats" / "phase2c_foldseek" / "foldseek.json"
EGGNOG_FILE = ROOT / "résultats" / "phase2d_eggnog" / "eggnog.json"
UNIPROT_FILE = ROOT / "résultats" / "phase2e_uniprot" / "uniprot.json"
CONSERV_FILE = ROOT / "résultats" / "phase2f_conservation" / "conservation.json"
PLDDT_FILE = ROOT / "résultats" / "phase2c_foldseek" / "plddt.json"
STRING_FILE = ROOT / "résultats" / "phase2h_string" / "string.json"
MCSA_FILE = ROOT / "résultats" / "phase2i_mcsa" / "mcsa.json"
PFAM_TENTATIVE_FILE = ROOT / "résultats" / "phase17_pfam_relaxed" / "pfam_relaxed.json"
LOCALIZATION_FILE = ROOT / "résultats" / "phase18_localization" / "localization.json"
FOLDSEEK_AF_FILE = ROOT / "résultats" / "phase16_alphafold" / "foldseek_af.json"
PLDDT_AF_FILE = ROOT / "résultats" / "phase16_alphafold" / "plddt_af.json"
ESSENTIALITY_FILE = ROOT / "résultats" / "phase14_essentiality" / "essentiality.json"
PROTEOMICS_FILE = ROOT / "résultats" / "phase19_proteomics" / "proteomics.json"
PROTPARAM_FILE = ROOT / "résultats" / "phase20_protparam" / "protparam.json"
RESISTANCE_FILE = ROOT / "résultats" / "phase21_resistance" / "resistance_by_gene.json"
FUNCCAT_FILE = ROOT / "résultats" / "phase22_funccat" / "funccat.json"
ORTHOLOGS_FILE = ROOT / "résultats" / "phase23_orthologs" / "orthologs.json"
MUTPHENO_FILE = ROOT / "résultats" / "phase24_mutant_phenotypes" / "mutant_phenotypes.json"
PDB_FILE = ROOT / "résultats" / "phase25_pdb" / "pdb.json"
OPERON_FILE = ROOT / "résultats" / "phase26_operon" / "operon.json"
REGULATION_FILE = ROOT / "résultats" / "phase27_regulation" / "regulation.json"
MYCO_CONC_FILE = ROOT / "résultats" / "phase28_concordance" / "concordance.json"
RD_FILE = ROOT / "résultats" / "phase29_rd" / "rd.json"
EC_OVERRIDE_FILE = ROOT / "résultats" / "phase31_ec_override" / "ec_override.json"
PHENO_LEAD_FILE = ROOT / "résultats" / "phase32_phenotype_requalification" / "leads.json"
CURATION = ROOT / "data" / "pilot_curation.json"
AUTOCURATION = ROOT / "data" / "autocuration.json"
OUTDIR = ROOT / "site" / "content" / "genes"

SCHEMA_VERSION = 1


def load_xref() -> dict[str, dict]:
    with open(XREF) as fh:
        return {r["rv"]: r for r in csv.DictReader(fh, delimiter="\t")}


def main(batch_file: str) -> None:
    xref = load_xref()
    curation = json.loads(CURATION.read_text())
    # Merge the rule-based autocuration underneath the manual one (manual wins).
    if AUTOCURATION.exists():
        for k, v in json.loads(AUTOCURATION.read_text()).items():
            curation.setdefault(k, v)
    # Hand-curations (P7.2 littérature RefSeq-behind, P7.1 HHpred) — hand-review, WINS over everything.
    # Seules les entrées AVEC verdict entrent dans la curation (les HHpred "no confident hit" n'y vont pas).
    for cur_file in [ROOT / "résultats" / "phase33_literature_curation" / "curation.json",
                     ROOT / "résultats" / "phase34_hhpred_curation" / "hhpred_curation.json"]:
        if cur_file.exists():
            for k, v in json.loads(cur_file.read_text()).items():
                if not v.get("verdict"):
                    continue
                curation[k] = {**curation.get(k, {}), **v, "auto": False, "needs_review": False}
    hhpred_all = json.loads((ROOT / "résultats" / "phase34_hhpred_curation" / "hhpred_curation.json").read_text()) \
        if (ROOT / "résultats" / "phase34_hhpred_curation" / "hhpred_curation.json").exists() else {}
    pfam_all = json.loads(PFAM_FILE.read_text()) if PFAM_FILE.exists() else {}
    foldseek_all = json.loads(FOLDSEEK_FILE.read_text()) if FOLDSEEK_FILE.exists() else {}
    eggnog_all = json.loads(EGGNOG_FILE.read_text()) if EGGNOG_FILE.exists() else {}
    uniprot_all = json.loads(UNIPROT_FILE.read_text()) if UNIPROT_FILE.exists() else {}
    conserv_all = json.loads(CONSERV_FILE.read_text()) if CONSERV_FILE.exists() else {}
    plddt_all = json.loads(PLDDT_FILE.read_text()) if PLDDT_FILE.exists() else {}
    string_all = json.loads(STRING_FILE.read_text()) if STRING_FILE.exists() else {}
    essent_all = json.loads(ESSENTIALITY_FILE.read_text()) if ESSENTIALITY_FILE.exists() else {}
    proteomics_all = json.loads(PROTEOMICS_FILE.read_text()) if PROTEOMICS_FILE.exists() else {}
    protparam_all = json.loads(PROTPARAM_FILE.read_text()) if PROTPARAM_FILE.exists() else {}
    resistance_all = json.loads(RESISTANCE_FILE.read_text()) if RESISTANCE_FILE.exists() else {}
    funccat_all = json.loads(FUNCCAT_FILE.read_text()) if FUNCCAT_FILE.exists() else {}
    orthologs_all = json.loads(ORTHOLOGS_FILE.read_text()) if ORTHOLOGS_FILE.exists() else {}
    mutpheno_all = json.loads(MUTPHENO_FILE.read_text()) if MUTPHENO_FILE.exists() else {}
    pdb_all = json.loads(PDB_FILE.read_text()) if PDB_FILE.exists() else {}
    operon_all = json.loads(OPERON_FILE.read_text()) if OPERON_FILE.exists() else {}
    regulation_all = json.loads(REGULATION_FILE.read_text()) if REGULATION_FILE.exists() else {}
    myco_conc_all = json.loads(MYCO_CONC_FILE.read_text()) if MYCO_CONC_FILE.exists() else {}
    rd_all = json.loads(RD_FILE.read_text()) if RD_FILE.exists() else {}
    ec_override_all = json.loads(EC_OVERRIDE_FILE.read_text()) if EC_OVERRIDE_FILE.exists() else {}
    pheno_lead_all = json.loads(PHENO_LEAD_FILE.read_text()) if PHENO_LEAD_FILE.exists() else {}
    mcsa_all = json.loads(MCSA_FILE.read_text()) if MCSA_FILE.exists() else {}
    pfam_tentative_all = json.loads(PFAM_TENTATIVE_FILE.read_text()) if PFAM_TENTATIVE_FILE.exists() else {}
    localization_all = json.loads(LOCALIZATION_FILE.read_text()) if LOCALIZATION_FILE.exists() else {}
    foldseek_af_all = json.loads(FOLDSEEK_AF_FILE.read_text()) if FOLDSEEK_AF_FILE.exists() else {}
    plddt_af_all = json.loads(PLDDT_AF_FILE.read_text()) if PLDDT_AF_FILE.exists() else {}
    rvs = [l.strip() for l in Path(batch_file).read_text().splitlines() if l.strip()]
    OUTDIR.mkdir(parents=True, exist_ok=True)

    for rv in rvs:
        row = xref.get(rv)
        if not row:
            print(f"[skip] {rv}: not in xref"); continue
        cur = curation.get(rv, {})
        esm_cond = {}
        ef = ESM_DIR / f"{rv}.json"
        if ef.exists():
            esm_cond = json.loads(ef.read_text()).get("condensed", {})
        homologs = esm_cond.get("esm_homologs", [])
        sim_max = max((h["similarity"] for h in homologs if h.get("similarity")), default=None)

        domains = [{"source": "pfam", **h} for h in pfam_all.get(rv, [])]
        struct_hits = foldseek_all.get(rv, [])
        eggnog = eggnog_all.get(rv, {})
        uniprot = uniprot_all.get(rv, {})
        conservation = conserv_all.get(rv, {})
        plddt = plddt_all.get(rv, {})
        string = string_all.get(rv, {})

        record = {
            "schema_version": SCHEMA_VERSION,
            "rv": rv,
            "mtbc0": row["mtbc0"],
            "np_id": row["np_id"],
            "gene": cur.get("gene") or row.get("gene") or "",
            "len_aa": int(row["len_aa"]),
            "strand": row["strand"],
            "start_mtbc0": int(row["start_mtbc0"]),
            "end_mtbc0": int(row["end_mtbc0"]),
            "product_h37rv": row["product_h37rv"],
            "product_mtbc0_pgap": row["product_mtbc0_pgap"],
            "function_revised": cur.get("function_revised", ""),
            "verdict": cur.get("verdict", ""),
            "confidence": cur.get("confidence", ""),
            "auto": cur.get("auto", False),
            "needs_review": cur.get("needs_review", False),
            "protein_mtbc0": row["protein_mtbc0"],
            "domains": domains,
            "struct_hits": struct_hits,
            "eggnog": eggnog,
            "uniprot": uniprot,
            "conservation": conservation,
            "plddt": plddt,
            "string": string,
            "plddt_af": plddt_af_all.get(rv, {}),
            "struct_af": {"hits": foldseek_af_all.get(rv, []),
                          "plddt_gated": bool(plddt_af_all.get(rv, {}).get("mean_plddt", 0) >= 70)},
            "mcsa": mcsa_all.get(rv, {}),
            "pfam_tentative": pfam_tentative_all.get(rv, {}),
            "localization": localization_all.get(rv, {}),
            "essentiality": essent_all.get(rv, {}),
            "proteomics": proteomics_all.get(rv, {}),
            "protparam": protparam_all.get(rv, {}),
            "resistance": resistance_all.get(rv, {}),
            "funccat": funccat_all.get(rv, {}),
            "orthologs": orthologs_all.get(rv, {}),
            "mutant_phenotypes": mutpheno_all.get(rv, {}),
            "pdb": pdb_all.get(rv, {}),
            "genomic_context": operon_all.get(rv, {}),
            "regulation": regulation_all.get(rv, {}),
            "mycobrowser": myco_conc_all.get(rv, {}),
            "rd": rd_all.get(rv, {}),
            "ec_override": ec_override_all.get(rv, {}),
            "phenotype_lead": pheno_lead_all.get(rv, {}),
            "hhpred": {k: hhpred_all[rv][k] for k in ("top_hits", "note", "no_confident_hit", "source")
                       if rv in hhpred_all and k in hhpred_all[rv]},
            "esm": {
                "protein_hash": esm_cond.get("protein_hash"),
                "mean_plddt": esm_cond.get("mean_plddt"),
                "n_homologs": len(homologs),
                "max_similarity": sim_max,
                "top_sae_features": esm_cond.get("top_sae_features", [])[:8],
            },
            "evidence": cur.get("evidence", []),
            "references": cur.get("references", []),
        }
        # surcharge EC curée (P5.11c) : hand-review traçable -> auto=false + référence
        ov = ec_override_all.get(rv)
        if ov:
            record["auto"] = False
            record["needs_review"] = False
            have = {r.get("title") for r in record["references"]}
            record["references"] += [r for r in ov.get("references", []) if r.get("title") not in have]
        (OUTDIR / f"{rv}.json").write_text(json.dumps(record, indent=2, ensure_ascii=False))
        print(f"  wrote site/content/genes/{rv}.json "
              f"[{record['verdict'] or 'uncurated'}] "
              f"{len(domains)} Pfam, {len(struct_hits)} struct hits")
    print(f"\nDone -> {OUTDIR}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/pilot_batch_01.txt")
