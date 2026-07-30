#!/usr/bin/env python3
"""phase33_literature_curation.py -- requalification des dark « RefSeq-behind » par littérature (P7.2).

Plongée littérature (tbmonitor) sur les 40 hypothétiques encore dark : plusieurs sont en fait
DÉJÀ caractérisés dans un article, mais restés « hypothetical » côté RefSeq/atlas (cas
« RefSeq-behind »). Ce script porte la HAND-CURATION traçable (`auto:false` + références) qui
les sort du dark. Requalification par lecture d'article, sans bioinfo de novo.

Sortie : résultats/phase33_literature_curation/curation.json + fusion (function_revised, verdict,
confidence, gene, références) dans les fiches + enregistrement phase4. Stdlib.
Run: python analyses/phase33_literature_curation.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "phase33_literature_curation"
GENES = ROOT / "site" / "content" / "genes"

R = {
    "yuan2013":  {"authors": "Yuan X, Chen L, Deng X, et al.", "year": 2013, "title": "Characterization of Rv0394c gene encoding hyaluronidase and chondrosulfatase from Mycobacterium tuberculosis", "journal": "Tuberculosis", "doi": "10.1016/j.tube.2013.01.004"},
    "grigg2023": {"authors": "Grigg JC, Copp JN, Krin JM, et al.", "year": 2023, "title": "Deciphering the biosynthesis of a novel lipid in Mycobacterium tuberculosis expands the known roles of the nitroreductase superfamily", "journal": "Journal of Biological Chemistry", "doi": "10.1016/j.jbc.2023.104924"},
    "han2026":   {"authors": "Han X, Arrowsmith TJ, Kardaszewska S, et al.", "year": 2026, "title": "Ribosomal RNA cleavage by the previously unidentified RelS-RelI toxin-antitoxin system controls growth of Mycobacterium tuberculosis", "journal": "Nucleic Acids Research", "doi": "10.1093/nar/gkag571"},
    "arora2020": {"authors": "Arora SK, Alam A, Naqvi N, et al.", "year": 2020, "title": "Immunodominant Mycobacterium tuberculosis Protein Rv1507A Elicits Th1 Response and Modulates Host Macrophage Effector Functions", "journal": "Frontiers in Immunology", "doi": "10.3389/fimmu.2020.01199"},
    "naqvi2025": {"authors": "Naqvi N, Alam A, Shariq M, et al.", "year": 2025, "title": "Leveraging Mycobacterium tuberculosis Rv1507A Protein for Improved Immunogenicity in Mycobacterium bovis BCG Vaccine", "journal": "Immunology", "doi": "10.1111/imm.70063"},
    "gideon2012": {"authors": "Gideon HP, Wilkinson KA, Rustad TR, et al.", "year": 2012, "title": "Bioinformatic and empirical analysis of novel hypoxia-inducible targets of the human antituberculosis T cell response", "journal": "Journal of Immunology", "doi": "10.4049/jimmunol.1202281"},
    "stupar2024": {"authors": "Stupar M, Tan L, Kerr ED, et al.", "year": 2024, "title": "TcrXY is an acid-sensing two-component transcriptional regulator of Mycobacterium tuberculosis required for persistent infection", "journal": "Nature Communications", "doi": "10.1038/s41467-024-45343-7"},
}

# curations expertes (lecture d'abstract tbmonitor)
CUR = {
    "Rv0394c": {"function_revised": "Hyaluronidase / chondroitin-sulfatase; lets M. tuberculosis use host hyaluronan as an alternative carbon source (Yuan 2013).",
                "verdict": "requalified", "confidence": "high", "refs": ["yuan2013"]},
    "Rv2336":  {"gene": "tyzA",
                "function_revised": "Amino-acid N-acyltransferase TyzA of the tyz cluster; N-acylates L-Tyr/L-Phe with lauroyl-CoA in acyl-oxazolone (tyrazolone) biosynthesis (Grigg 2023).",
                "verdict": "requalified", "confidence": "high", "refs": ["grigg2023"]},
    "Rv2663":  {"function_revised": "Ribosome-dependent toxin of the RelS-RelI toxin-antitoxin system; cleaves ribosomal RNA to arrest growth (Han 2026).",
                "verdict": "requalified", "confidence": "high", "refs": ["han2026"]},
    "Rv1507A": {"function_revised": "Immunodominant M. tuberculosis-specific antigen; elicits a Th1 T-cell response and modulates macrophage effector functions (BCG-boost vaccine candidate). Molecular biochemical function still undefined.",
                "verdict": "requalified", "confidence": "medium", "refs": ["arora2020", "naqvi2025"]},
    "Rv1954c": {"function_revised": "Hypoxia-induced protein and empirically validated human T-cell antigen of the anti-TB response (Gideon 2012). Molecular function still undefined.",
                "verdict": "requalified", "confidence": "medium", "refs": ["gideon2012"]},
    "Rv3705A": {"function_revised": "Acid-stress / persistence-associated protein; member of the TcrXY (acid-sensing two-component) regulon required for persistent infection (Stupar 2024). Precise molecular role undefined.",
                "verdict": "family_assigned", "confidence": "low", "refs": ["stupar2024"]},
}


def main():
    print("== phase33 : requalification RefSeq-behind par littérature (P7.2) ==")
    out = {rv: {**c, "references": [R[k] for k in c["refs"]]} for rv, c in CUR.items()}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "curation.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))

    n = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        c = CUR.get(d["rv"])
        if not c:
            continue
        if c.get("gene"):
            d["gene"] = c["gene"]
        d["function_revised"] = c["function_revised"]
        d["verdict"] = c["verdict"]
        d["confidence"] = c["confidence"]
        d["auto"] = False
        d["needs_review"] = False
        refs = d.get("references") or []
        have = {r.get("doi") for r in refs}
        for k in c["refs"]:
            if R[k]["doi"] not in have:
                refs.append(R[k])
        d["references"] = refs
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
        print(f"  {d['rv']} [{c['verdict']}/{c['confidence']}] {c.get('gene','')}: {c['function_revised'][:60]}")
    print(f"Fusionné dans {n} fiches (sortent du dark)")


if __name__ == "__main__":
    main()
