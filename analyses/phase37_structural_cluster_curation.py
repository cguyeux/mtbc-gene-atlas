#!/usr/bin/env python3
"""phase37_structural_cluster_curation.py -- P7.3a + P7.3b.

Requalifie des gènes dark à partir du CLUSTERING STRUCTURAL all-vs-all (phase36) : transfert de
fonction depuis un membre CARACTÉRISÉ du même cluster de repli (guilt-by-structure), hand-review.

P7.3a : 3 dark de la famille DUF732 (cluster 1), même repli que les effecteurs sécrétés
        caractérisés Rv1804c / Rv3354 -> family_assigned (rôle d'effecteur plausible, non figé).
P7.3b : Rv0122, repli SIGNIFICATIF (Foldseek E=2.4e-4, TM 0.77) de la toxine Rv2663 (RelS-RelI),
        antitoxine adjacente Rv0123 -> family_assigned/low (candidat toxine TA, activité non prouvée).

Garde-fous appliqués : repli != fonction (formulation « famille / candidat », pas d'activité figée) ;
cross-check littérature fait (Rv0122/Rv0123 non décrits comme TA -> prédiction neuve) ; pLDDT du modèle
signalé. Couche `struct_cluster` ; hand-curation `auto:false` + réfs.
Run: python analyses/phase37_structural_cluster_curation.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"

FOLDSEEK_REF = {"authors": "van Kempen M, Kim SS, Tumescheit C, et al.", "year": 2024,
                "title": "Fast and accurate protein structure search with Foldseek",
                "journal": "Nature Biotechnology", "doi": "10.1038/s41587-023-01773-0"}
RELS_REF = {"authors": "Han X, Arrowsmith TJ, Kardaszewska S, et al.", "year": 2026,
            "title": "Ribosomal RNA cleavage by the previously unidentified RelS-RelI toxin-antitoxin system controls growth of Mycobacterium tuberculosis",
            "journal": "Nucleic Acids Research", "doi": "10.1093/nar/gkag571"}

_DUF732 = ("Secreted protein of the mycobacterial DUF732 family. All-vs-all structural clustering "
           "(Foldseek on ESMFold models, phase36) places it in a tight fold family (median TM-score 0.67) "
           "that includes the characterised host-directed effectors Rv1804c (sequesters host IkappaB-alpha) "
           "and Rv3354 (targets the host COP9 signalosome), plus the secreted antigen Rv1271c. It shares the "
           "family's N-terminal signal peptide (DeepTMHMM) and the ~110-aa secreted fold. A host-interaction / "
           "secreted-effector role is therefore plausible, but the specific molecular target of THIS member is "
           "not experimentally established (the DUF732 family is functionally heterogeneous).")

CURATIONS = {
    "Rv1291c": {"verdict": "family_assigned", "confidence": "medium", "function_revised": _DUF732,
                "cluster": "phase36 cluster 1 (DUF732)", "refs": [FOLDSEEK_REF]},
    "Rv1810":  {"verdict": "family_assigned", "confidence": "medium", "function_revised": _DUF732,
                "cluster": "phase36 cluster 1 (DUF732)", "refs": [FOLDSEEK_REF]},
    "Rv3067":  {"verdict": "family_assigned", "confidence": "medium", "function_revised": _DUF732,
                "cluster": "phase36 cluster 1 (DUF732)", "refs": [FOLDSEEK_REF]},
    "Rv0122":  {"verdict": "family_assigned", "confidence": "low",
                "function_revised": (
                    "Candidate toxin of a type II toxin-antitoxin module. All-vs-all structural clustering "
                    "(Foldseek on ESMFold models, phase36) shows a SIGNIFICANT fold match (E=2.4e-4, TM-score 0.77) "
                    "to Rv2663, the ribosome-dependent ribonuclease toxin of the RelS-RelI system (Han 2026). The "
                    "immediately adjacent, overlapping gene Rv0123 (intergenic gap -4 bp) is a predicted RHH/copG-type "
                    "DNA-binding antitoxin, so Rv0122-Rv0123 is a plausible cognate toxin-antitoxin pair. No prior "
                    "literature describes this locus as a TA system (novel structural prediction). CAVEAT: a shared "
                    "toxin fold is not proof of ribonuclease activity, and the ESMFold model of Rv0122 is of moderate "
                    "confidence (pLDDT 59); toxin activity and TA function remain to be validated experimentally."),
                "cluster": "phase36 cluster 7 (RelS-RelI-like toxin fold)", "refs": [FOLDSEEK_REF, RELS_REF]},
}


def main():
    print("== phase37 : requalification par clustering structural (P7.3a/b) ==")
    n = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        c = CURATIONS.get(d["rv"])
        if not c:
            continue
        d["verdict"] = c["verdict"]
        d["confidence"] = c["confidence"]
        d["function_revised"] = c["function_revised"]
        d["auto"] = False
        d["needs_review"] = False
        d["struct_cluster"] = {"cluster": c["cluster"], "method": "Foldseek all-vs-all on ESMFold dark models (phase36)",
                               "source": "structural clustering (guilt-by-structure)"}
        refs = d.get("references") or []
        have = {r.get("doi") for r in refs}
        for r in c["refs"]:
            if r["doi"] not in have:
                refs.append(r)
        d["references"] = refs
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
        print(f"  {d['rv']} -> {c['verdict']}/{c['confidence']} [{c['cluster']}]")
    print(f"{n} fiche(s) requalifiée(s) par clustering structural.")


if __name__ == "__main__":
    main()
