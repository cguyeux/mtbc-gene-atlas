#!/usr/bin/env python3
"""phase42_operon_context_curation2.py -- P7.10, 2e lot (suite de phase41).

Instruit un 2e lot du gisement `operon_context_leads.tsv` (98 dark co-transcrits avec voisin nommé).
Même filtre strict que phase41 : co-transcription réelle + voisin de VOIE nommée + localisation
cohérente + purifiante ; « co-transcrit avec X » ≠ « est X » (candidate <voie>-associated, low).
ÉCARTÉS de ce lot : contextes mobiles/phage (obscurité attendue), loci TA (garde-fou anti-sur-appel),
voisins non caractérisés, et Rv3355c (membranaire co-transcrit avec folD soluble = localisation
incohérente).

Cas particulier Rv3527 : function_revised DÉJÀ raffinée (ferrédoxine/cholestérol) mais verdict resté
`dark` par incohérence, et le voisin kshA (hydroxylase du cholestérol) la corrobore → on met à jour le
VERDICT sans réécrire la fonction (keep_function).
Run: python analyses/phase42_operon_context_curation2.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"

CURATIONS = {
    "Rv1444c": {"operon": "Rv1444c-devB-opcA", "pathway": "oxidative pentose-phosphate pathway",
        "function_revised": (
            "Conserved protein co-transcribed with devB (6-phosphogluconolactonase) and opcA (the OpcA "
            "glucose-6-phosphate-dehydrogenase assembly/activator) in the oxidative pentose-phosphate operon, "
            "under purifying selection. Contextual candidate: associated with the oxidative pentose-phosphate "
            "pathway. Its molecular role is not established.")},
    "Rv2481c": {"operon": "Rv2481c-plsB2-plsC", "pathway": "phospholipid biosynthesis",
        "function_revised": (
            "Conserved protein co-transcribed with plsB2 (glycerol-3-phosphate acyltransferase) and plsC "
            "(1-acylglycerol-3-phosphate acyltransferase) in a phospholipid-biosynthesis operon, under purifying "
            "selection. Contextual candidate: associated with phospholipid / phosphatidic-acid biosynthesis. Its "
            "molecular role is not established.")},
    "Rv0544c": {"operon": "Rv0544c-pitA", "pathway": "phosphate transport",
        "function_revised": (
            "Membrane protein (2 predicted TM helices) co-transcribed with pitA (a low-affinity inorganic "
            "phosphate transporter) under purifying selection. Contextual candidate: associated with phosphate "
            "transport / homeostasis. Its molecular role is not established.")},
    "Rv2091c": {"operon": "Rv2091c-helY-tatC-tatA", "pathway": "Tat protein secretion",
        "function_revised": (
            "Membrane protein (1 predicted TM helix, DUF4333) co-transcribed with tatC and tatA — the twin-arginine "
            "translocation (Tat) protein-secretion system — under purifying selection. Contextual candidate: "
            "associated with the Tat secretion system. Its molecular role is not established.")},
    "Rv2929": {"operon": "tesA-Rv2929", "pathway": "lipid / thioesterase-linked metabolism",
        "function_revised": (
            "Secreted protein co-transcribed with tesA (a thioesterase involved in lipid / PDIM-related metabolism) "
            "under purifying selection. Contextual candidate: associated with lipid / thioesterase-linked "
            "metabolism. Its molecular role is not established.")},
    "Rv3527": {"operon": "kshA-Rv3527", "pathway": "cholesterol catabolism", "keep_function": True},
}


def main():
    print("== phase42 : requalification par contexte d'opéron, 2e lot (P7.10) ==")
    n = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        c = CURATIONS.get(d["rv"])
        if not c:
            continue
        d["verdict"] = "family_assigned"
        d["confidence"] = "low"
        if not c.get("keep_function"):
            d["function_revised"] = c["function_revised"]
        d["auto"] = False
        d["needs_review"] = False
        d["context_lead"] = {"basis": f"co-transcription in operon {c['operon']} + purifying selection",
                             "operon": c["operon"], "pathway": c["pathway"],
                             "source": "operon-context curation (P7.10, 2nd batch)"}
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
        kept = " (verdict fixed, function kept)" if c.get("keep_function") else ""
        print(f"  {d['rv']} -> family_assigned/low [{c['operon']}] :: {c['pathway']}{kept}")
    print(f"{n} fiche(s) requalifiée(s) (2e lot).")


if __name__ == "__main__":
    main()
