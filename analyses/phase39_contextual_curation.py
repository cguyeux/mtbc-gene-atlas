#!/usr/bin/env python3
"""phase39_contextual_curation.py -- P7.10a : requalification par CONTEXTE (co-transcription).

Instruit les dark à phenotype_lead. RÉSULTAT du filtre strict « fonction vs phénotype » : la plupart
des 35 leads sont « required in vivo » = un PHÉNOTYPE, pas une fonction moléculaire → NON requalifiables
(le phenotype_lead les documente déjà comme hypothèses ; verdict dark maintenu). SEULS passent ceux où
une CO-TRANSCRIPTION avec un gène nommé (opéron réel, pas simple STRING) donne une piste fonctionnelle
cohérente avec la localisation et la conservation :

- Rv1312 : gène TERMINAL de l'opéron F1F0-ATP synthase (Rv1303-atpBEFGHDC-Rv1312, littérature confirmée),
           membranaire, purifiante forte -> candidat composant accessoire de l'oxidative phosphorylation.
- Rv3231c : co-transcrit avec ppk2 (polyphosphate kinase 2), globulaire, purifiante forte -> candidat
            facteur associé au métabolisme du polyphosphate.

Garde-fou : « co-transcrit dans l'opéron X » ≠ « sous-unité/enzyme X » → formulation « candidate
associated component », fonction précise non figée. Couche `context_lead`, hand-review `auto:false`.
Run: python analyses/phase39_contextual_curation.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"

CURATIONS = {
    "Rv1312": {
        "verdict": "family_assigned", "confidence": "low",
        "function_revised": (
            "Membrane protein encoded as the terminal gene of the F1F0-ATP synthase operon "
            "(Rv1303-atpBEFGHDC-Rv1312; operon structure established in the published M. tuberculosis "
            "oxidative-phosphorylation literature), co-transcribed with the atp genes and under strong "
            "purifying selection. The convergent context — co-transcription within the ATP-synthase operon, "
            "predicted membrane localisation, and an in-vivo fitness requirement (Tn-seq) — makes it a "
            "candidate accessory / associated component of the mycobacterial ATP synthase (oxidative "
            "phosphorylation). Whether it is a structural subunit, an assembly factor, or a co-regulated "
            "membrane protein is not established. RefSeq leaves it 'hypothetical protein'."),
        "basis": "co-transcription in the atp operon + membrane localisation + in-vivo fitness + purifying selection",
        "operon": "Rv1303-atpBEFGHDC-Rv1312"},
    "Rv3231c": {
        "verdict": "family_assigned", "confidence": "low",
        "function_revised": (
            "Conserved cytoplasmic protein co-transcribed with ppk2 (polyphosphate kinase 2, Rv3232c) as a "
            "two-gene operon, under strong purifying selection and required for in-vivo fitness. The genomic "
            "context (operon with ppk2) points to a candidate role associated with polyphosphate metabolism, "
            "a persistence- and stress-linked pathway; the specific molecular function is not established."),
        "basis": "co-transcription with ppk2 + purifying selection + in-vivo fitness",
        "operon": "Rv3231c-ppk2"},
}


def main():
    print("== phase39 : requalification contextuelle (P7.10a) ==")
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
        d["context_lead"] = {"basis": c["basis"], "operon": c["operon"],
                             "source": "contextual multi-channel curation (P7.10a): operon co-transcription + phenotype + conservation"}
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
        print(f"  {d['rv']} -> {c['verdict']}/{c['confidence']} [{c['operon']}]")
    print(f"{n} fiche(s) requalifiée(s) par contexte. (33 autres leads = phénotype seul, restent dark.)")


if __name__ == "__main__":
    main()
