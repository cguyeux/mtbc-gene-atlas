#!/usr/bin/env python3
"""phase47_rv3222c_deepdive.py -- P7.6a : deep-dive + requalification contextuelle de Rv3222c.

Rv3222c = le gène le PLUS vulnérable du protéome hypothétique (CRISPRi VI -15.5) et pourtant DARK.
Deep-dive multi-couches + littérature :
- essentiel (Tn-seq ES), hautement vulnérable (VI -15.5), fortement conservé (purifiante, pN/pS 0.07) ;
- encodé DANS l'opéron sigH : chevauche (-4 bp de part et d'autre) le facteur sigma alternatif SigH
  (Rv3223c) et son anti-sigma cognate RshA (Rv3221A) ; activé transcriptionnellement par SigH ;
  partenaire STRING #1 = RshA (score 979, context-driven) ;
- séquence BASSE COMPLEXITÉ (183 aa, riche G/A/V, pLDDT 40, aucun hit Pfam/Foldseek) → peu ordonnée ;
- NON caractérisé en littérature (WebSearch + tbmonitor = 0 → nouveauté, pas RefSeq-behind).

Requalification par CONTEXTE (logique P7.10, co-transcription avec gène nommé) : dark -> family_assigned/low,
« candidate essential component of the SigH-RshA stress-response module », fonction non figée (garde-fou :
co-transcription != fonction ; basse complexité limite l'inférence). Hand-off HHpred noté (dernière carte,
mais la basse complexité peut ne rien donner). Run: python analyses/phase47_rv3222c_deepdive.py
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
F = ROOT / "site" / "content" / "genes" / "Rv3222c.json"

NEW_FUNC = (
    "Essential (Tn-seq ES) and HIGHLY VULNERABLE (CRISPRi vulnerability index -15.5, the most vulnerable of the "
    "atlas hypotheticals) conserved protein (strong purifying selection) encoded WITHIN the sigH stress-response "
    "operon: it overlaps (-4 bp on both sides) the alternative sigma factor SigH (Rv3223c) and its cognate "
    "anti-sigma factor RshA (Rv3221A), is transcriptionally activated by SigH, and its top STRING partner is RshA "
    "(score 979). This places it as a candidate essential component of the SigH-RshA oxidative/heat-stress "
    "regulatory module (a co-regulated factor or a modulator of the sigma / anti-sigma switch). Its molecular "
    "function is not established; the low-complexity, poorly-ordered sequence (183 aa, G/A/V-rich, pLDDT 40, no "
    "Pfam or Foldseek hit) limits structural inference. NOT previously characterised in the literature "
    "(WebSearch + tbmonitor negative) - a novel, essential, highly vulnerable, uncharacterised drug-target candidate.")

NOTE = ("2026-07-05 (P7.6a deep-dive): requalified from dark by genomic/regulatory context (sigH-rshA operon + "
        "SigH regulation + STRING-RshA), corroborated by essentiality + CRISPRi vulnerability. Function candidate, "
        "not established. HHpred web still worth trying (low-complexity sequence may yield nothing).")


def main():
    d = json.load(open(F))
    d["verdict"] = "family_assigned"
    d["confidence"] = "low"
    d["function_revised"] = NEW_FUNC
    d["auto"] = False
    d["needs_review"] = False
    d["curation_note"] = NOTE
    d["context_lead"] = {"basis": "co-transcription in the sigH-rshA operon (overlaps SigH and RshA) + SigH regulation "
                         "+ STRING-RshA + essentiality + CRISPRi vulnerability",
                         "operon": "rshA-Rv3222c-sigH", "pathway": "SigH-RshA oxidative/heat-stress response module",
                         "source": "deep-dive contextual curation (P7.6a)"}
    F.write_text(json.dumps(d, indent=2, ensure_ascii=False))
    print("== phase47 : Rv3222c requalifié (P7.6a) ==")
    print("  dark -> family_assigned/low :: candidate SigH-RshA stress-module component")


if __name__ == "__main__":
    main()
