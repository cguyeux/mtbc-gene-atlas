#!/usr/bin/env python3
"""phase44_string_context_curation.py -- P7.10b : requalification par CONTEXTE STRING (hors opéron).

Complète P7.10 (opéron) pour les dark SANS voisin d'opéron nommé, via le réseau STRING. Garde-fou
RENFORCÉ vs l'opéron (guilt-by-association plus faible que la co-transcription, KB) : on n'accepte
qu'un partenaire NOMMÉ de fonction connue, lié par un canal FIABLE (experimental OU database >=300,
PAS neighborhood/text-mining/cooccurence seuls) ET combined_no_tm >=700, AVEC un 2e signal concordant
(localisation cohérente OU régulon). Formulation « candidate functional partner of X » / guilt-by-
association, jamais la fonction de X littéralement.

Sur 6 candidats du filtre : 4 retenus ; ÉCARTÉS Rv2472 (partenaire = toxine parE1, contexte TA →
garde-fou anti-sur-appel) et Rv1158c (SP sécrété vs PNPase gpsI cytoplasmique = localisation
incompatible). Rv2468A porte un CAVEAT explicite : le supercomplexe respiratoire bcc:aa3 est un sink
de promiscuité de petites sous-unités TM (KB 2026-06-12) — le signal ici est experimental=997, mais
on cadre « candidate association à valider », pas « sous-unité assignée ».
Run: python analyses/phase44_string_context_curation.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"

STRING_REF = {"authors": "Szklarczyk D, Kirsch R, Koutrouli M, et al.", "year": 2023,
              "title": "The STRING database in 2023",
              "journal": "Nucleic Acids Research", "doi": "10.1093/nar/gkac1000"}

CURATIONS = {
    "Rv2468A": {"partner": "ctaD (cytochrome c oxidase subunit)", "process": "respiratory cytochrome c oxidase",
        "function_revised": (
            "Membrane protein with a strong STRING experimental association (score 997) to ctaD, a subunit of the "
            "cytochrome c (aa3) oxidase of the respiratory bcc:aa3 supercomplex. Contextual candidate: a small "
            "accessory / partner component of the cytochrome c oxidase. CAVEAT: the bcc:aa3 supercomplex carries "
            "many uncharacterised small TM subunits and is a known promiscuity sink, so this is a guilt-by-"
            "association candidate to validate, not an assigned subunit.")},
    "Rv2731": {"partner": "scpA (segregation-and-condensation protein)", "process": "chromosome segregation/condensation",
        "function_revised": (
            "Conserved protein with a strong STRING experimental association (score 788) to scpA (a "
            "segregation-and-condensation protein of the SMC/condensin chromosome-segregation complex), and itself "
            "part of a transcriptional regulon (second concordant channel). Contextual candidate: associated with "
            "chromosome segregation / condensation. Guilt-by-association, not an established function.")},
    "Rv3769": {"partner": "scpA (segregation-and-condensation protein)", "process": "chromosome segregation/condensation",
        "function_revised": (
            "Conserved protein with a strong STRING experimental association (score 788) to scpA (a "
            "segregation-and-condensation protein of the chromosome-segregation machinery). Contextual candidate: "
            "associated with chromosome segregation / condensation. Guilt-by-association, not an established function.")},
    "Rv1004c": {"partner": "ftsY (SRP receptor)", "process": "protein targeting / secretion",
        "function_revised": (
            "Signal-peptide protein with a strong STRING experimental association (score 829) to ftsY, the "
            "signal-recognition-particle receptor of the co-translational protein-targeting/insertion machinery. "
            "Contextual candidate: associated with protein targeting / secretion. Guilt-by-association, not an "
            "established function.")},
}


def main():
    print("== phase44 : requalification par contexte STRING, hors opéron (P7.10b) ==")
    n = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        c = CURATIONS.get(d["rv"])
        if not c:
            continue
        assert d.get("verdict") == "dark", f"{d['rv']} n'est plus dark"
        d["verdict"] = "family_assigned"
        d["confidence"] = "low"
        d["function_revised"] = c["function_revised"]
        d["auto"] = False
        d["needs_review"] = False
        d["context_lead"] = {"basis": f"STRING experimental association to {c['partner']} + concordant channel",
                             "partner": c["partner"], "pathway": c["process"],
                             "source": "STRING guilt-by-association curation (P7.10b)"}
        refs = d.get("references") or []
        if STRING_REF["doi"] not in {r.get("doi") for r in refs}:
            refs.append(STRING_REF)
        d["references"] = refs
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
        print(f"  {d['rv']} -> family_assigned/low :: partner {c['partner']} ({c['process']})")
    print(f"{n} fiche(s) requalifiée(s) par STRING (2 écartées : TA, localisation incohérente).")


if __name__ == "__main__":
    main()
