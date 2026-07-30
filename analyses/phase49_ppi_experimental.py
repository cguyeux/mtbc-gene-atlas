#!/usr/bin/env python3
"""phase49_ppi_experimental.py -- P7.7 : PPI expérimentales vs prédites (recadré).

CONSTAT : le canal `experimental` de STRING est DÉJÀ en base (1309 gènes, 6724 arêtes fortes, 548
dark/family_assigned à partenaire expérimental caractérisé) et déjà exploité (P7.10b). STRING intègre
BioGRID/IntAct comme sources → réimporter en brut = redondant ET réintroduit le bruit du grand
interactome B2H de Wang 2010. IntAct confirmé accessible (interactions M.tb curées, ex. sigA-whiB3)
mais largement dans STRING-experimental.

RÉSIDU utile = 8 dark à ancre expérimentale forte (exp≥400). GARDE-FOU CRUCIAL (leçon P7.7) : une PPI
expérimentale (surtout B2H) qui CONTREDIT le contexte (opéron / localisation) est un FAUX POSITIF.
Croisement des 8 :
- Rv0863 ~ctaD (987) : mais opéron Rv0863-moaC2-mog-moaE2 (molybdoptérine) → l'interaction ctaD (hub)
  CONTREDIT l'opéron = artefact B2H ; requalifier par l'OPÉRON (voie Moco), PAS par la PPI. ← seul gain
- Rv0943c ~atpA : SÉCRÉTÉ (SP) vs ATP synthase interne = compartiments incompatibles → écarté.
- Rv1363c ~pks2 : conservation relâchée + PPI seule → trop faible, écarté.
- Rv1158c/Rv2472/Rv3103c/Rv2077c/Rv0614 : localisation / contexte TA / PE-PGRS répétitif / tRNA-synthetase
  (hub non spécifique) → écartés.
Résultat : 1 requalification (Rv0863 par opéron Moco) + leçon « PPI expérimentale B2H à croiser avec le
contexte ». Run: python analyses/phase49_ppi_experimental.py
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
F = ROOT / "site" / "content" / "genes" / "Rv0863.json"


def main():
    d = json.load(open(F))
    assert d.get("verdict") == "dark"
    d["verdict"] = "family_assigned"
    d["confidence"] = "low"
    d["function_revised"] = (
        "Candidate component of molybdopterin (Moco) cofactor biosynthesis: first gene of the operon "
        "Rv0863-moaC2-mog-moaE2 (co-transcribed with the molybdopterin-synthesis genes moaC2/mog/moaE2), "
        "so contextually associated with Moco biosynthesis. Its precise molecular role is not established. "
        "NOTE: a strong STRING experimental interaction with the cytochrome-oxidase subunit ctaD (score 987) "
        "is discounted as a likely bacterial-two-hybrid false positive, as ctaD is a hub and the interaction "
        "contradicts the molybdopterin operon context (P7.7 guard: an experimental PPI that conflicts with "
        "the operon is treated as an artefact).")
    d["auto"] = False
    d["needs_review"] = False
    d["context_lead"] = {"basis": "co-transcription in the Rv0863-moaC2-mog-moaE2 operon (molybdopterin biosynthesis)",
                         "operon": "Rv0863-moaC2-mog-moaE2", "pathway": "molybdopterin (Moco) cofactor biosynthesis",
                         "source": "operon-context curation via PPI cross-check (P7.7)"}
    F.write_text(json.dumps(d, indent=2, ensure_ascii=False))
    print("== phase49 (P7.7) : Rv0863 requalifié par opéron Moco ; PPI ctaD écartée comme artefact B2H ==")


if __name__ == "__main__":
    main()
