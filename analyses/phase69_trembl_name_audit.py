#!/usr/bin/env python3
"""phase69_trembl_name_audit.py -- P16.14 : micro-audit « nom TrEMBL d'enzyme vs verdict=dark ».

Le scan P16.0 a révélé (et j'ai confirmé sur disque) 6 gènes `dark` portant un nom UniProt d'enzyme
spécifique DÉJÀ en base mais `reviewed=False` (TrEMBL automatique, à raison non flippé). Règle d'adjudication
(anti-survente) : ne flipper dark->family_assigned QUE si un signal INDÉPENDANT (eggNOG family/OG, Pfam,
Foldseek+M-CSA, conservation) corrobore le nom. Sinon rester dark en DOCUMENTANT pourquoi.

Résultat de l'adjudication manuelle (dossiers complets lus) :
- Rv0470A : CORROBORÉ (eggNOG COG1575 « MenA family » concorde avec le nom) -> FLIP family_assigned/low.
- Rv2699c, Rv0614, Rv3566A, Rv0787, Rv3190c : non corroborés / contredits / contamination de voisin / repli-seul
  -> restent dark, note d'adjudication écrite.

Édite EN PLACE site/content/genes/*.json (patron des scripts de curation ; NE PAS re-run phase4, destructif).
Run: python analyses/phase69_trembl_name_audit.py
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"

D = "P16.14 TrEMBL-name audit (2026-07-11): "

# rv -> dict de champs à écrire. `flip` porte verdict/confidence/function_revised ; sinon curation_note seule.
CURATIONS = {
    "Rv0470A": {
        "flip": True,
        "verdict": "family_assigned",
        "confidence": "low",
        "function_revised": (
            "Candidate UbiA/MenA-family prenyltransferase (menaquinone-biosynthesis-related). Family handle from "
            "eggNOG COG1575 ('Belongs to the MenA family. Type 1 subfamily'), concordant with the UniProt TrEMBL "
            "name '1,4-dihydroxy-2-naphthoate prenyltransferase'. Family-level assignment only, not experimentally "
            "validated; carries a pseudogene-candidate flag and relaxed intra-MTBC selection, consistent with a "
            "possibly degrading paralogue of the primary menA (Rv0534c)."
        ),
        "curation_note": (
            D + "unreviewed UniProt enzyme name corroborated by an independent signal (eggNOG COG1575, MenA/UbiA "
            "prenyltransferase family) -> dark -> family_assigned (low). Caveat: pseudogene-candidate + relaxed selection."
        ),
    },
    "Rv2699c": {
        "curation_note": (
            D + "UniProt TrEMBL name 'dUTPase' NOT corroborated and contradicted -- eggNOG assigns DUF4193 and Foldseek "
            "shows no dUTPase; the annotated dut (dUTPase) is the near neighbour Rv2697c, so the name is a "
            "mis-propagation. Kept dark. NB a real, important gene: essential (DeJesus ES), highly vulnerable "
            "(VI -6.8), DUF4193, conserved to Actinomycetia -- function still unknown."
        ),
    },
    "Rv0614": {
        "curation_note": (
            D + "UniProt TrEMBL name 'galactose-1-phosphate uridylyltransferase' is name-only -- no independent support "
            "(eggNOG description empty, no Pfam domain, Foldseek best a spurious MATE hit at prob 0.03). Kept dark."
        ),
    },
    "Rv3566A": {
        "curation_note": (
            D + "UniProt TrEMBL name 'arylamine N-acetyltransferase' is mis-propagated from the immediately adjacent, "
            "co-expressed nat (Rv3566c; STRING neighborhood 781 + coexpression 731); no independent support for Rv3566A "
            "itself. Kept dark -- an MTBC-specific small ORF at the nat locus, function unknown."
        ),
    },
    "Rv0787": {
        "curation_note": (
            D + "UniProt TrEMBL name 'metallo-beta-lactamase superfamily protein, putative' is a broad-fold hedge; "
            "Foldseek best is an uncharacterised protein (M. smegmatis Rv0999 ortholog) with only weak PBP3 hits, and "
            "there is no M-CSA active site -- fold does not imply an active enzyme. Kept dark (purifying intra-MTBC "
            "selection = a real gene, function unknown)."
        ),
    },
    "Rv3190c": {
        "curation_note": (
            D + "UniProt TrEMBL name 'restriction endonuclease' NOT corroborated -- Foldseek indicates a "
            "DNA-polymerase-X (pol beta/lambda) nucleotidyltransferase fold, not a restriction endonuclease, and eggNOG "
            "gives only an unannotated archaeal COG. Pseudogene-candidate + MTBC-specific + relaxed selection. Kept dark."
        ),
    },
}


def main() -> None:
    flips = 0
    for rv, c in CURATIONS.items():
        f = GENES / f"{rv}.json"
        d = json.loads(f.read_text())
        d["curation_note"] = c["curation_note"]
        if c.get("flip"):
            d["verdict"] = c["verdict"]
            d["confidence"] = c["confidence"]
            d["function_revised"] = c["function_revised"]
            d["auto"] = False
            d["needs_review"] = False
            flips += 1
            print(f"  FLIP {rv}: dark -> {c['verdict']} ({c['confidence']})")
        else:
            print(f"  KEEP {rv}: dark (adjudication documentée)")
        f.write_text(json.dumps(d, indent=2, ensure_ascii=False))
    print(f"\n{len(CURATIONS)} fiches auditées, {flips} flip(s) (Rv0470A), {len(CURATIONS)-flips} maintenues dark documentées.")
    print("Garde-fou P6.2 d'ingest : recale un verdict<->fonction incohérent à l'ingestion (ici cohérent).")


if __name__ == "__main__":
    main()
