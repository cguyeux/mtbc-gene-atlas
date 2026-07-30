#!/usr/bin/env python3
"""phase38_rv3032a_correction.py -- P7.3c : rétrograder l'annotation sur-vendue de Rv3032A.

Drapeau levé pendant le parcours des clusters (phase36) : Rv3032A était annoté « membrane
metalloprotease-like (YxkI-like) » alors qu'AUCUN signal ne le soutient — pas de motif catalytique
HExxH, prédit GLOBULAIRE (0 hélice TM, donc pas « membrane »), eggNOG/UniProt/Pfam muets,
struct_af non significatif (pointe Spt16 chromatine). Garde-fou du projet : repli ≠ enzyme active.

Ce qui RESTE vrai (préservé) : gène réel sous forte sélection purifiante (pN/pS=0, 145 209 souches) ;
repli défini partagé avec Rv1914c (petite famille structurale intra-MTBC, Foldseek TM 0,78) ;
ressemblance LOINTAINE au repli YxkI, sans valeur fonctionnelle. Fonction réécrite, metalloprotease
retirée. Hand-review `auto:false`. Run: python analyses/phase38_rv3032a_correction.py
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
F = ROOT / "site" / "content" / "genes" / "Rv3032A.json"

NEW_FUNC = (
    "Conserved protein under strong purifying selection (pN/pS=0 over 145,209 strains), hence a genuine "
    "functional gene, adopting a defined fold shared with Rv1914c (a small intra-MTBC structural family, "
    "Foldseek TM-score 0.78). The fold is only remotely reminiscent of the YxkI family; Rv3032A LACKS the "
    "HExxH zinc-metalloprotease catalytic motif and is predicted cytoplasmic (globular, no transmembrane "
    "helix). The previously proposed 'membrane metalloprotease-like' assignment is therefore NOT supported "
    "and is withdrawn. Molecular function remains undetermined.")

CORR_NOTE = ("2026-07-04 (P7.3c): over-call withdrawn — 'membrane metalloprotease-like (YxkI-like)' was not "
             "supported (no HExxH motif; predicted cytoplasmic; eggNOG/UniProt/Pfam silent; struct_af non-significant). "
             "Kept as fold-only / function-unknown; gene confirmed real by strong purifying selection.")


def main():
    d = json.load(open(F))
    old = d.get("function_revised")
    d["function_revised"] = NEW_FUNC
    d["verdict"] = "family_assigned"      # repli défini conservé, mais fonction indéterminée
    d["confidence"] = "low"
    d["auto"] = False
    d["needs_review"] = False
    d["curation_note"] = CORR_NOTE
    F.write_text(json.dumps(d, indent=2, ensure_ascii=False))
    print("== phase38 : Rv3032A corrigé (P7.3c) ==")
    print("  AVANT:", (old or "")[:80])
    print("  APRÈS:", NEW_FUNC[:80], "...")


if __name__ == "__main__":
    main()
