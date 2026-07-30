#!/usr/bin/env python3
"""phase48_rv3222c_hhpred.py -- P7.6a-curation : réaction aux hits HHpred de Rv3222c.

Hit HHpred (path web, CG) : top hit = DUF7993 (Pfam PF25956) à **Prob 9,5 % / E-value 1300** = BRUIT TOTAL.
CORRECTION (les probabilités reçues invalident l'interprétation initiale) : ce n'est PAS une assignation de
famille — à 9,5 %/E=1300, DUF7993 est du bruit, Rv3222c ne matche RIEN de significatif. C'est un VRAI NO-HIT
HHpred, cohérent avec la séquence basse complexité / désordonnée. On NE peut PAS dire qu'il « appartient à
DUF7993 ». La fonction repose ENTIÈREMENT sur le contexte SigH-RshA (P7.6a). Leçon : sans la probabilité, ne
pas conclure ; une DUF au TOP de la liste peut n'être que le moins mauvais des non-hits.
Run: python analyses/phase48_rv3222c_hhpred.py
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
F = ROOT / "site" / "content" / "genes" / "Rv3222c.json"
HHPRED_REF = {"authors": "Zimmermann L, Stephens A, Nam SZ, et al.", "year": 2018,
              "title": "A Completely Reimplemented MPI Bioinformatics Toolkit with a New HHpred Server at its Core",
              "journal": "Journal of Molecular Biology", "doi": "10.1016/j.jmb.2017.12.007"}
ADD = (" HHpred profile-profile gives NO confident hit (best is DUF7993 / Pfam PF25956 at only 9.5% probability, "
       "E-value 1300 = noise), consistent with the low-complexity, poorly-ordered sequence; the SigH-RshA "
       "contextual hypothesis remains the sole functional handle.")
# retire une éventuelle mention erronée d'une passe antérieure
STALE = (" HHpred profile-profile: the top hit is the unknown-function family DUF7993 (Pfam PF25956) — this "
         "places Rv3222c in a conserved but functionally uncharacterised family and does NOT resolve its "
         "molecular role (consistent with the low-complexity sequence). The SigH-RshA contextual hypothesis "
         "remains the best functional handle.")


def main():
    d = json.load(open(F))
    fr = (d.get("function_revised") or "").replace(STALE, "")
    if "NO confident hit" not in fr:
        d["function_revised"] = fr.rstrip() + ADD
    else:
        d["function_revised"] = fr
    d["hhpred"] = {"top_hits": ["PF25956 DUF7993 (Family of unknown function) 9.5%, E=1300"],
                   "source": "HHpred (MPI Toolkit web)", "no_confident_hit": True,
                   "note": "No confident hit: best is DUF7993 at 9.5% / E=1300 (noise). Confirms the low-complexity / "
                           "disordered nature; function rests entirely on the SigH-RshA context."}
    refs = d.get("references") or []
    if HHPRED_REF["doi"] not in {r.get("doi") for r in refs}:
        refs.append(HHPRED_REF)
    d["references"] = refs
    F.write_text(json.dumps(d, indent=2, ensure_ascii=False))
    print("== phase48 : Rv3222c — hit HHpred DUF7993 enregistré (fonction non résolue, contexte SigH-RshA gardé) ==")


if __name__ == "__main__":
    main()
