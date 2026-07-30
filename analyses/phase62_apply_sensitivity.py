#!/usr/bin/env python3
"""phase62_apply_sensitivity.py -- P11.2c : appliquer l'audit de sensibilité à la couche outgroup.ntm.

Le test de sensibilité (phase61, tblastn relâché) a montré que certains « MTBC-specific » ont en fait un
homologue faible/divergent raté par le seuil strict. On INSCRIT le verdict de sensibilité dans la fiche
(outgroup.ntm.sensitivity) pour que l'atlas ne sur-vende pas : un gène présent-divergent n'est PAS une vraie
innovation spécifique-MTBC. Écrit EN PLACE. (Ne touche que les 17 candidats durcis audités.)
Run: python analyses/phase62_apply_sensitivity.py
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
SENS = ROOT / "résultats" / "phase61" / "sensitivity.tsv"


def main():
    verdicts = {}
    for l in open(SENS):
        if l.startswith("rv"):
            continue
        p = l.rstrip("\n").split("\t")
        verdicts[p[0]] = {"verdict": p[1], "detail": p[2]}
    n = 0
    for rv, v in verdicts.items():
        f = GENES / f"{rv}.json"
        if not f.exists():
            continue
        d = json.load(open(f))
        og = d.get("outgroup") or {}
        ntm = og.get("ntm") or {}
        ntm["sensitivity"] = v["verdict"]
        ntm["sensitivity_detail"] = v["detail"]
        if v["verdict"] != "robust-MTBC-specific":
            # rétrograder la lecture : présent-divergent, pas une innovation
            ntm["interpretation"] = ("SENSITIVITY DOWNGRADE: a divergent homolog is detectable at relaxed thresholds "
                                     f"({v['detail']}); NOT a robust MTBC-specific innovation — present but divergent. "
                                     "The strict tblastn absence was a coverage/identity-threshold artefact.")
        og["ntm"] = ntm
        d["outgroup"] = og
        f.write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
    robust = [rv for rv, v in verdicts.items() if v["verdict"] == "robust-MTBC-specific"]
    print(f"sensibilité inscrite sur {n} fiches. {len(robust)} robustes, {len(verdicts)-len(robust)} rétrogradés.")
    print("Robustes :", ", ".join(robust))


if __name__ == "__main__":
    main()
