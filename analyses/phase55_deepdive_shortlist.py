#!/usr/bin/env python3
"""phase55_deepdive_shortlist.py -- P8 pivot : shortlist des dark méritant un deep-dive dédié.

READ-ONLY. La leçon des batchs 1-2 : la requalification de masse ne marche pas (0 flip / 10). La vraie
valeur pour les dark RÉELS non assignables = identifier ceux qui méritent une ÉTUDE DÉDIÉE type projet
Rv3222c : gène réel (conservé), IMPORTANT (essentiel et/ou hautement vulnérable CRISPRi), TRACTABLE
(bien replié = structure exploitable), et INCONNU (verdict dark). Ce sont des cibles thérapeutiques
candidates + des questions biologiques ouvertes. Classe et priorise ; n'écrit qu'un rapport.
Run: python analyses/phase55_deepdive_shortlist.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "résultats" / "phase55_deepdive"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        if d.get("verdict") != "dark":
            continue
        ess = bool((d.get("essentiality") or {}).get("essential"))
        vi = (d.get("vulnerability") or {}).get("vi")
        vuln = vi is not None and vi <= -5
        plaf = (d.get("plddt_af") or {}).get("mean_plddt") or 0
        ple = (d.get("plddt") or {}).get("mean_plddt") or 0
        wellfolded = max(plaf, ple) >= 80
        pdb = bool((d.get("pdb") or {}).get("n_structures"))
        ms = bool((d.get("proteomics") or {}).get("detected"))
        cons = d.get("conservation") or {}
        snp = cons.get("snp_sites")
        conserved = (snp is not None and snp <= 20)
        seq = d.get("protein_mtbc0") or ""
        L = len(seq)
        # score de priorité deep-dive : importance (essentiel/vulnérable) + tractabilité (structure) + réalité
        importance = (3 if ess else 0) + (3 if vuln else 0) + (1 if (vi is not None and vi <= -3) else 0)
        tractable = (2 if (wellfolded or pdb) else 0) + (1 if ms else 0) + (1 if conserved else 0)
        if importance == 0:
            continue  # deep-dive réservé aux gènes IMPORTANTS (essentiel/vulnérable)
        score = importance + tractable
        rows.append({"rv": d["rv"], "score": score, "ess": ess, "vi": vi, "vuln": vuln,
                     "plddt": round(max(plaf, ple)), "pdb": pdb, "ms": ms, "snp": snp, "len": L,
                     "func": (d.get("function_revised") or d.get("product_h37rv") or "")[:55]})
    rows.sort(key=lambda r: (-r["score"], r["vi"] if r["vi"] is not None else 0))

    lines = ["# Shortlist deep-dive — dark IMPORTANTS et tractables (P8 pivot)", "",
             "Gènes `verdict=dark` réels, essentiels et/ou hautement vulnérables (CRISPRi VI≤-5), bien repliés",
             "(structure exploitable) : cibles thérapeutiques candidates + questions biologiques ouvertes, chacune",
             "méritant une étude dédiée type projet `mtbc/Rv3222c/` (HHpred profil-profil, AF-Multimer, littérature).",
             "NON une requalification.", "",
             f"**{len(rows)} cibles** (score = importance[essentiel+3/vulnérable+3] + tractabilité[structure+2/MS+1/conservé+1]).", "",
             "| score | gène | essentiel | VI | pLDDT | PDB | MS | snp_sites | len | fonction actuelle |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['score']} | {r['rv']} | {'oui' if r['ess'] else '-'} | "
                     f"{r['vi'] if r['vi'] is not None else '?'} | {r['plddt']} | {'oui' if r['pdb'] else '-'} | "
                     f"{'oui' if r['ms'] else '-'} | {r['snp']} | {r['len']} | {r['func']} |")
    (OUT / "deepdive_shortlist.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Deep-dive shortlist : {len(rows)} dark importants (essentiel/vulnérable). Top :")
    for r in rows[:12]:
        tags = ",".join(t for t, ok in [("ess", r["ess"]), ("vuln", r["vuln"]), ("struct", r["plddt"] >= 80 or r["pdb"])] if ok)
        print(f"  {r['rv']:<9} score={r['score']} VI={r['vi']} [{tags}] :: {r['func']}")
    print(f"Rapport : {OUT/'deepdive_shortlist.md'}")


if __name__ == "__main__":
    main()
