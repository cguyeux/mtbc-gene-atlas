#!/usr/bin/env python3
"""phase52_druggable_structure.py -- P7.6d : croiser vulnérabilité CRISPRi × structure.

Une cible thérapeutique est exploitable pour du structure-based drug design (docking) SI elle est à la
fois (a) HAUTEMENT VULNÉRABLE au knockdown (CRISPRi VI ≤ seuil) et (b) STRUCTURELLEMENT résolue. Ce
script croise la couche `vulnerability` (P7.6, idéalement étendue à tout le protéome par P7.6c) avec les
couches structurales déjà en base : PDB expérimentale (`pdb`), modèle AlphaFold fiable (`plddt_af`), hit
Foldseek AFDB (`struct_af`). Produit une shortlist priorisée de cibles druggables structurellement
exploitables. NON une requalification : livrable hand-off drug-discovery (verdicts inchangés).
Run: python analyses/phase52_druggable_structure.py [VI_MAX]   (défaut VI_MAX=-6)
"""
from __future__ import annotations
import json, glob, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "résultats" / "phase52_druggable_structure"
VI_MAX = float(sys.argv[1]) if len(sys.argv) > 1 else -6.0


def structural_support(d):
    """(niveau, détail) — la meilleure évidence structurale disponible."""
    pdb = d.get("pdb") or {}
    if isinstance(pdb, dict) and pdb.get("n_structures"):
        return 3, f"experimental PDB ({pdb['n_structures']})"
    plaf = d.get("plddt_af") or {}
    if isinstance(plaf, dict) and (plaf.get("mean_plddt") or 0) >= 70:
        return 2, f"AlphaFold model pLDDT {round(plaf['mean_plddt'])}"
    sa = d.get("struct_af") or {}
    if isinstance(sa, dict) and sa.get("hits"):
        return 1, "AlphaFold-DB Foldseek hit"
    return 0, "no reliable structure"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        v = (d.get("vulnerability") or {}).get("vi")
        if v is None or v > VI_MAX:
            continue
        lvl, detail = structural_support(d)
        rows.append({"rv": d["rv"], "vi": v, "struct_lvl": lvl, "struct": detail,
                     "verdict": d.get("verdict"), "gene": d.get("gene_name"),
                     "func": (d.get("function_revised") or d.get("product_h37rv") or "")[:60]})
    # priorité : structure expérimentale d'abord, puis VI le plus négatif
    rows.sort(key=lambda r: (-r["struct_lvl"], r["vi"]))
    exploitable = [r for r in rows if r["struct_lvl"] >= 2]  # PDB exp. ou bon modèle AF = docking-ready

    lines = [f"# Cibles druggables structurellement exploitables (CRISPRi VI ≤ {VI_MAX}) — P7.6d", "",
             "Vulnérable au knockdown ET structure résolue = candidat structure-based drug design (docking).",
             "Livrable hand-off drug-discovery (NON une requalification).", "",
             f"Total vulnérables (VI ≤ {VI_MAX}) : {len(rows)} ; dont **docking-ready** (PDB exp. ou modèle AF pLDDT≥70) : "
             f"**{len(exploitable)}**.", "",
             "| VI | gène | verdict | structure | fonction |", "|---|---|---|---|---|"]
    for r in rows:
        g = r["gene"] or r["rv"]
        lines.append(f"| {r['vi']:.2f} | {g} ({r['rv']}) | {r['verdict']} | {r['struct']} | {r['func']} |")
    # gisement NEUF : les docking-ready encore HYPOTHÉTIQUES (verdict dark/family_assigned) = cibles peu caractérisées
    hyp = [r for r in exploitable if r["verdict"] in ("dark", "family_assigned")]
    still_dark = [r for r in hyp if r["verdict"] == "dark"]
    lines += ["", f"## Gisement NEUF : {len(hyp)} cibles docking-ready ENCORE HYPOTHÉTIQUES (peu/pas caractérisées)",
              f"Dont **{len(still_dark)} totalement DARK** (cibles druggables inconnues = candidats deep-dive type Rv3222c) : "
              f"{', '.join(r['rv'] for r in still_dark)}.", "",
              "| VI | gène | verdict | structure | fonction |", "|---|---|---|---|---|"]
    for r in hyp:
        g = r["gene"] or r["rv"]
        lines.append(f"| {r['vi']:.2f} | {g} ({r['rv']}) | {r['verdict']} | {r['struct']} | {r['func']} |")
    (OUT / "druggable_with_structure.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Vulnérables VI≤{VI_MAX} : {len(rows)} | docking-ready : {len(exploitable)} | "
          f"dont hypothétiques : {len(hyp)} (dark : {len(still_dark)} -> {[r['rv'] for r in still_dark]})")
    print(f"Rapport : {OUT/'druggable_with_structure.md'}")


if __name__ == "__main__":
    main()
