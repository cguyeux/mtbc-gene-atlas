#!/usr/bin/env python3
"""phase6_review_structural.py -- Human review of the 26 structure-only promotions.

Each of the 26 genes that phase5 promoted dark->family on a Foldseek hit alone
(needs_review=True) was inspected against its best structural hit. The decisions
below move them into the MANUAL curation (data/pilot_curation.json), which is
stable and takes priority, so they lose the needs_review flag.

Decision key per gene: verdict + a reviewed function note. The evidence is
rebuilt from PGAP + Pfam + the Foldseek best hit, tagged as reviewed.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
XREF = ROOT / "data" / "gene_xref.tsv"
PFAM = ROOT / "résultats" / "phase2b_pfam" / "pfam.json"
FOLD = ROOT / "résultats" / "phase2c_foldseek" / "foldseek.json"
MANUAL = ROOT / "data" / "pilot_curation.json"

# verdict, confidence, reviewed function note
DECISIONS = {
    "Rv0398c": ("family_assigned", "low", "Secreted protein with a solved structure (PDB 7YD4, Rv0398c itself); fold characterised, precise function not established."),
    "Rv0455c": ("family_assigned", "low", "DUF5078; structure solved (PDB 8DRI, Rv0455c); fold characterised, function not established."),
    "Rv0481c": ("family_assigned", "low", "DUF2505; ESMFold matches a START-like / Ups1 lipid-transfer fold (PDB 5JQM); putative lipid-binding protein."),
    "Rv0540":  ("family_assigned", "low", "DUF2064; fold of a putative nucleotide-diphospho-sugar transferase (PDB 3CGX); putative sugar/nucleotidyl transferase."),
    "Rv0603":  ("family_assigned", "low", "Secreted protein with a solved NMR structure (PDB 2LRA, Rv0603); fold characterised, function not established."),
    "Rv0633c": ("dark", "low", "Conserved hypothetical; Foldseek best hit (KirBac channel, TM 0.43, 4.1A) not conclusive. Function unknown."),
    "Rv0965c": ("family_assigned", "low", "ESX/WXG100-like secreted fold (matches EsxB, PDB 4J7K, TM 0.82); putative ESX-secreted protein."),
    "Rv0966c": ("dark", "low", "DUF1707; Foldseek hits non-conclusive (moderate TM). Function unknown."),
    "Rv1006":  ("family_assigned", "low", "Glycoside-hydrolase-like fold (GH2, PDB 5T99); putative glycosidase."),
    "Rv1097c": ("family_assigned", "low", "Apa-like fold (fibronectin-binding protein Apa, PDB 5ZX9); putative adhesin / secreted protein."),
    "Rv1269c": ("family_assigned", "low", "DUF4189; same fold as Rv1813c (PDB 7NHZ); secreted DUF4189-family protein."),
    "Rv1312":  ("dark", "low", "DUF2550; Foldseek best hit (PH domain, TM 0.46) not conclusive. Function unknown."),
    "Rv1480":  ("family_assigned", "low", "DUF58; fold of the BatB/DUF58 family (PDB 3IBS); putative integrity/assembly protein."),
    "Rv1502":  ("family_assigned", "low", "Glycoside-hydrolase fold (PDB 5MUI); putative glycosidase."),
    "Rv1510":  ("family_assigned", "low", "MATE multidrug-transporter fold (PDB 6FHZ); putative membrane transporter."),
    "Rv1632c": ("family_assigned", "low", "DUF402; structure of the DUF402 family (PDB 2P12, TM 0.96); putative nucleotide-binding/phosphatase fold."),
    "Rv1691":  ("family_assigned", "low", "TPR-repeat fold (PDB 2PL2); putative protein-protein-interaction scaffold."),
    "Rv1724c": ("dark", "low", "Conserved hypothetical; Foldseek hits non-conclusive (moderate TM). Function unknown."),
    "Rv1754c": ("requalified", "medium", "Endo-alpha-D-arabinanase / glycoside hydrolase (DUF4185). Strong structural match to arabinanase EndoMA1 (PDB 8IC1, E=4e-36, TM 0.84); putative arabinan-degrading glycosidase."),
    "Rv1778c": ("dark", "low", "Conserved hypothetical; Foldseek best hit (RbmA, TM 0.38) not conclusive. Function unknown."),
    "Rv1780":  ("family_assigned", "low", "BPI/SPLUNC1-like lipid-binding fold (PDB 5I7L); putative lipid-binding protein."),
    "Rv1813c": ("family_assigned", "low", "DUF4189; structure solved (PDB 7NHZ, Rv1813c); secreted DUF4189-family protein."),
    "Rv1815":  ("family_assigned", "low", "Protease-like fold (L5 protease, PDB 5MRT); putative peptidase."),
    "Rv1887":  ("family_assigned", "low", "Exopolyphosphatase-like fold (PDB 3HI0); putative phosphatase/hydrolase."),
    "Rv1948c": ("family_assigned", "low", "Suppressor-of-Fused-like fold (PDB 4KMA); fold characterised, function not established."),
    "Rv2033c": ("dark", "low", "DUF3097; Foldseek best hit (OLD protein, TM 0.49) not conclusive. Function unknown."),
}


def main() -> None:
    xref = {r["rv"]: r for r in csv.DictReader(open(XREF), delimiter="\t")}
    pfam = json.loads(PFAM.read_text())
    fold = json.loads(FOLD.read_text())
    manual = json.loads(MANUAL.read_text())

    n_fam = n_dark = n_req = 0
    for rv, (verdict, conf, note) in DECISIONS.items():
        row = xref[rv]
        ev = [f"MTBC0 PGAP product: {row['product_mtbc0_pgap']}"]
        if pfam.get(rv):
            ev.append("Pfam: " + ", ".join(f"{d['pfam_name']} {d['pfam_acc']}" for d in pfam[rv]))
        if fold.get(rv):
            h = fold[rv][0]
            ev.append(f"Foldseek best: {h['description'][:80]} (prob {h.get('prob',0):.2f}, "
                      f"E={h.get('evalue',0):.0e}, TM={h.get('tmscore',0):.2f})")
        ev.append("(structure-only promotion reviewed by hand, 2026-06-01)")
        manual[rv] = {"verdict": verdict, "confidence": conf, "function_revised": note,
                      "evidence": ev, "references": []}
        n_fam += verdict == "family_assigned"
        n_dark += verdict == "dark"
        n_req += verdict == "requalified"

    MANUAL.write_text(json.dumps(manual, indent=2, ensure_ascii=False))
    print(f"Reviewed {len(DECISIONS)} structural promotions -> pilot_curation.json")
    print(f"  kept family={n_fam}  requalified={n_req}  demoted to dark={n_dark}")


if __name__ == "__main__":
    main()
