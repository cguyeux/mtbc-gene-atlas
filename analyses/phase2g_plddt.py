#!/usr/bin/env python3
"""phase2g_plddt.py -- Model-confidence (pLDDT) gating for the Foldseek layer.

The Foldseek structural neighbours (phase2c) were searched on per-gene ESMFold
models. ESMFold writes its per-residue confidence (pLDDT, 0-100) into the PDB
B-factor column. A fold hypothesis drawn from a low-confidence model is weak:
a Foldseek hit against a mostly-disordered or poorly-predicted model can be an
artefact. This phase computes the mean pLDDT of each model so the site can gate
the structural claim (and a curator can discount low-confidence hits).

Bands (AlphaFold/ESMFold convention): >=90 very high, 70-90 confident,
50-70 low, <50 very low.

Inputs:  résultats/phase2c_foldseek/<rv>.pdb
Output:  résultats/phase2c_foldseek/plddt.json  (keyed by rv; read by phase4)
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PDB_DIR = ROOT / "résultats" / "phase2c_foldseek"


def mean_plddt(pdb: Path) -> tuple[float, int]:
    vals = []
    for line in pdb.read_text().split("\n"):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                continue
    if not vals:
        return 0.0, 0
    # ESMFold may write pLDDT on a 0-1 scale (B-factor ~0.5) or 0-100. Normalise to 0-100.
    if max(vals) <= 1.5:
        vals = [v * 100 for v in vals]
    return sum(vals) / len(vals), len(vals)


def band(p: float) -> str:
    if p >= 90:
        return "very high"
    if p >= 70:
        return "confident"
    if p >= 50:
        return "low"
    return "very low"


def main() -> None:
    out = {}
    for pdb in sorted(PDB_DIR.glob("*.pdb")):
        rv = pdb.stem
        mp, n = mean_plddt(pdb)
        out[rv] = {"mean_plddt": round(mp, 1), "n_residues": n, "band": band(mp)}
    (PDB_DIR / "plddt.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    lo = sum(1 for v in out.values() if v["mean_plddt"] < 70)
    print(f"[plddt] {len(out)} models | {lo} below 70 (low-confidence) "
          f"| mean of means {sum(v['mean_plddt'] for v in out.values())/max(1,len(out)):.1f}")
    print(f"Wrote {PDB_DIR / 'plddt.json'}")


if __name__ == "__main__":
    main()
