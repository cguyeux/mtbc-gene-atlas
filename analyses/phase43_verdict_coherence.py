#!/usr/bin/env python3
"""phase43_verdict_coherence.py -- P7.10c : aligner verdict↔function_revised.

Audit (né du cas Rv3527) : des fiches `verdict=dark` portent une function_revised DÉJÀ curée
(auto:false), issue des modules compagnon (conserved_orphan_modules : cholestérol/phoP/arabinane/
Mce/cobalamine) ou d'antigènes caractérisés — mais le verdict n'a jamais suivi (bug de propagation,
cf. P6). On corrige le VERDICT uniquement (function_revised inchangée, elle est déjà bonne).

Rv0966c EXCLU : sa fonction dit « Function unknown » → dark légitime.
Run: python analyses/phase43_verdict_coherence.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"

# rv -> confidence : 'medium' pour rôle/antigène établi, 'low' pour « candidate accessory »
FIX = {
    "Rv1974": "low", "Rv0531": "low", "Rv3572": "low", "Rv0756c": "low",
    "Rv2206": "low", "Rv2342": "low",          # candidate accessory of a named module
    "Rv2645": "medium", "Rv2348c": "medium",   # characterised T-cell antigen
}
NOTE = ("2026-07-04 (P7.10c): verdict aligned to the pre-existing hand-curated function_revised "
        "(a functional role was already documented but the verdict had remained 'dark').")


def main():
    print("== phase43 : cohérence verdict↔fonction (P7.10c) ==")
    n = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        conf = FIX.get(d["rv"])
        if not conf:
            continue
        assert d.get("verdict") == "dark", f"{d['rv']} n'est plus dark ({d.get('verdict')})"
        d["verdict"] = "family_assigned"
        d["confidence"] = conf
        d["needs_review"] = False
        d["curation_note"] = NOTE
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
        print(f"  {d['rv']} -> family_assigned/{conf} :: {(d.get('function_revised') or '')[:60]}")
    print(f"{n} verdict(s) corrigé(s).")


if __name__ == "__main__":
    main()
