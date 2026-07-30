#!/usr/bin/env python3
"""phase31_ec_override.py -- surcharges EC curées (couche eggNOG divergente d'UniProt) (P5.11c).

Applique les correctifs EC identifiés par l'arbitrage P5.11b où l'EC PRÉDIT (eggNOG) diverge
de la source CURÉE (UniProt). Règle (KB) : quand eggNOG et UniProt divergent sur un EC,
préférer UniProt. Couche de HAND-CURATION traçable (`auto:false` + référence).

Aujourd'hui 1 seule surcharge (Rv0112/gca, seul vrai correctif atlas des 25 désaccords) ; la
structure est prête à en accueillir d'autres (mêmes conflits eggNOG↔UniProt ailleurs).

Sortie : résultats/phase31_ec_override/ec_override.json (keyed Rv). Fusion `ec_override` dans
les fiches + `auto:false` + référence (idempotent) + enregistrement phase4. Stdlib.
Run: python analyses/phase31_ec_override.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "phase31_ec_override"
GENES = ROOT / "site" / "content" / "genes"

OVERRIDES = {
    "Rv0112": {
        "curated_ec": ["4.2.1.47"],
        "superseded_ec": ["1.1.1.281"],
        "superseded_source": "eggNOG",
        "resolution": {
            "category": "eggnog_vs_uniprot",
            "label": "eggNOG EC diverges from UniProt — favour UniProt",
            "note": "eggNOG EC 1.1.1.281 (GDP-4-dehydro-6-deoxy-D-mannose reductase) is superseded by curation to "
                    "EC 4.2.1.47 (GDP-mannose 4,6-dehydratase), aligning with UniProt (Gca) and Mycobrowser",
        },
        "references": [
            {"authors": "UniProt Consortium", "year": 2024,
             "title": "UniProtKB O53634 (Gca, possible GDP-mannose 4,6-dehydratase, Mycobacterium tuberculosis)",
             "journal": "UniProt", "doi": "10.1093/nar/gkac1052"},
        ],
    },
}


def main():
    print("== phase31 : surcharges EC curées eggNOG->UniProt (P5.11c) ==")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "ec_override.json").write_text(json.dumps(OVERRIDES, ensure_ascii=False, indent=1))

    n = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        ov = OVERRIDES.get(d["rv"])
        if not ov:
            continue
        d["ec_override"] = ov
        d["auto"] = False            # désormais hand-reviewed sur l'EC
        d["needs_review"] = False
        refs = d.get("references") or []
        have = {r.get("title") for r in refs}
        for r in ov["references"]:
            if r["title"] not in have:
                refs.append(r)
        d["references"] = refs
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
        print(f"  {d['rv']} ({d.get('gene')}) : EC {ov['superseded_ec']} ({ov['superseded_source']}) "
              f"-> {ov['curated_ec']} (curated) ; auto=false")
    print(f"Fusionné 'ec_override' dans {n} fiche(s)")


if __name__ == "__main__":
    main()
