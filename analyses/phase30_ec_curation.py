#!/usr/bin/env python3
"""phase30_ec_curation.py -- arbitrage manuel des désaccords EC vs Mycobrowser (P5.11b).

phase28 a isolé 25 « vrais désaccords » EC (accord ni exact ni au niveau classe) entre
l'atlas et Mycobrowser. Ce script porte l'ARBITRAGE MANUEL (curation experte, ancrée sur
UniProt SwissProt = source curée des EC de l'atlas, vérifiée en direct pour les cas ambigus).

Verdict global : les 25 ne sont quasi jamais des erreurs de l'atlas. Catégories :
  - reclassified   : EC renuméroté (même enzyme), 1er chiffre changé → raté par le test classe
                     (ex. pncB 2.4.2.11->6.3.4.21 ; cytochrome oxidase 1.9.3.1->7.1.1.9 ;
                     hélicase 3.6.1.->5.6.2.4).
  - atlas_current  : Mycobrowser obsolète/incomplet/misannoté, l'atlas porte l'EC curé UniProt
                     actuel (ex. gcp « protéase 3.4.24.57 » -> TsaD 2.3.1.234 ; gpm2 -> acid
                     phosphatase ; thiG 4.-.-.- -> 2.8.1.10).
  - subunit_level  : granularité différente, les deux défendables (accD : holoenzyme 6.4.1.-
                     côté Myco vs sous-unité carboxyltransférase 2.1.3.- côté atlas).
  - review_eggnog  : SEUL vrai point à corriger — l'EC atlas vient d'eggNOG et diverge d'UniProt
                     (Rv0112 gca : atlas 1.1.1.281 vs UniProt/Myco GDP-mannose 4,6-déshydratase
                     4.2.1.47) → préférer l'annotation UniProt.

Sortie : résultats/phase30_ec_curation/{ec_curation.json, ec_curation.md}. Fusionne
`mycobrowser.ec_resolution` dans les fiches concernées (le template affiche le verdict au lieu
de « flagged for curation »). Stdlib.
Run: python analyses/phase30_ec_curation.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "phase30_ec_curation"
GENES = ROOT / "site" / "content" / "genes"

# arbitrage expert par gène : (catégorie, note courte)
VERDICTS = {
    "Rv0075":  ("atlas_current", "cysteine-S-conjugate beta-lyase (EC 4.4.1.13, UniProt); Mycobrowser's generic aminotransferase (2.6.1.-) is a fold-based guess"),
    "Rv0112":  ("review_eggnog", "atlas EC 1.1.1.281 is an eggNOG call; UniProt and Mycobrowser both point to GDP-mannose 4,6-dehydratase (4.2.1.47) — prefer the UniProt annotation"),
    "Rv0391":  ("atlas_current", "O-succinylhomoserine sulfhydrylase MetZ (transferase 2.5.1.-); Mycobrowser's 4.2.99.- is superseded"),
    "Rv0417":  ("atlas_current", "thiazole synthase ThiG (EC 2.8.1.10); Mycobrowser's 4.-.-.- is incomplete"),
    "Rv0573c": ("reclassified",  "nicotinate phosphoribosyltransferase: EC 2.4.2.11 was reclassified to 6.3.4.21 (ATP-dependent)"),
    "Rv0958":  ("atlas_current", "magnesium chelatase (EC 6.6.1.1); Mycobrowser's ferrochelatase-class 4.99.1.- is superseded"),
    "Rv0974c": ("subunit_level", "AccD2: Mycobrowser gives the holoenzyme (carboxylase 6.4.1.-), the atlas the carboxyltransferase subunit (2.1.3.-)"),
    "Rv1005c": ("atlas_current", "aminodeoxychorismate synthase (EC 2.6.1.85); Mycobrowser's 4.1.3.- is an older assignment"),
    "Rv1151c": ("atlas_current", "NAD-dependent sirtuin deacylase CobB (EC 2.3.1.286); Mycobrowser's 3.5.1.- predates the modern EC"),
    "Rv1327c": ("atlas_current", "maltosyltransferase GlgE (EC 2.4.99.16, created 2011); Mycobrowser's glycosidase 3.2.1.- is obsolete"),
    "Rv1330c": ("reclassified",  "nicotinate phosphoribosyltransferase: EC 2.4.2.11 was reclassified to 6.3.4.21 (ATP-dependent)"),
    "Rv1602":  ("atlas_current", "imidazole-glycerol-phosphate synthase HisH glutaminase (EC 3.5.1.2 / 4.3.2.10); Mycobrowser's 2.4.2.- is a mis-class"),
    "Rv2199c": ("reclassified",  "cytochrome c oxidase: EC 1.9.3.1 became 7.1.1.9 in the 2018 translocase reclassification"),
    "Rv2604c": ("atlas_current", "glutamine amidotransferase SnoP (EC 3.5.1.2 / 4.3.3.6); Mycobrowser's 2.6.-.- is generic"),
    "Rv3025c": ("atlas_current", "cysteine desulfurase IscS (EC 2.8.1.7); Mycobrowser's 4.4.1.- is the wrong class"),
    "Rv3202c": ("reclassified",  "ATP-dependent DNA helicase: EC 3.6.1.-/3.6.4.12 moved to 5.6.2.4 (motor activity) in 2018"),
    "Rv3214":  ("atlas_current", "acid phosphatase (EC 3.1.3.2/3.1.3.11, UniProt); Mycobrowser's phosphoglycerate mutase 5.4.2.1 is a family-based mis-assignment"),
    "Rv3275c": ("atlas_current", "class II PurE N5-CAIR mutase (EC 5.4.99.18); Mycobrowser's 4.1.1.21 is the old lumped assignment"),
    "Rv3276c": ("atlas_current", "PurK N5-CAIR synthetase (EC 6.3.4.18); Mycobrowser's 4.1.1.21 lumped purK with purE"),
    "Rv3389c": ("atlas_current", "3-hydroxyacyl-thioester dehydratase HtdY (EC 4.2.1.119); Mycobrowser's 1.-.-.- is wrong"),
    "Rv3419c": ("atlas_current", "tRNA threonylcarbamoyltransferase Gcp/TsaD (EC 2.3.1.234); Mycobrowser's 'O-sialoglycoprotein endopeptidase' 3.4.24.57 is a decades-old misannotation"),
    "Rv3538":  ("atlas_current", "enoyl-CoA hydratase ChsH3 (EC 4.2.1.-, UniProt); Mycobrowser's 1.-.-.- (dehydrogenase) is wrong"),
    "Rv3540c": ("atlas_current", "cholesterol-catabolism lyase Ltp2 (EC 4.1.3.-, UniProt); Mycobrowser's thiolase 2.3.1.16 is obsolete"),
    "Rv3552":  ("atlas_current", "cholesterol ring-cleaving hydrolase IpdB (EC 4.1.99.-, UniProt); Mycobrowser's CoA-transferase 2.8.3.- is obsolete"),
    "Rv3799c": ("subunit_level", "AccD4: Mycobrowser gives the holoenzyme (propionyl-CoA carboxylase 6.4.1.3), the atlas the carboxyltransferase subunit (2.1.3.-)"),
}

CAT_LABEL = {
    "reclassified": "EC re-numbering (same enzyme, current class differs)",
    "atlas_current": "atlas uses the current curated EC; Mycobrowser is obsolete/mis-assigned",
    "subunit_level": "subunit vs holoenzyme (both defensible at different granularity)",
    "review_eggnog": "atlas eggNOG EC diverges from UniProt — favour the UniProt annotation",
}


def main():
    print("== phase30 : arbitrage des désaccords EC vs Mycobrowser (P5.11b) ==")
    # récupérer legacy/atlas EC pour le rapport
    info = {}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        if d["rv"] in VERDICTS:
            m = d.get("mycobrowser") or {}
            info[d["rv"]] = {"gene": d.get("gene") or "", "product": d.get("product_h37rv") or "",
                             "legacy_ec": m.get("legacy_ec") or [], "atlas_ec": m.get("atlas_ec") or []}

    from collections import Counter
    cats = Counter(v[0] for v in VERDICTS.values())
    OUT.mkdir(parents=True, exist_ok=True)

    out = {}
    lines = ["# Arbitrage des désaccords EC atlas vs Mycobrowser (P5.11b)", "",
             f"25 désaccords isolés par phase28. Verdict : {dict(cats)}.", "",
             "| Rv | gène | Mycobrowser EC | atlas EC | catégorie | note |",
             "|---|---|---|---|---|---|"]
    for rv in sorted(VERDICTS):
        cat, note = VERDICTS[rv]
        i = info.get(rv, {})
        out[rv] = {"category": cat, "label": CAT_LABEL[cat], "note": note,
                   "legacy_ec": i.get("legacy_ec", []), "atlas_ec": i.get("atlas_ec", [])}
        lines.append(f"| {rv} | {i.get('gene','')} | {','.join(i.get('legacy_ec',[]))} | "
                     f"{','.join(i.get('atlas_ec',[]))} | {cat} | {note} |")
    (OUT / "ec_curation.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
    (OUT / "ec_curation.md").write_text("\n".join(lines) + "\n")
    print(f"Écrit {OUT/'ec_curation.md'} ({len(out)} arbitrés)")
    print("Catégories :", dict(cats))
    print(" -> atlas correct/à jour (reclassified+atlas_current+subunit) :",
          cats["reclassified"] + cats["atlas_current"] + cats["subunit_level"],
          "| à corriger côté atlas (review_eggnog) :", cats["review_eggnog"])

    # fusion ec_resolution dans les fiches
    n = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        if d["rv"] in VERDICTS and d.get("mycobrowser"):
            cat, note = VERDICTS[d["rv"]]
            d["mycobrowser"]["ec_resolution"] = {"category": cat, "label": CAT_LABEL[cat], "note": note}
            Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
            n += 1
    print(f"Fusionné 'ec_resolution' dans {n} fiches")


if __name__ == "__main__":
    main()
