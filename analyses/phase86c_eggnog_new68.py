#!/usr/bin/env python3
"""phase86c-eggnog — parse eggNOG-mapper output into the `eggnog` layer for the 68
Mycobrowser genes (P16.5d-cont). Reads résultats/phase86c_eggnog/new68.emapper.annotations
(produced by emapper diamond mode) and writes the same-shape `eggnog` layer as phase2d:
OG / COG category / description / preferred name / EC / KEGG KO / pathway / module / CAZy / GO.
Upgrades a still-dark fiche to family_assigned when eggNOG yields a real functional handle
(a named OG/description or EC/KO), setting function_revised from the eggNOG description if empty.
Idempotent.
"""
from __future__ import annotations
import json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANN = os.path.join(ROOT, "résultats", "phase86c_eggnog", "new68.emapper.annotations")
CG = os.path.join(ROOT, "site", "content", "genes")
COLS = ["query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
        "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
        "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
        "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"]
REF = {"authors": "Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P, Huerta-Cepas J",
       "year": 2021, "title": "eggNOG-mapper v2", "journal": "Mol Biol Evol 38:5825-5829"}


def lst(v):
    v = (v or "").strip()
    return [] if v in ("", "-") else [x.strip() for x in v.split(",") if x.strip() and x.strip() != "-"]


def clean(v):
    v = (v or "").strip()
    return None if v in ("", "-") else v


def pick_og(ogs):
    """eggNOG_OGs like 'COG2897@1|root,2A3F1@2|Bacteria,...' -> a COG id if present, else first token."""
    toks = [t.split("@")[0] for t in lst(ogs)]
    for t in toks:
        if t.startswith("COG"):
            return t
    return toks[0] if toks else None


def main():
    if not os.path.exists(ANN):
        print(f"[phase86c-eggnog] annotations absentes: {ANN} (emapper non abouti ?)")
        return
    rows = {}
    for line in open(ANN):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < len(COLS):
            f += [""] * (len(COLS) - len(f))
        d = dict(zip(COLS, f))
        rows[d["query"]] = d

    written = upgraded = 0
    for rv, d in rows.items():
        fn = os.path.join(CG, f"{rv}.json")
        if not os.path.exists(fn):
            continue
        rec = json.load(open(fn))
        if not str(rec.get("curation_note", "")).startswith("Added P16.5d"):
            continue
        layer = {
            "og": pick_og(d["eggNOG_OGs"]),
            "cog_cat": clean(d["COG_category"]),
            "description": clean(d["Description"]),
            "preferred_name": clean(d["Preferred_name"]),
            "ec": lst(d["EC"]),
            "kegg_ko": lst(d["KEGG_ko"]),
            "kegg_pathway": lst(d["KEGG_Pathway"]),
            "kegg_module": lst(d["KEGG_Module"]),
            "cazy": lst(d["CAZy"]),
            "go": lst(d["GOs"]),
            "source": "eggNOG-mapper v2 (diamond), eggNOG 5.0",
            "refs": [REF],
        }
        rec["eggnog"] = layer
        # real functional handle -> upgrade a still-dark fiche (curated/orthology evidence)
        handle = layer["ec"] or layer["kegg_ko"] or (layer["cog_cat"] and layer["cog_cat"] not in ("S",)) \
            or (layer["description"] and "unknown" not in layer["description"].lower())
        if rec.get("verdict") == "dark" and handle:
            rec["verdict"] = "family_assigned"
            rec["confidence"] = rec.get("confidence") or "low"
            if not rec.get("function_revised") and layer["description"]:
                rec["function_revised"] = layer["description"]
            upgraded += 1
        json.dump(rec, open(fn, "w"), ensure_ascii=False, indent=2)
        written += 1

    n_og = sum(1 for d in rows.values() if pick_og(d["eggNOG_OGs"]))
    n_ec = sum(1 for d in rows.values() if lst(d["EC"]))
    n_ko = sum(1 for d in rows.values() if lst(d["KEGG_ko"]))
    print(f"phase86c-eggnog: couche eggnog écrite sur {written} gènes (sur 68)")
    print(f"  avec OG: {n_og} | avec EC: {n_ec} | avec KEGG KO: {n_ko}")
    print(f"  dark -> family_assigned (handle eggNOG réel): {upgraded}")


if __name__ == "__main__":
    main()
