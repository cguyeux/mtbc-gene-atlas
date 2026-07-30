#!/usr/bin/env python3
"""phase86f — STRING v12 functional-interaction layer for the 68 Mycobrowser genes (P16.5d-cont).

Correction: STRING was wrongly marked "blocked" for these genes. STRING's M. tuberculosis
proteome (taxon 83332) has 4054 protein IDs — it INCLUDES the A-suffixed / Mycobrowser loci,
and all 68 P16.5d genes are present. So STRING is feasible, reusing the already-downloaded bulk
(résultats/phase2h_string/83332.protein.links.detailed.v12.0.txt.gz) and phase2h's validated
noisy-OR / anchor logic — no re-download.

Reuses phase2h.{load_xref, parse_info, parse_links, build_record} and merges the `string` layer
into the 68 fiches only. STRING is guilt-by-association: the anchor is a functional LEAD, not a
function; it never changes a verdict here (consistent with the whole-proteome STRING pass).
Idempotent.
"""
from __future__ import annotations
import importlib.util
import json, glob, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CG = os.path.join(ROOT, "site", "content", "genes")


def main():
    spec = importlib.util.spec_from_file_location(
        "p2h", os.path.join(ROOT, "analyses", "phase2h_string.py"))
    p2h = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(p2h)

    if not p2h.LINKS_GZ.exists():
        raise SystemExit(f"bulk STRING absent: {p2h.LINKS_GZ} — lancer phase2h d'abord")
    xref = p2h.load_xref()
    info = p2h.parse_info()
    edges_by_rv, n_edges, prior = p2h.parse_links()
    print(f"phase86f: STRING bulk parsé ({n_edges} arêtes, prior {prior:.3f})")

    targets = [json.load(open(f))["rv"] for f in glob.glob(os.path.join(CG, "Rv*.json"))
               if str(json.load(open(f)).get("curation_note", "")).startswith("Added P16.5d")]
    written = n_anchor = 0
    for rv in targets:
        edges = edges_by_rv.get(rv, [])
        if not edges:
            continue
        rec = p2h.build_record(rv, edges, xref, info)
        fn = os.path.join(CG, f"{rv}.json")
        d = json.load(open(fn))
        d["string"] = rec                       # STRING = guilt-by-association; verdict UNCHANGED
        json.dump(d, open(fn, "w"), ensure_ascii=False, indent=2)
        written += 1
        if rec.get("anchor"):
            n_anchor += 1

    print(f"phase86f: couche string écrite sur {written}/{len(targets)} gènes")
    print(f"  avec anchor guilt-by-association (lead, pas flip): {n_anchor}")
    # show a few anchors for the still-dark ones
    for rv in targets:
        d = json.load(open(os.path.join(CG, f"{rv}.json")))
        a = (d.get("string") or {}).get("anchor")
        if d.get("verdict") == "dark" and a:
            print(f"    [dark] {rv} -> anchor {a['rv']} ({a.get('gene') or a['product'][:30]}), "
                  f"no_tm={a['combined_no_tm']}")


if __name__ == "__main__":
    main()
