#!/usr/bin/env python3
"""phase50_ptm.py -- P7.8 : couche modifications post-traductionnelles (phosphoprotéome + PTM).

Dimension jusqu'ici ABSENTE. Un site PTM (surtout phospho, médié par les 11 STPK de M.tb) sur un
hypothétique = preuve d'existence protéique (orthogonale à la MS) + signe de régulation (substrat de
signalisation). Ce n'est PAS une fonction (garde-fou : « régulé par phosphorylation » ≠ fonction moléculaire).

Source : **UniProt** (feature « Modified residue »), API stream par organisme (org 83332), mapping
accession→Rv via la couche uniprot déjà en base. dbPTM (source complète Prisic 2010) est bloqué (403) →
la couche complète du phosphoprotéome reste à faire (P7.8b, parser suppléments Prisic/Fortuin). Rendement
UniProt : 265 protéines / 413 sites (170 phospho), dont 71 hypothétiques.
Couche `ptm` sur `Gene`. Run: python analyses/phase50_ptm.py
"""
from __future__ import annotations
import json, glob, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
REF = {"source": "UniProt (Modified residue features; PTM sites curated from the M. tuberculosis literature)",
       "db": "UniProtKB", "url": "https://www.uniprot.org"}


def main():
    acc2rv = {}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        a = (d.get("uniprot") or {}).get("acc")
        if a:
            acc2rv[a] = d["rv"]

    url = ("https://rest.uniprot.org/uniprotkb/stream?query=organism_id:83332+AND+ft_mod_res:*"
           "&fields=accession,ft_mod_res&format=json")
    data = json.loads(urllib.request.urlopen(url, timeout=120).read().decode())

    n = 0
    for r in data.get("results", []):
        rv = acc2rv.get(r.get("primaryAccession"))
        if not rv:
            continue
        sites = []
        for ft in (r.get("features") or []):
            if ft.get("type") == "Modified residue":
                sites.append({"pos": ft.get("location", {}).get("start", {}).get("value"),
                              "type": ft.get("description", "")})
        if not sites:
            continue
        nphos = sum(1 for s in sites if "phospho" in s["type"].lower())
        f = GENES / f"{rv}.json"
        d = json.load(open(f))
        d["ptm"] = {"n_sites": len(sites), "n_phospho": nphos, "sites": sites,
                    "note": ("Experimentally reported post-translational modification(s). A phosphosite indicates the "
                             "protein is expressed and is a substrate of the M. tuberculosis Ser/Thr/Tyr kinase "
                             "signalling network — a regulatory context, NOT a molecular function."),
                    "source": REF["source"], "reference": REF}
        f.write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
    print(f"P7.8 : couche `ptm` écrite sur {n} fiches (UniProt Modified residue).")


if __name__ == "__main__":
    main()
