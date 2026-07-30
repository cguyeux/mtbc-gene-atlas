#!/usr/bin/env python3
"""phase28_mycobrowser_concordance.py -- vérification champ-par-champ vs Mycobrowser (P5.11).

Cœur de la mission (« partir des dires actuels de Mycobrowser, les VÉRIFIER »). Compare, par
gène, l'enregistrement legacy Mycobrowser (Product, Function, EC, Functional_Category) aux
annotations de l'atlas, et produit :
  (a) le DELTA de requalification : parmi les « conserved hypotheticals » de Mycobrowser, ceux
      auxquels l'atlas donne un handle fonctionnel (= atlas DEVANT l'original) ;
  (b) la concordance EC, en distinguant honnêtement l'accord, la RECLASSIFICATION de
      nomenclature (les EC de Mycobrowser sont périmés : il n'est plus maintenu, et la classe
      EC 7 « translocases » (2018) n'existait pas → 3.6.3.x -> 7.x = même enzyme, numéro à jour),
      et le rare VRAI désaccord (cible de curation).

GARDE-FOU (leçon KB scientific-writing) : ne PAS présenter l'accord sur gènes déjà caractérisés
comme une « validation » des requalifications. Le delta hypothétique est la contribution ; l'accord
EC est un simple contrôle de cohérence, avec le caveat reclassification explicite.

Source : txt de release Mycobrowser rapatrié data/mycobrowser/H37Rv.txt (P5.4).
Sortie : résultats/phase28_concordance/{concordance.json (keyed Rv), summary.json}. Fusion du
champ `mycobrowser` dans les fiches (idempotent) + enregistrement phase4. Stdlib.
Run: python analyses/phase28_mycobrowser_concordance.py
"""
from __future__ import annotations
import json, csv, glob
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent
MYCO = ROOT / "data" / "mycobrowser" / "H37Rv.txt"
OUT = ROOT / "résultats" / "phase28_concordance"
GENES = ROOT / "site" / "content" / "genes"

SOURCE = "Field-by-field comparison against the Mycobrowser legacy record (release 5)"
REFS = [
    {"authors": "Kapopoulou A, Lew JM, Cole ST", "year": 2011,
     "title": "The MycoBrowser portal: a comprehensive and manually annotated resource for mycobacterial genomes",
     "journal": "Tuberculosis", "doi": "10.1016/j.tube.2010.09.006"},
]


def atlas_handle(d):
    """Handles fonctionnels forts de l'atlas (mêmes signaux que 'ahead of UniProt', discipliné)."""
    up = d.get("uniprot") or {}
    eg = d.get("eggnog") or {}
    t = []
    if up.get("function"):
        t.append("curated function (UniProt)")
    if up.get("ec") or eg.get("ec"):
        t.append("EC number")
    cog = eg.get("cog_cat") or eg.get("cog") or ""
    if cog and cog not in ("S", "", None):
        t.append("COG category")
    if d.get("function_revised") and d.get("verdict") == "requalified":
        t.append("requalified function")
    return t


def ec_status(myco_ec: set, atlas_ec: set):
    if not myco_ec or not atlas_ec:
        return None
    if myco_ec & atlas_ec:
        return "agree"
    mc = {x.split(".")[0] for x in myco_ec}
    ac = {x.split(".")[0] for x in atlas_ec}
    # reclassification connue : ATPases translocases 3.6.3.x -> classe 7 (2018)
    if mc & ac or ("3" in mc and "7" in ac) or ("7" in mc and "3" in ac):
        return "reclassified"
    return "differ"


def main():
    print("== phase28 : concordance champ-par-champ vs Mycobrowser (P5.11) ==")
    myco = {r["Locus"]: r for r in csv.DictReader(open(MYCO, encoding="utf-8"), delimiter="\t")}

    rec = {}
    agg = Counter()
    ec_agg = Counter()
    ahead_examples = []
    ec_differ = []
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        rv = d["rv"]
        r = myco.get(rv)
        if not r:
            continue
        fc = (r.get("Functional_Category") or "").strip()
        is_ch = fc.lower() == "conserved hypotheticals"
        handle = atlas_handle(d)
        ahead = bool(is_ch and handle)

        myco_ec = set(x.strip() for x in (r.get("Enzyme Classification") or "").split(",") if x.strip())
        atlas_ec = set((d.get("uniprot") or {}).get("ec") or []) | set((d.get("eggnog") or {}).get("ec") or [])
        ecs = ec_status(myco_ec, atlas_ec)

        func = (r.get("Function") or "").strip()
        rec[rv] = {
            "legacy_product": (r.get("Product") or "").strip(),
            "legacy_function": func[:400] if func and func.lower() != "unknown" else "",
            "legacy_ec": sorted(myco_ec),
            "conserved_hypothetical": is_ch,
            "atlas_ahead": ahead,
            "ahead_via": handle if ahead else [],
            "ec_status": ecs,
            "atlas_ec": sorted(atlas_ec),
            "source": SOURCE, "refs": REFS,
        }

        agg["compared"] += 1
        if is_ch:
            agg["conserved_hypothetical"] += 1
            if ahead:
                agg["ahead"] += 1
                if len(ahead_examples) < 8:
                    ahead_examples.append((rv, handle[0]))
        if ecs:
            ec_agg[ecs] += 1
            if ecs == "differ" and len(ec_differ) < 12:
                ec_differ.append((rv, sorted(myco_ec), sorted(atlas_ec)))

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "concordance.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))

    summary = {
        "compared": agg["compared"],
        "conserved_hypothetical": agg["conserved_hypothetical"],
        "atlas_ahead_on_hypothetical": agg["ahead"],
        "atlas_ahead_pct": round(100 * agg["ahead"] / max(agg["conserved_hypothetical"], 1), 1),
        "ahead_examples": ahead_examples,
        "ec_compared": sum(ec_agg.values()),
        "ec_agree": ec_agg["agree"],
        "ec_reclassified": ec_agg["reclassified"],
        "ec_differ": ec_agg["differ"],
        "ec_differ_examples": ec_differ,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1))

    print(f"Comparés : {summary['compared']}")
    print(f"Conserved hypotheticals Mycobrowser : {summary['conserved_hypothetical']} ; "
          f"atlas DEVANT (handle fonctionnel) : {summary['atlas_ahead_on_hypothetical']} "
          f"({summary['atlas_ahead_pct']} %)")
    print(f"EC comparés : {summary['ec_compared']} | accord {summary['ec_agree']} | "
          f"reclassification {summary['ec_reclassified']} | VRAI désaccord {summary['ec_differ']}")
    print("Désaccords EC (curation) :", ec_differ[:6])

    # fusion dans les fiches
    n = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        d["mycobrowser"] = rec.get(d["rv"], {})
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
    print(f"Fusionné 'mycobrowser' dans {n} fiches")


if __name__ == "__main__":
    main()
