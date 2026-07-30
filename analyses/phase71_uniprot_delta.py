#!/usr/bin/env python3
"""phase71_uniprot_delta.py -- P16.2b : delta UniProt (volet « péremption » de l'audit de fraîcheur).

UniProt bouge ; notre couche `uniprot` (phase2e) est un SNAPSHOT. Ce script re-télécharge l'état COURANT du
protéome H37Rv (UP000001584) en UN appel (endpoint `stream`, pas 3883 requêtes), le compare fiche à fiche, et
rafraîchit chirurgicalement les 4 champs qui peuvent avoir bougé : `reviewed`, `protein_name`, `ec`, `function`.
Les autres champs du snapshot (go, kegg, evidence, acc) sont PRÉSERVÉS.

    curl "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000001584\
&fields=accession,reviewed,protein_name,ec,cc_function&format=tsv"

RÉSULTAT du premier passage (2026-07-13) :
  - **0 des 219 gènes `dark` n'a gagné la moindre curation UniProt** depuis le snapshot → négatif propre : même
    les curateurs UniProt n'ont pas bougé dessus, le résidu dark est génuinement bloqué.
  - Delta réel sur les gènes DÉJÀ annotés : 5 promus SwissProt, 11 fonctions curées nouvelles, 6 EC nouveaux.
    L'atlas n'était PAS en retard sur la famille (il l'avait dérivée seul) mais l'était sur l'ACTIVITÉ précise :
      Rv2969c +EC 1.8.4.-   (disulfide oxidase DsbA, plus précis que « mycothiol-dependent reductase »)
      Rv2968c +EC 1.17.4.-  (VKOR)
      Rv2864c +EC 3.4.16.4  (D-Ala-D-Ala carboxypeptidase, plus précis que « penicillin-binding protein »)
      Rv1912c +EC 1.6.5.5   (NADPH:quinone reductase, plus précis que « famille MDR »)
      Rv2092c  EC PÉRIMÉ    3.6.4.- -> 5.6.2.6 (reclassification des hélicases)

GARDE-FOU : on rafraîchit une couche de DONNÉES (le snapshot d'une base externe), on ne réécrit AUCUN verdict ni
aucune fonction curée à la main. Un gain d'EC est signalé (`ec_gained`) pour relecture, pas propagé en verdict.

Écrit EN PLACE. Run: python analyses/phase71_uniprot_delta.py  (le TSV doit être téléchargé au préalable)
"""
from __future__ import annotations
import csv, glob, json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
TSV = ROOT / "résultats" / "phase71_uniprot_delta" / "uniprot_current.tsv"
TODAY = "2026-07-13"

ECO = re.compile(r"\s*\{ECO:[^}]*\}")
PREFIX = re.compile(r"^FUNCTION:\s*")


def clean_function(s: str) -> str:
    """« FUNCTION: ... {ECO:0000255|HAMAP...} » -> texte nu."""
    s = PREFIX.sub("", (s or "").strip())
    s = ECO.sub("", s)
    return s.strip()


def split_ec(s: str) -> list[str]:
    return [e.strip() for e in (s or "").split(";") if e.strip()]


def main() -> None:
    if not TSV.exists():
        raise SystemExit(f"{TSV} absent — télécharger d'abord l'état UniProt courant (voir docstring).")
    cur = {}
    with TSV.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cur[r["Entry"]] = r

    changed = 0
    report = {"newly_reviewed": [], "new_function": [], "ec_gained": [], "ec_changed": [], "dark_moved": []}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.loads(Path(f).read_text())
        up = d.get("uniprot") or {}
        acc = up.get("acc")
        if not acc or acc not in cur:
            continue
        c = cur[acc]
        rv, verdict = d["rv"], d.get("verdict")

        new_rev = c["Reviewed"] == "reviewed"
        new_name = (c["Protein names"] or "").strip()
        new_ec = split_ec(c["EC number"])
        new_func = clean_function(c["Function [CC]"])

        old_rev = bool(up.get("reviewed"))
        old_ec = up.get("ec") or []
        old_func = (up.get("function") or "").strip()

        delta = False
        if new_rev != old_rev:
            report["newly_reviewed"].append((rv, verdict, new_name[:50])) if new_rev else None
            delta = True
        if new_func and not old_func:
            report["new_function"].append((rv, verdict, new_func[:60]))
            delta = True
        if new_ec and not old_ec:
            report["ec_gained"].append((rv, verdict, ";".join(new_ec)))
            delta = True
        elif new_ec and old_ec and set(new_ec) != set(old_ec):
            report["ec_changed"].append((rv, verdict, f"{old_ec} -> {new_ec}"))
            delta = True
        if delta and verdict == "dark":
            report["dark_moved"].append((rv, new_name[:50]))

        if delta:
            up["reviewed"] = new_rev
            up["protein_name"] = new_name
            up["ec"] = new_ec
            up["function"] = new_func
            up["refreshed"] = TODAY
            d["uniprot"] = up
            Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
            changed += 1

    print(f"couche `uniprot` rafraîchie sur {changed} fiches (snapshot -> UniProt du {TODAY}).\n")
    for k in ("newly_reviewed", "new_function", "ec_gained", "ec_changed"):
        print(f"{k}: {len(report[k])}")
        for x in report[k][:8]:
            print(f"    {x[0]:9} ({x[1]}) {x[2]}")
    print(f"\n>>> gènes DARK ayant bougé chez UniProt : {len(report['dark_moved'])}")
    if not report["dark_moved"]:
        print("    AUCUN. Négatif propre : même les curateurs UniProt n'ont rien pu dire des 219 dark.")
    print("\nGarde-fou : couche de DONNÉES rafraîchie ; aucun verdict ni fonction curée n'est réécrit.")


if __name__ == "__main__":
    main()
