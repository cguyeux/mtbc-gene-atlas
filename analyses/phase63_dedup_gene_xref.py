#!/usr/bin/env python3
"""P14.2.1 — Dédoublonnage idempotent de site/content/gene_xref.tsv.

Contexte : le catalogue brut gene_xref.tsv contient 31 locus en DOUBLON (une ligne
fragment/copie + la vraie), cause racine du bug de coordonnées corrigé en v53 (Rv2082).
Le site servi est déjà correct (la fiche JSON enrichie prime à l'ingest), mais le fichier
catalogue reste sale. Ce script nettoie la SOURCE.

Règle de choix (déterministe) : pour chaque rv, garder la ligne qui correspond à la fiche
JSON enrichie (source de vérité), en priorisant l'accord de COORDONNÉES (start+end) sur
l'accord de SÉQUENCE protéique (utile quand deux copies génomiques partagent la protéine).
En cas d'égalité, garder la première rencontrée. Idempotent : relancé sur un fichier déjà
dédoublonné, il ne change rien.

Usage : python3 phase63_dedup_gene_xref.py [--write]
Sans --write : dry-run (n'écrit rien, affiche le plan). Avec --write : sauvegarde l'original
en .predup puis réécrit le fichier dédoublonné.
"""
import csv
import json
import sys
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
XREF = ROOT / "site" / "content" / "gene_xref.tsv"
GENES = ROOT / "site" / "content" / "genes"


def json_truth(rv):
    p = GENES / f"{rv}.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    return d.get("start_mtbc0"), d.get("end_mtbc0"), d.get("protein_mtbc0")


def score(row, truth):
    """2 = accord coordonnées, 1 = accord protéine, 0 = ni l'un ni l'autre."""
    if not truth:
        return 0
    js, je, jp = truth
    rs = int(row["start_mtbc0"]) if row.get("start_mtbc0") else None
    re_ = int(row["end_mtbc0"]) if row.get("end_mtbc0") else None
    if js is not None and rs == js and re_ == je:
        return 2
    if jp and row.get("protein_mtbc0") and row["protein_mtbc0"] == jp:
        return 1
    return 0


def main():
    write = "--write" in sys.argv
    with XREF.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames
        rows = list(reader)

    # Regrouper par rv en conservant l'ordre d'apparition.
    by_rv = OrderedDict()
    for row in rows:
        by_rv.setdefault(row["rv"], []).append(row)

    kept = []
    dropped = []
    n_dup = 0
    for rv, group in by_rv.items():
        if len(group) == 1:
            kept.append(group[0])
            continue
        n_dup += 1
        truth = json_truth(rv)
        best_i = max(range(len(group)), key=lambda i: score(group[i], truth))
        kept.append(group[best_i])
        for i, row in enumerate(group):
            if i != best_i:
                dropped.append((rv, row.get("start_mtbc0"), row.get("len_aa"),
                                score(row, truth)))

    print(f"lignes avant     : {len(rows)}")
    print(f"rv uniques       : {len(by_rv)}")
    print(f"rv dédoublonnés  : {n_dup}")
    print(f"lignes supprimées: {len(dropped)}")
    print(f"lignes après     : {len(kept)}")
    # Garde-fou : aucun rv ne doit rester en double, et on ne perd aucun gène.
    assert len({r["rv"] for r in kept}) == len(kept), "collision rv après dédup !"
    assert len(kept) == len(by_rv), "on a perdu un gène !"
    if dropped[:8]:
        print("\nexemples de lignes supprimées (rv, start, len_aa, score) :")
        for d in dropped[:8]:
            print("  ", d)

    if not write:
        print("\n[dry-run] relancer avec --write pour appliquer.")
        return

    backup = XREF.with_suffix(".tsv.predup")
    if not backup.exists():  # ne pas écraser une sauvegarde antérieure
        backup.write_text(XREF.read_text())
        print(f"\nsauvegarde -> {backup}")
    with XREF.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(kept)
    print(f"écrit -> {XREF} ({len(kept)} lignes)")


if __name__ == "__main__":
    main()
