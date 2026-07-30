#!/usr/bin/env python3
"""phase22_funccat.py -- catégorie fonctionnelle TubercuList/Mycobrowser (P5.4).

Champ de parité Mycobrowser : les 11 catégories fonctionnelles communautaires du schéma
TubercuList (Cole et al. 1998, maintenu par Mycobrowser). C'est le classement fonctionnel
de haut niveau familier aux tuberculologues (« cell wall and cell processes », « lipid
metabolism », « virulence, detoxification, adaptation », « conserved hypotheticals »...).

Source : release 5 de Mycobrowser (fichier plat `txt` H37Rv, colonne `Functional_Category`),
téléchargée dans data/mycobrowser/H37Rv.txt. Keyed par locus (= locus tag Rv), jointure
directe sur nos fiches. Curation de la base de RÉFÉRENCE reprise telle quelle et ATTRIBUÉE
(pas une valeur qu'on requalifie ; c'est le champ de parité lui-même).

Sortie : résultats/phase22_funccat/funccat.json (keyed Rv). Fusion directe dans les fiches
(idempotent) + enregistrement dans phase4.
Run: python analyses/phase22_funccat.py
"""
from __future__ import annotations
import json, csv, glob
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent
MYCO = ROOT / "data" / "mycobrowser" / "H37Rv.txt"
OUT = ROOT / "résultats" / "phase22_funccat"
GENES = ROOT / "site" / "content" / "genes"

SOURCE = "TubercuList functional category (via Mycobrowser release 5)"
REFS = [
    {"authors": "Cole ST, Brosch R, Parkhill J, et al.", "year": 1998,
     "title": "Deciphering the biology of Mycobacterium tuberculosis from the complete genome sequence",
     "journal": "Nature", "doi": "10.1038/31159"},
    {"authors": "Kapopoulou A, Lew JM, Cole ST", "year": 2011,
     "title": "The MycoBrowser portal: a comprehensive and manually annotated resource for mycobacterial genomes",
     "journal": "Tuberculosis", "doi": "10.1016/j.tube.2010.09.006"},
]


def main():
    print("== phase22 : catégorie fonctionnelle TubercuList (P5.4) ==")
    cat_by_locus = {}
    with open(MYCO, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            loc = (row.get("Locus") or "").strip()
            fc = (row.get("Functional_Category") or "").strip()
            if loc and fc:
                cat_by_locus[loc] = fc
    print(f"Mycobrowser : {len(cat_by_locus)} loci avec catégorie")
    print("Distribution :", Counter(cat_by_locus.values()).most_common())

    rec = {}
    n_written = n_hit = 0
    miss = []
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        rv = d["rv"]
        fc = cat_by_locus.get(rv)
        if fc:
            rec[rv] = {"category": fc, "source": SOURCE, "refs": REFS}
            d["funccat"] = rec[rv]
            n_hit += 1
        else:
            d["funccat"] = {}
            miss.append(rv)
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n_written += 1

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "funccat.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))
    print(f"Écrit {OUT/'funccat.json'} ({len(rec)} catégories)")
    print(f"Fusionné 'funccat' dans {n_written} fiches ({n_hit} avec catégorie, {len(miss)} sans)")
    if miss:
        print(f"Sans catégorie (ex.) : {miss[:10]}")


if __name__ == "__main__":
    main()
