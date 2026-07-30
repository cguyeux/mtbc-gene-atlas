#!/usr/bin/env python3
"""phase51_expression_imodulon.py -- P7.9 (volet expression conditionnelle) : couche iModulon.

Second volet de P7.9 (le GSMM est différé, redondant). L'atlas d'EXPRESSION CONDITIONNELLE de M.tb
sous forme interprétable = les **iModulons** (Yoo et al. 2022 mSphere, PMC9044949 ; compendium de 647
profils RNA-seq / 231 conditions décomposé par ICA en 80 gene sets co-régulés indépendamment, 41 =
régulons connus). Source machine-lisible : repo `Reosu/modulome_mtb` (iModulonDB).

Pour un gène, l'appartenance à un iModulon NOMMÉ (SigH, DosR, IniR, cholestérol...) = contexte de
CO-EXPRESSION conditionnelle interprétable. Pour un dark : co-régulé avec un ensemble sous une condition
donnée = handle (guilt-by-coexpression structuré). GARDE-FOU : co-expression ≠ fonction moléculaire →
couche de contexte (comme regulation/phenotype), PAS un levier de requalification.
Couche `expression` sur `Gene`. Run: python analyses/phase51_expression_imodulon.py
"""
from __future__ import annotations
import json, glob, csv, io, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "résultats" / "phase51_imodulon"
BASE = ("https://raw.githubusercontent.com/Reosu/modulome_mtb/HEAD/data/iModulonDB/"
        "organisms/m_tuberculosis/modulome/data_files")
REF = {"authors": "Yoo R, Rychel K, Poudel S, et al.", "year": 2022,
       "title": "Machine Learning of All Mycobacterium tuberculosis H37Rv RNA-seq Data Reveals a Structured "
                "Interplay between Metabolism, Stress Response, and Infection",
       "journal": "mSphere", "doi": "10.1128/msphere.00033-22"}


def fetch(name):
    return urllib.request.urlopen(f"{BASE}/{name}", timeout=60).read().decode()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # 1. table des iModulons : k -> {name, regulator, function, category}
    im = {}
    for row in csv.DictReader(io.StringIO(fetch("iM_table.csv"))):
        im[row["k"]] = {"name": row["Name"], "regulator": row.get("Regulator") or None,
                        "function": (row.get("Function") or "").strip(), "category": row.get("Category") or None}
    # 2. gene x iModulon (True/False)
    rows = list(csv.reader(io.StringIO(fetch("gene_presence_matrix.csv"))))
    header = rows[0][1:]   # indices d'iModulon
    gene_im = {}
    for r in rows[1:]:
        rv = r[0]
        ks = [header[i] for i, v in enumerate(r[1:]) if v.strip().lower() == "true"]
        if ks:
            gene_im[rv] = ks
    (OUT / "gene_imodulons.json").write_text(json.dumps(gene_im))
    print(f"iModulons : {len(im)} sets, {len(gene_im)} gènes avec ≥1 iModulon.")

    # 3. écrire la couche `expression` sur les fiches
    n = ndark = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        ks = gene_im.get(d["rv"])
        if not ks:
            continue
        mods = [{"name": im[k]["name"], "function": im[k]["function"],
                 "category": im[k]["category"], "regulator": im[k]["regulator"]}
                for k in ks if k in im]
        d["expression"] = {"n_imodulons": len(mods), "imodulons": mods,
                           "note": ("iModulon membership (independently-modulated gene sets from a 647-sample RNA-seq "
                                    "compendium): the conditional co-expression context. Co-expression is a regulatory "
                                    "context, NOT a molecular function."),
                           "source": "iModulonDB / modulome_mtb (Yoo 2022)", "reference": REF}
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
        if d.get("verdict") == "dark":
            ndark += 1
    print(f"Couche `expression` (iModulon) écrite sur {n} fiches, dont {ndark} dark.")


if __name__ == "__main__":
    main()
