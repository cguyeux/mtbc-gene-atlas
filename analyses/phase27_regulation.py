#!/usr/bin/env python3
"""phase27_regulation.py -- régulation transcriptionnelle / réseau de TF (P5.6).

Champ de parité Mycobrowser (contexte régulateur) : par gène, les facteurs de transcription
(TF) qui le régulent (activation/répression), et si le gène EST un TF, la taille de son régulon.

Source : « MTB Signed TRN » du MTB Network Portal (Institute for Systems Biology), edge list
curée `TF | Target | Sign` (±1) dérivée de la ChIP-seq de 143 TF (Minch et al. 2015) + du
réseau d'induction par surexpression de TF (TFOE, Rustad et al. 2014). Locus tags Rv, jointure
directe. Fichier `data/tfoe/MTB-Signed-TRN.xlsx` (converti du .xls d'origine via libreoffice).

Sortie : résultats/phase27_regulation/regulation.json (keyed Rv). Fusion dans les fiches
(idempotent) + enregistrement phase4.
Run: python analyses/phase27_regulation.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path
from collections import defaultdict
import openpyxl

ROOT = Path(__file__).resolve().parent.parent
TRN = ROOT / "data" / "tfoe" / "MTB-Signed-TRN.xlsx"
OUT = ROOT / "résultats" / "phase27_regulation"
GENES = ROOT / "site" / "content" / "genes"
CAP = 15  # régulateurs affichés max par gène

SOURCE = "MTB signed transcriptional regulatory network (ISB MTB Network Portal)"
REFS = [
    {"authors": "Minch KJ, Rustad TR, Peterson EJR, et al.", "year": 2015,
     "title": "The DNA-binding network of Mycobacterium tuberculosis",
     "journal": "Nature Communications", "doi": "10.1038/ncomms6829"},
    {"authors": "Rustad TR, Minch KJ, Ma S, et al.", "year": 2014,
     "title": "Mapping and manipulating the Mycobacterium tuberculosis transcriptome using a transcription factor overexpression-derived regulatory network",
     "journal": "Genome Biology", "doi": "10.1186/gb-2014-15-11-502"},
]


def main():
    print("== phase27 : régulation transcriptionnelle / TF (P5.6) ==")
    # noms de gènes rv -> gene
    rv2name = {}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        if d.get("gene"):
            rv2name[d["rv"]] = d["gene"]

    ws = openpyxl.load_workbook(TRN, read_only=True)[0 if False else "MTB Signed TRN"]
    regulators = defaultdict(list)   # target -> [(tf, sign)]
    targets = defaultdict(list)      # tf -> [(target, sign)]
    n_edges = 0
    for i, row in enumerate(ws.iter_rows(values_only=True)):
        if i == 0:
            continue  # header TF/Target/Sign
        tf, tgt, sign = (str(row[0]).strip() if row[0] else ""), (str(row[1]).strip() if row[1] else ""), row[2]
        if not tf or not tgt:
            continue
        try:
            s = int(float(str(sign)))
        except (ValueError, TypeError):
            s = 0
        regulators[tgt].append((tf, s))
        targets[tf].append((tgt, s))
        n_edges += 1
    print(f"{n_edges} arêtes | {len(targets)} TF (régulateurs) | {len(regulators)} gènes cibles")

    def eff(s):
        return "activates" if s > 0 else ("represses" if s < 0 else "regulates")

    rec = {}
    all_rv = set(regulators) | set(targets)
    for rv in all_rv:
        regs = regulators.get(rv, [])
        reg_objs = [{"tf": t, "tf_gene": rv2name.get(t), "effect": eff(s)} for t, s in regs]
        is_tf = rv in targets
        rec[rv] = {
            "regulated_by": reg_objs[:CAP],
            "n_regulators": len(regs),
            "is_tf": is_tf,
            "regulon_size": len(targets.get(rv, [])),
            "source": SOURCE, "refs": REFS,
        }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "regulation.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))
    print(f"Écrit {OUT/'regulation.json'} ({len(rec)} loci)")

    # fusion dans les fiches
    n_written = n_reg = n_tf = hyp_reg = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        r = rec.get(d["rv"])
        d["regulation"] = r or {}
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n_written += 1
        if r and (r["n_regulators"] or r["is_tf"]):
            n_reg += 1
            if r["is_tf"]:
                n_tf += 1
            if r["n_regulators"] and "hypothetical" in (d.get("product_h37rv") or "").lower():
                hyp_reg += 1
    print(f"Fusionné 'regulation' dans {n_written} fiches ({n_reg} avec régulation, {n_tf} TF, "
          f"{hyp_reg} hypothétiques régulés)")


if __name__ == "__main__":
    main()
