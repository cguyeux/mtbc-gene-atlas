#!/usr/bin/env python3
"""phase24_mutant_phenotypes.py -- phénotypes de mutants conditionnels (P5.2).

Champ de parité Mycobrowser AU-DELÀ de l'essentialité binaire (P3.4) : la CONSÉQUENCE de la
disruption d'un gène selon la CONDITION (in vivo souris, cholestérol, stress acide/NO/fer,
carence, exposition médicament...). Un gène non essentiel in vitro peut être requis in vivo
(facteur de virulence) — c'est ce profil conditionnel qu'on ajoute.

Source : MtbTnDB (Jinich et al. 2025, Mol Microbiol), compendium standardisé de ~146 conditions
Tn-seq agrégées et ré-analysées uniformément (repo ajinich/mtb_tn_db_demo). Format long
`Rv_ID / Expt / log2FC / q-val` : log2FC<0 significatif = mutant APPAUVRI = gène requis dans
la condition ; log2FC>0 = disruption avantageuse. Métadonnées de condition dans col_desc.tsv.

Garde-fous : (1) on EXCLUT les 22 comparaisons d'interaction génétique (KO-vs-KO, épistasie,
pas phénotype de délétion simple). (2) Seuil q<=0.05 ET |log2FC|>=1 (2x). (3) « disruption »
(insertion Tn), pas délétion propre — libellé prudent.

Sortie : résultats/phase24_mutant_phenotypes/mutant_phenotypes.json (keyed Rv). Fusion directe
dans les fiches (idempotent) + enregistrement dans phase4. Stdlib pur.
Run: python analyses/phase24_mutant_phenotypes.py
"""
from __future__ import annotations
import json, csv, glob
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent
MTB = ROOT / "data" / "mtbtndb"
OUT = ROOT / "résultats" / "phase24_mutant_phenotypes"
GENES = ROOT / "site" / "content" / "genes"

QMAX = 0.05
LFC_MIN = 1.0
CAP = 12  # phénotypes affichés max par gène

DRUGS = {"isoniazid", "ethambutol", "rifampicin", "vancomycin", "streptomycin", "moxifloxacin",
         "bedaquiline", "linezolid", "meropenem", "rifabutin"}
SOURCE = "MtbTnDB (standardized Tn-seq compendium, ~146 conditions)"
REFS = [
    {"authors": "Jinich A, Zaveri A, DeJesus MA, et al.", "year": 2025,
     "title": "The Mycobacterium tuberculosis Transposon Sequencing Database (MtbTnDB): A Large-Scale Guide to Genetic Conditional Essentiality",
     "journal": "Molecular Microbiology", "doi": "10.1111/mmi.15370"},
]


def build_condition_map():
    """column_ID_std -> {label, category, is_gi} pour les conditions interprétables."""
    cmap = {}
    for r in csv.DictReader(open(MTB / "col_desc.tsv", encoding="utf-8"), delimiter="\t"):
        cid = (r.get("column_ID_std") or "").strip()
        if not cid:
            continue
        is_gi = bool((r.get("GI_RvID") or "").strip())
        meaning = (r.get("meaning") or "").strip()
        invivo = (r.get("in vitro/cell/in vivo") or "").strip()
        stress = (r.get("stress") or "").strip()
        carbon = (r.get("carbon source") or "").strip().lower()
        mouse = (r.get("mouse strain") or "").strip()
        cell = (r.get("cell type") or "").strip()
        low = cid.lower()

        cat = label = None
        if stress and stress != "-":
            if stress.lower() in DRUGS:
                cat, label = "drug exposure", f"altered fitness under {stress}"
            else:
                cat, label = "stress", f"altered fitness under {stress}"
        elif carbon == "cholesterol":
            cat, label = "carbon source", "fitness on cholesterol (vs glycerol)"
        elif invivo == "in_vivo" and (mouse or "mouse" in low or "vivo" in low):
            det = ""
            if "d45" in low: det = ", day 45"
            elif "d10" in low: det = ", day 10"
            if "mhcii" in low: det = ", immunodeficient (MHC-II-/-)" + det
            cat, label = "in vivo", f"fitness in mouse infection{det}"
        elif invivo == "cell" or cell:
            cat, label = "macrophage", f"fitness in {cell or 'macrophage'}"
        elif "day" in low and "vs" in low:
            cat, label = "in vitro passage", "fitness after prolonged in vitro passage"
        elif meaning:
            cat, label = "other", meaning
        # sinon : comparaison de normalisation/batch → non classée, ignorée

        if cat:
            cmap[cid] = {"label": label, "category": cat, "is_gi": is_gi}
    return cmap


def main():
    print("== phase24 : phénotypes de mutants conditionnels MtbTnDB (P5.2) ==")
    cmap = build_condition_map()
    usable = {k: v for k, v in cmap.items() if not v["is_gi"]}
    print(f"Conditions interprétables : {len(cmap)} (dont {len(cmap)-len(usable)} interactions génétiques exclues) "
          f"-> {len(usable)} utilisées")

    # parcours du fichier long, groupé par gène
    per_gene = defaultdict(list)
    n_rows = 0
    with open(MTB / "standardized.tsv", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n_rows += 1
            expt = row["Expt"]
            cond = usable.get(expt)
            if not cond:
                continue
            try:
                lf = float(row["log2FC"]); q = float(row["q-val"])
            except (ValueError, KeyError):
                continue
            if q <= QMAX and abs(lf) >= LFC_MIN:
                per_gene[row["Rv_ID"]].append({
                    "condition": cond["label"], "category": cond["category"],
                    "log2fc": round(lf, 2), "qval": float(f"{q:.2g}"),
                    "effect": "required (mutant depleted)" if lf < 0 else "disruption advantageous (mutant enriched)",
                })
    print(f"Parcouru {n_rows} lignes (gène×condition) ; {len(per_gene)} gènes avec >=1 phénotype significatif")

    rec = {}
    for rv, phs in per_gene.items():
        phs.sort(key=lambda p: -abs(p["log2fc"]))
        n_invivo = sum(1 for p in phs if p["category"] in ("in vivo", "macrophage"))
        rec[rv] = {
            "n_significant": len(phs),
            "n_in_vivo": n_invivo,
            "phenotypes": phs[:CAP],
            "source": SOURCE,
            "refs": REFS,
        }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "mutant_phenotypes.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))
    print(f"Écrit {OUT/'mutant_phenotypes.json'} ({len(rec)} gènes)")

    # fusion dans les fiches
    n_written = n_hit = hyp_hit = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        d["mutant_phenotypes"] = rec.get(d["rv"], {})
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n_written += 1
        if d["rv"] in rec:
            n_hit += 1
            if "hypothetical" in (d.get("product_h37rv") or "").lower():
                hyp_hit += 1
    print(f"Fusionné 'mutant_phenotypes' dans {n_written} fiches ({n_hit} avec phénotype, dont {hyp_hit} hypothétiques)")

    # bilan
    from collections import Counter
    cat = Counter()
    for r in rec.values():
        for p in r["phenotypes"]:
            cat[p["category"]] += 1
    print("Phénotypes affichés par catégorie :", cat.most_common())
    invivo_genes = sum(1 for r in rec.values() if r["n_in_vivo"] > 0)
    print(f"Gènes avec >=1 phénotype IN VIVO (virulence) : {invivo_genes}")


if __name__ == "__main__":
    main()
