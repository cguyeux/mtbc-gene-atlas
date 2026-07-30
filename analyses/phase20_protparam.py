#!/usr/bin/env python3
"""phase20_protparam.py -- propriétés physico-chimiques par protéine (P5.3).

Champ de parité Mycobrowser trivial mais utile : poids moléculaire, point isoélectrique
théorique, GRAVY (hydropathie moyenne), indice d'instabilité, aromaticité, indice
aliphatique. Purement calculable depuis la séquence (Biopython ProtParam, méthode ExPASy
de Gasteiger et al. 2005). Sert de contexte de base sur chaque fiche et, pour les
« hypothetical », donne des indices grossiers (un GRAVY très positif + pas de TM DeepTMHMM
= incohérence à regarder ; un pI extrême oriente la fonction).

Séquence utilisée : `protein_mtbc0` (protéine ancestrale MTBC0 de la fiche). Résidus non
standard (X, U/sélénocystéine, *, B, Z, O, J) retirés avant calcul (patron anti-crash déjà
appris sur phase2i/BLOSUM) ; la longueur rapportée reste celle de la séquence d'origine.

Sortie : résultats/phase20_protparam/protparam.json (keyed Rv). Fusion directe dans les
fiches (patron phase8/9/14/19, idempotent) + enregistrement dans phase4.
Run: python analyses/phase20_protparam.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "phase20_protparam"
GENES = ROOT / "site" / "content" / "genes"

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
SOURCE = "Biopython ProtParam (méthode ExPASy ProtParam)"
REFS = [
    {"authors": "Gasteiger E, Hoogland C, Gattiker A, et al.", "year": 2005,
     "title": "Protein Identification and Analysis Tools on the ExPASy Server",
     "journal": "The Proteomics Protocols Handbook (Humana Press)", "doi": "10.1385/1-59259-890-0:571"},
]


def aliphatic_index(aa_percent: dict) -> float:
    """Indice aliphatique d'Ikai (1980) : AI = %Ala + 2.9*%Val + 3.9*(%Ile + %Leu)."""
    a = aa_percent.get("A", 0) * 100
    v = aa_percent.get("V", 0) * 100
    i = aa_percent.get("I", 0) * 100
    l = aa_percent.get("L", 0) * 100
    return round(a + 2.9 * v + 3.9 * (i + l), 1)


def analyse(seq: str) -> dict | None:
    raw_len = len(seq)
    clean = "".join(c for c in seq.upper() if c in STANDARD)
    if len(clean) < 5:
        return None
    pa = ProteinAnalysis(clean)
    inst = pa.instability_index()
    return {
        "length_aa": raw_len,
        "mw_kda": round(pa.molecular_weight() / 1000.0, 1),
        "pi": round(pa.isoelectric_point(), 2),
        "gravy": round(pa.gravy(), 3),
        "instability_index": round(inst, 1),
        "instability_class": "unstable" if inst > 40 else "stable",
        "aromaticity": round(pa.aromaticity(), 3),
        "aliphatic_index": aliphatic_index(pa.amino_acids_percent),
        "source": SOURCE,
        "refs": REFS,
    }


def main():
    print("== phase20 : propriétés physico-chimiques ProtParam (P5.3) ==")
    rec = {}
    n_written = n_ok = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        seq = (d.get("protein_mtbc0") or "").strip()
        pp = analyse(seq) if seq else None
        d["protparam"] = pp or {}
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n_written += 1
        if pp:
            rec[d["rv"]] = pp
            n_ok += 1

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "protparam.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))
    print(f"Écrit {OUT/'protparam.json'} ({n_ok} protéines)")
    print(f"Fusionné 'protparam' dans {n_written} fiches ({n_ok} calculées)")

    # ── bilan distributionnel (sanity) ──
    import statistics as st
    mw = [v["mw_kda"] for v in rec.values()]
    pi = [v["pi"] for v in rec.values()]
    gr = [v["gravy"] for v in rec.values()]
    unstable = sum(1 for v in rec.values() if v["instability_class"] == "unstable")
    print(f"MW kDa   : médiane {st.median(mw):.1f} (min {min(mw):.1f}, max {max(mw):.1f})")
    print(f"pI       : médiane {st.median(pi):.2f} (acides <7 : {sum(1 for x in pi if x<7)}, basiques >7 : {sum(1 for x in pi if x>7)})")
    print(f"GRAVY    : médiane {st.median(gr):.3f} (hydrophobes >0 : {sum(1 for x in gr if x>0)})")
    print(f"Instables (indice >40) : {unstable}/{len(rec)} ({100*unstable/len(rec):.0f} %)")


if __name__ == "__main__":
    main()
