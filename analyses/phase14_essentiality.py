#!/usr/bin/env python3
"""phase14_essentiality.py -- couche d'essentialité par mutagenèse de transposon (P3.4).

Ajoute un signal fonctionnel décisif et orthogonal au pN/pS : un gène "essentiel" en Tn-seq
est presque certainement un gène réel et fonctionnel (utile pour requalifier les "hypothetical").

Sources (déposées sous data/essentiality/, CC BY) :
  - DeJesus et al. 2017, mBio 8:e02133-16 (Himar1 saturé, H37Rv). Appel final par ORF :
    ES (essential) / ESD (essential domain) / GD (growth-defect) / GA (growth-advantage) /
    NE (non-essential) / Uncertain. Table S3.
  - Griffin et al. 2011, PLoS Pathog 7:e1002251. Table S4 : gènes spécifiquement requis pour la
    croissance sur cholestérol (signal métabolique distinctif du MTBC), utilisé ici comme drapeau.

Sortie : résultats/phase14_essentiality/essentiality.json, keyed par locus Rv :
  {rv: {"dejesus2017": "ES", "essential": true, "cholesterol_required": false, "refs": [...]}}
Stdlib + openpyxl.
Run: python analyses/phase14_essentiality.py
"""
from __future__ import annotations
import json, re
from pathlib import Path
import openpyxl

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "essentiality"
DEJESUS = DATA / "dejesus2017_tableS3_essentiality.xlsx"
GRIFFIN = DATA / "griffin2011_tableS4_cholesterol.xlsx"
OUT = ROOT / "résultats" / "phase14_essentiality"
GENES = ROOT / "site" / "content" / "genes"

RV = re.compile(r"^Rv\d{4}[A-Bc]?$")
# gènes "essentiels" au sens fonctionnel fort (requis pour la croissance in vitro)
ESSENTIAL_CALLS = {"ES", "ESD", "GD"}
REFS = [
    {"authors": "DeJesus MA, Gerrick ER, Xu W, et al.", "year": 2017,
     "title": "Comprehensive Essentiality Analysis of the Mycobacterium tuberculosis Genome via Saturating Transposon Mutagenesis",
     "journal": "mBio", "doi": "10.1128/mBio.02133-16"},
    {"authors": "Griffin JE, Gawronski JD, DeJesus MA, Ioerger TR, Akerley BJ, Sassetti CM", "year": 2011,
     "title": "High-resolution phenotypic profiling defines genes essential for mycobacterial growth and cholesterol catabolism",
     "journal": "PLoS Pathogens", "doi": "10.1371/journal.ppat.1002251"},
]


def load_dejesus():
    wb = openpyxl.load_workbook(DEJESUS, read_only=True, data_only=True)
    ws = wb[wb.sheetnames[0]]
    out = {}
    for r in ws.iter_rows(values_only=True):
        orf = r[0]
        if isinstance(orf, str) and RV.match(orf.strip()):
            call = str(r[12]).strip() if r[12] is not None else ""
            out[orf.strip()] = call
    return out


def load_griffin_cholesterol():
    wb = openpyxl.load_workbook(GRIFFIN, read_only=True, data_only=True)
    ws = wb[wb.sheetnames[0]]
    loci = set()
    for r in ws.iter_rows(values_only=True):
        for c in r:
            if isinstance(c, str) and RV.match(c.strip()):
                loci.add(c.strip())
    return loci


def main():
    dej = load_dejesus()
    chol = load_griffin_cholesterol()
    print(f"DeJesus 2017 : {len(dej)} ORF avec appel | Griffin cholestérol : {len(chol)} loci")

    rec = {}
    for rv, call in dej.items():
        rec[rv] = {
            "dejesus2017": call,
            "essential": call in ESSENTIAL_CALLS,
            "cholesterol_required": rv in chol,
            "refs": REFS,
        }
    # loci cholestérol absents de DeJesus (rare) : les ajouter quand même
    for rv in chol - set(dej):
        rec[rv] = {"dejesus2017": "", "essential": False, "cholesterol_required": True, "refs": REFS}

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "essentiality.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))
    print(f"Écrit {OUT/'essentiality.json'} ({len(rec)} loci)")

    # ── fusion directe dans les fiches existantes (idempotent, patron phase8/9) ──
    # (évite un re-run complet de phase4 qui écraserait les curations TA/dark_enzymes)
    import glob
    n_written = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        d["essentiality"] = rec.get(d["rv"], {})
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n_written += 1
    print(f"Fusionné 'essentiality' dans {n_written} fiches site/content/genes/")

    # ── stats + intersection avec les hypothétiques ──
    from collections import Counter
    print("Distribution DeJesus :", Counter(v["dejesus2017"] for v in rec.values()).most_common())
    print("Essentiels (ES/ESD/GD) :", sum(1 for v in rec.values() if v["essential"]))
    if GENES.exists():
        import glob
        hyp = {}
        for f in glob.glob(str(GENES / "*.json")):
            d = json.load(open(f))
            if "hypothetical" in (d.get("product_h37rv") or "").lower():
                hyp[d["rv"]] = d
        h_ess = sum(1 for rv in hyp if rec.get(rv, {}).get("essential"))
        h_es = sum(1 for rv in hyp if rec.get(rv, {}).get("dejesus2017") == "ES")
        print(f"Hypothétiques : {len(hyp)} | essentiels (ES/ESD/GD) : {h_ess} | strictement ES : {h_es}")
        print(f"Couverture DeJesus sur les 3906 fiches : {sum(1 for f in glob.glob(str(GENES/'*.json')) if json.load(open(f))['rv'] in rec)}")


if __name__ == "__main__":
    main()
