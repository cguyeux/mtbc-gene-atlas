#!/usr/bin/env python3
"""phase59_synteny_mtbc_specific.py -- P10.2c : confirmation par synténie des candidats spécifiques-MTBC.

Garde-fou (KB) : un gène classé « MTBC-specific » par tblastn (absent des 53 NTM) peut être (a) une vraie
INNOVATION insérée à un locus conservé, (b) un élément MOBILE (IS/prophage/PE-PGRS), (c) une région
plus large perdue chez les NTM. On distingue par la SYNTÉNIE : si les deux voisins génomiques immédiats
sont conservés à travers le genre (présents dans les NTM) mais que le gène est absent, c'est une insertion
propre au MTBC à un locus conservé = candidat innovation SOLIDE. Sinon on rétrograde.

Entrées : couche `outgroup.ntm` (phase58) + `genomic_context.left/right` + `funccat` de chaque fiche.
Sortie : résultats/phase59_synteny/mtbc_specific_confirmed.tsv (+ md).
Run: python analyses/phase59_synteny_mtbc_specific.py
"""
from __future__ import annotations
import json, glob, re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "résultats" / "phase59_synteny"
FLANK_MIN_NTM = 10           # un voisin est "conservé" s'il est présent dans >=10 NTM
MOBILE_RE = re.compile(r"transpos|integrase|phage|insertion seq|IS\d|PE.?PGRS|PPE|resolvase|recombinase", re.I)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    G = {json.load(open(f))["rv"]: json.load(open(f)) for f in glob.glob(str(GENES / "*.json"))}

    def ntm_cls(rv):
        return ((G.get(rv, {}).get("outgroup") or {}).get("ntm") or {}).get("classification")

    def ntm_n(rv):
        return ((G.get(rv, {}).get("outgroup") or {}).get("ntm") or {}).get("n_present", 0)

    def is_mobile(d):
        cat = (d.get("funccat") or {}).get("category") or ""
        txt = f"{cat} {d.get('product_h37rv','')} {d.get('gene_name') or ''} {d.get('function_revised','')}"
        return "insertion seq" in cat.lower() or bool(MOBILE_RE.search(txt))

    rows = []
    for rv, d in G.items():
        cls = ntm_cls(rv)
        if cls not in ("MTBC-specific", "MTBC-near-specific"):
            continue
        gc = d.get("genomic_context") or {}
        left = (gc.get("left") or {}).get("locus")
        right = (gc.get("right") or {}).get("locus")
        lc = ntm_n(left) >= FLANK_MIN_NTM if left else False
        rc = ntm_n(right) >= FLANK_MIN_NTM if right else False
        mob = is_mobile(d)
        rd = bool(d.get("rd"))
        if mob:
            verdict = "mobile-element"           # IS/phage/PE-PGRS : innovation non retenue (contexte mobile)
        elif lc and rc:
            verdict = "confirmed-insertion"       # voisins conservés + gène absent = insertion propre au MTBC SOLIDE
        elif not lc and not rc:
            verdict = "regional-island"           # voisins aussi absents = région MTBC-specific plus large
        else:
            verdict = "partial-synteny"           # un seul voisin conservé
        rows.append({"rv": rv, "ntm_class": cls, "verdict": verdict, "dark": d.get("verdict") == "dark",
                     "len": len(d.get("protein_mtbc0") or ""), "left": left, "left_n": ntm_n(left),
                     "right": right, "right_n": ntm_n(right), "mobile": mob, "rd": rd,
                     "ess": bool((d.get("essentiality") or {}).get("essential")),
                     "vi": (d.get("vulnerability") or {}).get("vi"),
                     "prod": (d.get("product_h37rv") or "")[:45]})

    rows.sort(key=lambda r: (r["verdict"] != "confirmed-insertion", not r["dark"], -r["len"]))
    with open(OUT / "mtbc_specific_confirmed.tsv", "w") as fh:
        fh.write("rv\tntm_class\tsynteny_verdict\tdark\tlen\tleft\tleft_ntm\tright\tright_ntm\tmobile\trd\tess\tvi\tproduct\n")
        for r in rows:
            fh.write(f"{r['rv']}\t{r['ntm_class']}\t{r['verdict']}\t{int(r['dark'])}\t{r['len']}\t{r['left']}\t"
                     f"{r['left_n']}\t{r['right']}\t{r['right_n']}\t{int(r['mobile'])}\t{int(r['rd'])}\t"
                     f"{int(r['ess'])}\t{r['vi'] if r['vi'] is not None else ''}\t{r['prod']}\n")

    vc = Counter(r["verdict"] for r in rows)
    dark_conf = [r for r in rows if r["dark"] and r["verdict"] == "confirmed-insertion"]
    print(f"{len(rows)} MTBC-specific/near total. Verdicts synténie : {dict(vc)}")
    print(f"DARK + confirmed-insertion (innovations MTBC solides, inconnues) : {len(dark_conf)}")
    for r in sorted(dark_conf, key=lambda r: -r["len"])[:12]:
        print(f"  {r['rv']} ({r['len']}aa) flanks {r['left']}[{r['left_n']}]/{r['right']}[{r['right_n']}] "
              f"{'ESS' if r['ess'] else ''} :: {r['prod']}")
    print(f"Rapport : {OUT/'mtbc_specific_confirmed.tsv'}")


if __name__ == "__main__":
    main()
