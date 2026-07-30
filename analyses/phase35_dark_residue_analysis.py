#!/usr/bin/env python3
"""phase35_dark_residue_analysis.py -- P7.1d + P7.1e (post-HHpred).

P7.1d : corrélation taille d'ORF vs confiance HHpred sur les 28 gènes du guide
        (best probability parsée des top_hits de phase34, verdict, longueur).
P7.1e : audit des mentions "coiled-coil" dans TOUT l'atlas (site/content/genes/*.json)
        pour repérer d'éventuelles requalifications sur ce seul signal (garde-fou anti-artefact,
        cf. KB tuberculosis.md 2026-06-10 + rétrogradation Rv3430a 2026-07-04).
Run: python analyses/phase35_dark_residue_analysis.py
"""
from __future__ import annotations
import json, re, glob, statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
GUIDE = ROOT / "data" / "dark_residue" / "hhpred_guide.md"

# import du dict HITS de phase34 sans exécuter son main()
import importlib.util
spec = importlib.util.spec_from_file_location("phase34", ROOT / "analyses" / "phase34_hhpred_curation.py")
assert spec is not None and spec.loader is not None
phase34 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(phase34)
HITS = phase34.HITS


def best_prob(entry):
    """meilleure Probability HHpred parsée du 1er top_hit (ex. '... (80.5%, E=4.3)')."""
    hits = entry.get("top_hits") or []
    probs = []
    for h in hits:
        for m in re.findall(r"(\d+(?:\.\d+)?)\s*%", h):
            probs.append(float(m))
    return max(probs) if probs else None


def parse_guide_lengths():
    """{rv: length_aa} depuis les headers '>RvXXXX NNaa' du guide."""
    out = {}
    for line in open(GUIDE):
        if line.startswith(">"):
            m = re.match(r">(\S+)\s+(\d+)aa", line)
            if m:
                out[m.group(1)] = int(m.group(2))
    return out


def p71d():
    print("=" * 70)
    print("P7.1d — taille d'ORF vs confiance HHpred (28 gènes du guide)")
    print("=" * 70)
    lengths = parse_guide_lengths()
    rows = []
    for rv, entry in HITS.items():
        try:
            d = json.load(open(GENES / f"{rv}.json"))
        except FileNotFoundError:
            continue
        verdict = d.get("verdict")
        requalified = verdict in ("requalified", "family_assigned")
        rows.append({
            "rv": rv,
            "len": lengths.get(rv),
            "best_prob": best_prob(entry),
            "requalified": requalified,
        })
    rows = [r for r in rows if r["len"] is not None]
    rows.sort(key=lambda r: r["len"])

    small = [r for r in rows if r["len"] < 65]
    large = [r for r in rows if r["len"] >= 65]

    def summ(g, label):
        if not g:
            return
        n_req = sum(r["requalified"] for r in g)
        probs = [r["best_prob"] for r in g if r["best_prob"] is not None]
        med = statistics.median(probs) if probs else float("nan")
        print(f"  {label:<18} n={len(g):<3} requalifiés={n_req} ({100*n_req/len(g):.0f}%)  "
              f"best-prob médiane={med:.1f}%")

    print(f"\nEnsemble du guide : {len(rows)} gènes, "
          f"{sum(r['requalified'] for r in rows)} requalifiés")
    summ(small, "micro-ORF <65 aa")
    summ(large, "ORF >=65 aa")

    print("\n  Détail (trié par taille) :")
    print(f"  {'gène':<9}{'len':>5}{'best%':>8}  verdict")
    for r in rows:
        bp = f"{r['best_prob']:.1f}" if r["best_prob"] is not None else "-"
        print(f"  {r['rv']:<9}{r['len']:>5}{bp:>8}  {'REQUAL' if r['requalified'] else 'dark'}")

    # corrélation length vs requalification
    if len(rows) > 2:
        import math
        xs = [r["len"] for r in rows]
        ys = [1 if r["requalified"] else 0 for r in rows]
        mx, my = statistics.mean(xs), statistics.mean(ys)
        num = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
        den = math.sqrt(sum((x-mx)**2 for x in xs) * sum((y-my)**2 for y in ys))
        r_pb = num/den if den else float("nan")
        print(f"\n  Corrélation point-bisériale longueur↔requalification : r = {r_pb:.2f}")


def p71e():
    print("\n" + "=" * 70)
    print("P7.1e — audit des mentions 'coiled-coil' dans l'atlas complet")
    print("=" * 70)
    hits = []
    for f in glob.glob(str(GENES / "*.json")):
        raw = open(f).read()
        if re.search(r"coiled[- ]coil", raw, re.I):
            d = json.loads(raw)
            fields = []
            for k in ("product", "function_revised"):
                v = d.get(k)
                if isinstance(v, str) and re.search(r"coiled[- ]coil", v, re.I):
                    fields.append(k)
            # champs imbriqués (pfam, hhpred, string...)
            for k in ("pfam", "pfam_tentative", "hhpred", "eggnog", "uniprot"):
                if re.search(r"coiled[- ]coil", json.dumps(d.get(k) or ""), re.I):
                    fields.append(k)
            hits.append({"rv": d["rv"], "verdict": d.get("verdict"),
                         "auto": d.get("auto"), "fields": sorted(set(fields))})
    hits.sort(key=lambda h: (h["verdict"] != "requalified", h["rv"]))
    print(f"\n  {len(hits)} fiche(s) mentionnent un coiled-coil.")
    # celles requalifiées où le coiled-coil apparaît dans function_revised = à surveiller
    suspect = [h for h in hits if h["verdict"] in ("requalified", "family_assigned")
               and "function_revised" in h["fields"]]
    print(f"  Requalifiées AVEC coiled-coil dans function_revised (à re-vérifier) : {len(suspect)}")
    for h in hits:
        flag = "  <-- REQUAL sur function_revised" if h in suspect else ""
        print(f"    {h['rv']:<9} verdict={h['verdict']:<14} champs={h['fields']}{flag}")
    if not suspect:
        print("\n  => aucune requalification reposant sur un coiled-coil en function_revised "
              "(Rv3430a déjà rétrogradé). Garde-fou anti-artefact respecté.")


if __name__ == "__main__":
    p71d()
    p71e()
