#!/usr/bin/env python3
"""phase45_vaccine_shortlist.py -- P7.5 (cadrage vaccin) : shortlist de candidats vaccins.

NON une requalification (le levier IEDB de requalification dark est épuisé, cf. cahier). Ici =
LIVRABLE de VALORISATION / hand-off labo : croiser les hypothétiques SÉCRÉTÉS porteurs d'un lead
phénotype IN VIVO (P5.2b) avec l'antigénicité IEDB, pour proposer des candidats vaccins.

Profil vaccin TB idéal : SÉCRÉTÉ (accessible à l'immunité) + REQUIS IN VIVO (virulence/persistence,
cible pertinente) + CONSERVÉ (une cible diversifiante = échappement immunitaire) + ANTIGÉNIQUE
(épitopes T mesurés, IEDB). Score = priorité phénotype + bonus conservation (strong purifying) +
bonus antigénicité ; PÉNALITÉ diversifying (mauvais candidat vaccin).

Contrôle positif attendu : Rv1860 (Apa/FbpD, antigène sécrété immunodominant connu, ~61 épitopes T)
doit sortir en tête → valide le scoring. Le VRAI gisement = les candidats fortement conservés SANS
antigénicité connue (sous-explorés, à tester).
Run: python analyses/phase45_vaccine_shortlist.py
"""
from __future__ import annotations
import json, glob, urllib.request, urllib.parse
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "résultats" / "phase45_vaccine_candidates"


def iedb_tcell(acc):
    if not acc:
        return 0
    iri = urllib.parse.quote(f"UNIPROT:{acc}", safe='')
    url = (f"https://query-api.iedb.org/tcell_search?parent_source_antigen_iri=eq.{iri}"
           f"&select=structure_id,qualitative_measure")
    try:
        rows = json.loads(urllib.request.urlopen(url, timeout=25).read().decode())
        pos = [x for x in rows if 'positive' in (x.get('qualitative_measure') or '').lower()]
        return len(set(x['structure_id'] for x in pos))
    except Exception:
        return -1


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cands = []
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        if (d.get("localization") or {}).get("type") not in ("SP", "SP+TM"):
            continue
        pl = d.get("phenotype_lead") or {}
        lead = pl.get("lead") or "" if isinstance(pl, dict) else ""
        if not any(k in lead.lower() for k in ("in vivo", "virulence", "persistence")):
            continue
        fc = d.get("funccat") or {}
        cat = fc.get("category") if isinstance(fc, dict) else None
        if not (d.get("is_hypothetical") or cat == "conserved hypotheticals"
                or d.get("verdict") in ("dark", "family_assigned")):
            continue
        acc = (d.get("uniprot") or {}).get("acc")
        cons = (d.get("conservation") or {}).get("selection", "")
        nepi = iedb_tcell(acc)
        s = (pl.get("priority") or 0)
        if "strong purifying" in cons:
            s += 2
        elif "purifying" in cons:
            s += 1
        elif "diversifying" in cons:
            s -= 2
        if nepi > 0:
            s += 2 + min(nepi, 3)
        cands.append({"rv": d["rv"], "loc": (d.get("localization") or {}).get("type"),
                      "acc": acc, "cons": cons, "prio": pl.get("priority"),
                      "nepi": nepi, "verdict": d.get("verdict"), "score": round(s, 1),
                      "lead": lead[:80]})
    cands.sort(key=lambda c: -c["score"])

    lines = ["# Candidats vaccins MTBC — hypothétiques sécrétés, requis in vivo (P7.5)", "",
             "Livrable hand-off labo (NON une requalification). Profil : sécrété + requis in vivo + conservé + antigénique.",
             "Score = priorité phénotype (P5.2b) + conservation (strong purifying +2) + antigénicité IEDB (+2..5) ; diversifying -2.",
             "Contrôle positif : Rv1860 (Apa/FbpD, antigène connu) doit ressortir en tête.", "",
             "| score | gène | loc | conservation | prio phéno | épitopes T (IEDB) | verdict |",
             "|---|---|---|---|---|---|---|"]
    for c in cands:
        ep = c["nepi"] if c["nepi"] >= 0 else "err"
        lines.append(f"| {c['score']} | {c['rv']} | {c['loc']} | {c['cons']} | {c['prio']} | {ep} | {c['verdict']} |")
    # sous-explorés = fortement conservés, requis in vivo, SANS antigénicité connue (le vrai gisement)
    subexpl = [c for c in cands if c["nepi"] == 0 and "purifying" in c["cons"]]
    lines += ["", "## Candidats SOUS-EXPLORÉS (conservés, in vivo, PAS d'antigénicité connue → à tester)",
              ", ".join(c["rv"] for c in subexpl)]
    (OUT / "shortlist.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"{len(cands)} candidats. Top: {cands[0]['rv']} (score {cands[0]['score']}, {cands[0]['nepi']} épitopes) = contrôle.")
    print(f"Sous-explorés (à tester): {len(subexpl)} -> {[c['rv'] for c in subexpl]}")
    print(f"Rapport: {OUT/'shortlist.md'}")


if __name__ == "__main__":
    main()
