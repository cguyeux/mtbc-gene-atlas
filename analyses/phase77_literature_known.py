#!/usr/bin/env python3
"""phase77_literature_known.py -- P16.2a-quater : couche `literature` pour les 2795 gènes CONNUS + « étudié ailleurs ».

Complète phase70 (219 dark) et phase75 (892 hypothétiques non-dark) : après ce script, **tout le protéome
(3906 gènes) est ancré dans la littérature**. Source : `phase76_literature_proteome.py`
(`résultats/phase76_literature/proteome_literature.tsv`).

DEUX APPORTS SUR LA FICHE :
 1. **Citations** : les publications qui discutent le gène (recherche par tag Rv + NOM DE GÈNE + identifiants
    d'orthologues, sous filtre de contexte mycobactérien, vérifiée sur le texte quand le compte est petit).
    Rappel du piège : la littérature n'appelle JAMAIS un gène connu par son locus tag (`katG` = 3 papiers sous
    `Rv1908c`, mais 1857 sous `katG`).
 2. **Encart « biologie documentée AILLEURS que chez Mtb »** quand `n_other >= n_tb` (et >= 3 papiers) : la
    génétique mycobactérienne se fait massivement chez *M. smegmatis*, et une part du protéome de *M. tuberculosis*
    n'est comprise que **par procuration**, via ses orthologues. L'atlas doit le DIRE au lecteur — c'est l'inverse
    exact de l'angle mort qu'on vient de corriger (ne regarder que la littérature Mtb).

GARDE-FOU ÉCRIT SUR LA FICHE (impératif) : **« mieux étudié ailleurs » ≠ « fonction établie chez Mtb ».**
Un résultat obtenu chez *M. smegmatis* (espèce NON pathogène, croissance rapide, milieu et régulation différents)
ne se transpose PAS automatiquement à *M. tuberculosis*. Ces travaux sont un CONTEXTE à vérifier, pas un acquis.
Cette couche CITE et CONTEXTUALISE ; elle ne modifie AUCUN verdict ni aucune fonction.

Écrit EN PLACE. Run: python analyses/phase77_literature_known.py   (lancer NU, jamais dans un pipe : cf. CLAUDE.md)
"""
from __future__ import annotations
import csv, json, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
TSV = ROOT / "résultats" / "phase76_literature" / "proteome_literature.tsv"
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL, EMAIL = "annotation_mtbc", "guyeux@gmail.com"
MAX_LISTED = 5
ELSEWHERE_MIN = 3          # seuil de papiers non-TB pour lever l'encart « étudié ailleurs »
SOURCE = ("PubMed (whole): H37Rv locus tag + GENE NAME + ortholog identifiers (Mb…, MMAR_…, MSMEG_…, ML…, MAB_…), "
          "under a mycobacterial context filter; hits verified against the abstract text. Species-context counts "
          "distinguish M. tuberculosis literature from literature on other mycobacteria. phase76/phase77, 2026-07-13")

CAVEAT_ELSEWHERE = (
    "IMPORTANT — 'better studied elsewhere' does NOT mean 'function established in M. tuberculosis'. "
    "Findings obtained in M. smegmatis (a non-pathogenic, fast-growing species with a different lifestyle and "
    "regulation), or in M. marinum / M. leprae / M. abscessus, do NOT transfer automatically to M. tuberculosis. "
    "Treat this body of work as CONTEXT to verify, not as settled knowledge.")
CAVEAT_PLAIN = ("This layer CITES the literature and adds context; it does not change the verdict or the function "
                "stated elsewhere in this fiche.")


def efetch(pmids: list[str]) -> dict[str, dict]:
    out = {}
    for i in range(0, len(pmids), 180):
        lot = pmids[i:i + 180]
        q = urllib.parse.urlencode({"db": "pubmed", "retmode": "xml", "id": ",".join(lot),
                                    "tool": TOOL, "email": EMAIL})
        try:
            with urllib.request.urlopen(f"{EU}efetch.fcgi?{q}", timeout=60) as r:
                root = ET.fromstring(r.read())
        except Exception:
            time.sleep(2)
            continue
        for art in root.iter("PubmedArticle"):
            pmid = art.findtext(".//PMID") or ""
            # itertext() obligatoire : `.text` s'arrête au premier balisage (bug qui jetait 92 % du texte)
            title = "".join(next(art.iter("ArticleTitle"), ET.Element("x")).itertext()).strip()
            year = art.findtext(".//PubDate/Year") or art.findtext(".//PubDate/MedlineDate") or ""
            doi = next((i.text for i in art.iter("ArticleId") if i.get("IdType") == "doi"), None)
            out[pmid] = {"title": title, "doi": doi,
                         "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", "date": year}
        time.sleep(0.4)
        print(f"  efetch {min(i+180, len(pmids))}/{len(pmids)}")
    return out


def main() -> None:
    if not TSV.exists():
        raise SystemExit(f"{TSV} introuvable — lancer phase76_literature_proteome.py d'abord.")
    rows = [r for r in csv.DictReader(TSV.open(), delimiter="\t") if r["n_lit"] not in ("NA", "")]

    # PMIDs à récupérer : les papiers généraux + les papiers « ailleurs » (prioritaires pour l'encart)
    pmids = []
    for r in rows:
        pmids += [p for p in r["pmids"].split(",") if p][:MAX_LISTED]
        pmids += [p for p in r["other_pmids"].split(",") if p][:3]
    meta = efetch(sorted(set(pmids)))
    print(f"\n{len(rows)} gènes connus ; {len(meta)} papiers récupérés")

    n_cited = n_none = 0
    elsewhere = []
    for r in rows:
        rv = r["rv"]
        f = GENES / f"{rv}.json"
        if not f.exists():
            continue
        d = json.loads(f.read_text())
        n_lit = int(r["n_lit"])
        n_tb = int(r["n_tb"] or 0)
        n_other = int(r["n_other"] or 0)
        n_ortho = int(r["n_ortho"] or 0)
        by_sp = dict(x.split(":") for x in (r["other_species"] or "").split(";") if ":" in x)
        papers = [meta[p] for p in [x for x in r["pmids"].split(",") if x][:MAX_LISTED] if p in meta]
        opapers = [meta[p] for p in [x for x in r["other_pmids"].split(",") if x][:3] if p in meta]

        is_elsewhere = n_other >= ELSEWHERE_MIN and n_other >= n_tb
        if n_lit == 0:
            layer = {
                "n_papers": 0, "category": "unstudied",
                "category_label": "no publication mentions this gene",
                "summary": ("No publication mentions this gene in its title or abstract — not under its H37Rv locus "
                            "tag, not under its gene name, and not under any ortholog identifier. Its annotation "
                            "rests on sequence/structure evidence, with no primary study behind it."),
                "papers": [], "caveat": CAVEAT_PLAIN, "source": SOURCE,
            }
            n_none += 1
        else:
            sp_txt = ", ".join(f"{k} ({v})" for k, v in sorted(by_sp.items(), key=lambda x: -int(x[1])))
            if is_elsewhere:
                summary = (f"{n_lit} publication(s) discuss this gene. **Its biology is documented at least as much "
                           f"OUTSIDE M. tuberculosis as within it** ({n_other} papers in a non-TB mycobacterial "
                           f"context — {sp_txt} — vs {n_tb} in a TB context). Mycobacterial genetics is largely done "
                           f"in M. smegmatis, so part of what is 'known' about this gene is known by proxy.")
                elsewhere.append((n_other, n_tb, rv, d.get("gene") or "-", sp_txt))
            else:
                summary = (f"{n_lit} publication(s) discuss this gene "
                           f"({n_tb} in a M. tuberculosis context"
                           + (f", {n_other} in other mycobacteria — {sp_txt}" if n_other else "") + ").")
            layer = {
                "n_papers": n_lit, "n_papers_listed": len(papers),
                "category": "studied_elsewhere" if is_elsewhere else "cited",
                "category_label": ("biology documented at least as much outside M. tuberculosis"
                                   if is_elsewhere else "cited in the literature"),
                "n_tb": n_tb, "n_other": n_other, "n_ortholog_tag": n_ortho,
                "other_species": by_sp,
                "summary": summary,
                "papers": papers,
                "papers_elsewhere": opapers if is_elsewhere else [],
                "caveat": (CAVEAT_ELSEWHERE if is_elsewhere else CAVEAT_PLAIN),
                "source": SOURCE,
            }
            n_cited += 1
        d["literature"] = layer
        f.write_text(json.dumps(d, indent=2, ensure_ascii=False))

    print(f"\ncouche `literature` écrite : {n_cited} cités + {n_none} sans mention = {n_cited+n_none} gènes connus.")
    print(f"  -> le protéome ENTIER (3906) est désormais ancré dans la littérature.")
    print(f"\n>>> {len(elsewhere)} gènes dont la BIOLOGIE EST DOCUMENTÉE AU MOINS AUTANT HORS Mtb. Top 30 :")
    for no, tb, rv, name, sp in sorted(elsewhere, reverse=True)[:30]:
        print(f"    {rv:9} {name:9} non-TB={no:3} vs TB={tb:4}   [{sp}]")
    print("\nGARDE-FOU rappelé sur chaque fiche concernée : « mieux étudié ailleurs » != « fonction établie chez Mtb ».")


if __name__ == "__main__":
    main()
