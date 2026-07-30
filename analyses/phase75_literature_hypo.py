#!/usr/bin/env python3
"""phase75_literature_hypo.py -- P16.2a : couche `literature` pour les 892 hypothétiques NON-dark.

Suite de phase70 (qui couvrait les 219 dark). Source : le balayage PubMed VÉRIFIÉ multi-alias de `phase73`
(`résultats/phase73_pubmed/pubmed_verified_hypo.tsv`) — tag H37Rv + tous les identifiants d'orthologues
(Mb…, MMAR_…, MSMEG_…, ML…, MAB_…), chaque hit confirmé sur le TEXTE du résumé.

CE QUE CE BALAYAGE A RÉVÉLÉ (le trou est BIEN plus large que le cas Rv2660c/H56) :
  - 384 des 892 ont de la littérature vérifiée ; **342 d'entre eux n'ont AUCUNE référence dans leur fiche**.
  - 14 angles morts sévères (>=8 papiers, 0 citation), dont **Rv0678 (57 papiers) = MmpR5, LE régulateur de la
    résistance à la bédaquiline**, **Rv1860 (15) = Apa/MPT32, antigène immunodominant**, **Rv1813c (10) =
    composant d'un vaccin multistage**.
  - **25 gènes ne sont cités QUE sous le nom de leur orthologue** (surtout MSMEG_ : la génétique mycobactérienne
    se fait chez *M. smegmatis*), donc totalement invisibles à une recherche par tag Rv.

HONNÊTETÉ (cadrage) : l'atlas donne le plus souvent la BONNE FONCTION (Rv0678 est bien annoté « MmpR »). Ce qui
manquait n'est pas la fonction mais l'**ancrage littérature** : les citations, et le CONTEXTE (« ce gène est le
déterminant de résistance à la bédaquiline »). Cette couche ne change AUCUN verdict : elle cite.

Écrit EN PLACE. Run: python analyses/phase75_literature_hypo.py   (lancer NU, jamais dans un pipe : cf. CLAUDE.md)
"""
from __future__ import annotations
import csv, json, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
TSV = ROOT / "résultats" / "phase73_pubmed" / "pubmed_verified_hypo.tsv"
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL, EMAIL = "annotation_mtbc", "guyeux@gmail.com"
MAX_LISTED = 5          # papiers listés sur la fiche (les plus récents)
SOURCE = ("PubMed (whole), MULTI-ALIAS sweep: H37Rv locus tag AND every ortholog identifier (M. bovis Mb…, "
          "M. marinum MMAR_…, M. smegmatis MSMEG_…, M. leprae ML…, M. abscessus MAB_…), each hit VERIFIED "
          "against the abstract text (word-boundary regex). phase73/phase75, 2026-07-13")


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
            # itertext() obligatoire : `.text` s'arrête au premier balisage (cf. bug corrigé en phase73)
            title = "".join(next(art.iter("ArticleTitle"), ET.Element("x")).itertext())
            year = art.findtext(".//PubDate/Year") or art.findtext(".//PubDate/MedlineDate") or ""
            doi = next((i.text for i in art.iter("ArticleId") if i.get("IdType") == "doi"), None)
            out[pmid] = {"title": title.strip(), "doi": doi,
                         "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", "date": year}
        time.sleep(0.4)
        print(f"  efetch {min(i+180, len(pmids))}/{len(pmids)}")
    return out


def main() -> None:
    rows = [r for r in csv.DictReader(TSV.open(), delimiter="\t") if r["n_verified"] not in ("NA", "0")]
    all_pmids = []
    for r in rows:
        all_pmids += [p for p in r["pmids"].split(",") if p][:MAX_LISTED]
    print(f"{len(rows)} gènes avec littérature ; {len(all_pmids)} PMIDs à récupérer")
    meta = efetch(sorted(set(all_pmids)))

    n = 0
    blind = []
    for r in rows:
        f = GENES / f"{r['rv']}.json"
        d = json.loads(f.read_text())
        pm = [p for p in r["pmids"].split(",") if p][:MAX_LISTED]
        papers = [meta[p] for p in pm if p in meta]
        by = dict(x.split(":") for x in r["by_alias"].split(";") if ":" in x)
        ortho_only = "H37Rv" not in by
        nv = int(r["n_verified"])
        n_refs = len(d.get("references") or [])
        if nv >= 8 and n_refs == 0:
            blind.append((nv, r["rv"], (d.get("function_revised") or "")[:45]))

        if ortho_only:
            sp = ", ".join(by)
            summary = (f"Invisible under its H37Rv locus tag: this gene appears in the literature ONLY under its "
                       f"ortholog identifier(s) ({sp}). Searching '{r['rv']}' alone finds nothing.")
        else:
            summary = (f"{nv} publication(s) mention this gene (title/abstract), verified against the text. "
                       "The atlas states the function; these are the primary sources that discuss it.")
        d["literature"] = {
            "n_papers": nv,
            "n_papers_listed": len(papers),
            "category": "cited_ortholog_only" if ortho_only else "cited",
            "category_label": ("cited in the literature ONLY under an ortholog name"
                               if ortho_only else "cited in the literature"),
            "by_alias": by,
            "summary": summary,
            "papers": papers,
            "caveat": ("This layer CITES the literature, it does not change the verdict or the function. "
                       "A gene may be heavily studied (as an antigen, a resistance determinant, a drug target) "
                       "while its molecular function is settled elsewhere in the fiche."),
            "source": SOURCE,
        }
        f.write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1

    # les 508 sans mention : on le DIT explicitement (une absence vérifiée est une information)
    seen = {r["rv"] for r in csv.DictReader(TSV.open(), delimiter="\t")}
    cited = {r["rv"] for r in rows}
    m = 0
    for rv in sorted(seen - cited):
        f = GENES / f"{rv}.json"
        d = json.loads(f.read_text())
        d["literature"] = {
            "n_papers": 0, "n_papers_listed": 0, "category": "unstudied",
            "category_label": "no publication mentions this gene",
            "summary": ("No publication in PubMed mentions this gene in its title or abstract — not under its H37Rv "
                        "locus tag, nor under any of its ortholog identifiers (M. bovis, M. marinum, M. smegmatis, "
                        "M. leprae, M. abscessus). The atlas annotation rests on sequence/structure evidence, not on "
                        "a primary study of this gene."),
            "papers": [], "by_alias": {},
            "caveat": "A verified absence of literature is itself information: it flags an annotation with no primary study behind it.",
            "source": SOURCE,
        }
        f.write_text(json.dumps(d, indent=2, ensure_ascii=False))
        m += 1

    print(f"\ncouche `literature` écrite : {n} gènes cités + {m} gènes sans mention = {n+m} fiches.")
    print(f"\n>>> ANGLES MORTS corrigés (>=8 papiers, 0 référence dans la fiche) : {len(blind)}")
    for nv, rv, fn in sorted(blind, reverse=True):
        print(f"    {rv:9} {nv:3} papiers  {fn}")


if __name__ == "__main__":
    main()
