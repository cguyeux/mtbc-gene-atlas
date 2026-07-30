#!/usr/bin/env python3
"""phase73_pubmed_sweep.py -- P16.2a : balayage littérature PubMed, MULTI-ALIAS et VÉRIFIÉ SUR TEXTE.

=== POURQUOI CE SCRIPT A ÉTÉ RÉÉCRIT TROIS FOIS (chaque version a produit un résultat FAUX) ===

v1 — cherchait le SEUL locus tag `Rv####`.
   FAUX NÉGATIFS : un gène « hypothétique » chez H37Rv peut être caractérisé chez une espèce voisine sous un AUTRE
   identifiant (la génétique mycobactérienne se fait massivement chez *M. smegmatis*). -> ajout des orthologues.

v2 — interrogeait `ALIAS[tiab]` pour chaque alias.
   FAUX POSITIFS MASSIFS : **PubMed DÉCOUPE les identifiants à tiret bas.** `MAB_2018[tiab]` renvoyait **230 papiers**
   … sur les anticorps monoclonaux anti-CGRP : PubMed le lit comme `MAB` (= mAb, anticorps monoclonal, ultra-fréquent)
   ET `2018` (= une année). Idem `MAB_1996` -> 68, `MAB_1983` -> 23. Les formats `MSMEG_`/`MMAR_` semblaient sains
   uniquement parce que leur préfixe est un token RARE (pas de co-occurrence) : le découpage a lieu quand même.
   LEÇON : valider UNE famille d'alias ne valide PAS les autres.

v2bis — les guillemets ? `"MAB_2018"[tiab]` -> 0 (les faux positifs tombent)… MAIS `"MSMEG_0232"[tiab]` -> 0 alors que
   la version nue en trouve 1 (papier réel !). **Les guillemets détruisent aussi de VRAIS positifs.** Ni la requête
   nue (faux +) ni la requête entre guillemets (faux −) n'est fiable.

v3 — vérification sur texte… mais avec un BUG ElementTree : `findtext()`/`.text` ne renvoient que le texte AVANT le
   premier élément enfant. Les résumés PubMed contiennent du balisage (`<i>Mycobacterium</i>`, `<sup>`), donc on
   n'extrayait que **162 caractères sur 1989** (≈ 92 % du texte jeté) → la vérification REJETAIT de VRAIS positifs
   (Rv2082 est bien dans son résumé : « …fadE18, Rv0988, and Rv2082 variants… »). Détecté UNIQUEMENT parce que le
   résultat (25 gènes) était INFÉRIEUR à celui de tbmonitor (32) alors qu'il aurait dû être supérieur : **une
   incohérence avec un contre-jeu de référence est le seul garde-fou qui attrape une vérification silencieusement
   cassée.** -> `itertext()`.

v4 (CE SCRIPT) — la seule méthode saine : **requête SENSIBLE + filtre de contexte, puis VÉRIFICATION SUR LE TEXTE**
   (texte extrait avec `itertext()`, jamais `.text`).
   1. `esearch` sur l'UNION des alias, ET un filtre de contexte mycobactérien (écrase le bruit : MAB_2018 230 -> 3,
      MAB_1996 68 -> 0, tout en PRÉSERVANT les vrais : MSMEG_0232/0241/5827 intacts).
   2. `efetch` des candidats, puis **regex `\\bALIAS\\b` sur titre+résumé** : un papier n'est compté QUE si
      l'identifiant exact y figure vraiment. La regex donne AUSSI l'attribution par alias (gratuite).

ALIAS : Rv (H37Rv) + orthologues Mb (M. bovis), MMAR_ (M. marinum), MSMEG_ (M. smegmatis), ML (M. leprae),
MAB_ (M. abscessus), depuis la couche `orthologs` (RBH DIAMOND). *M. orygis* EXCLU (locus `RJtmp_` = tag interne
jamais publié).

Sortie : résultats/phase73_pubmed/pubmed_verified_<set>.tsv
Run: python analyses/phase73_pubmed_sweep.py --set dark|hypo   (lancer NU, jamais dans un pipe : cf. CLAUDE.md)
"""
from __future__ import annotations
import glob, json, re, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

try:
    sys.stdout.reconfigure(line_buffering=True)   # type: ignore[attr-defined]  (stub TextIO incomplet ; vérifié OK à l'exécution)
except Exception:
    pass

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "résultats" / "phase73_pubmed"
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL, EMAIL = "annotation_mtbc", "guyeux@gmail.com"
SLEEP = 0.4
MAX_VERIFY = 60          # candidats vérifiés par gène (au-delà : signalé, pas silencieusement tronqué)
SKIP_SPECIES = {"M. orygis"}
CTX = ("(mycobacteri*[tiab] OR mycobacterium[tiab] OR tuberculosis[tiab] OR abscessus[tiab] "
       "OR smegmatis[tiab] OR marinum[tiab] OR leprae[tiab] OR bovis[tiab])")
CONTROL_EVERY, CONTROL_TERM, CONTROL_MIN = 120, "katG", 500


def _get(url: str, retries: int = 4):
    for a in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=40) as r:
                return r.read()
        except Exception:
            time.sleep(1.5 * (a + 1))
    return None


def esearch(term: str, retmax: int = 200):
    q = urllib.parse.urlencode({"db": "pubmed", "retmode": "json", "retmax": retmax,
                                "term": term, "tool": TOOL, "email": EMAIL})
    b = _get(f"{EU}esearch.fcgi?{q}")
    if b is None:
        return None
    d = json.loads(b)["esearchresult"]
    return int(d["count"]), d.get("idlist", [])


def efetch_texts(pmids: list[str]) -> dict[str, str]:
    """PMID -> titre + résumé (texte brut) pour vérification regex."""
    if not pmids:
        return {}
    q = urllib.parse.urlencode({"db": "pubmed", "retmode": "xml", "id": ",".join(pmids),
                                "tool": TOOL, "email": EMAIL})
    b = _get(f"{EU}efetch.fcgi?{q}")
    if b is None:
        return {}
    out = {}
    try:
        root = ET.fromstring(b)
    except ET.ParseError:
        return {}
    for art in root.iter("PubmedArticle"):
        pmid = art.findtext(".//PMID") or ""
        # PIÈGE ElementTree : `.text` / `findtext()` ne renvoient QUE le texte AVANT le premier enfant.
        # Un résumé PubMed contient du balisage (<i>Mycobacterium</i>, <sup>, …) -> on perdait ~92 % du texte
        # (162 car. au lieu de 1989 !) et la vérification rejetait des VRAIS positifs. `itertext()` prend tout.
        title = "".join(next(art.iter("ArticleTitle"), ET.Element("x")).itertext())
        abst = " ".join("".join(t.itertext()) for t in art.iter("AbstractText"))
        out[pmid] = f"{title} {abst}"
    return out


def aliases(d: dict) -> dict[str, str]:
    al = {"H37Rv": d["rv"]}
    for sp, v in (((d.get("orthologs") or {}).get("orthologs")) or {}).items():
        if sp in SKIP_SPECIES:
            continue
        loc = (v or {}).get("locus")
        if loc and not loc.startswith("RJtmp"):
            al[sp] = loc
    return al


def targets(which: str) -> list[dict]:
    out = []
    for f in glob.glob(str(GENES / "*.json")):
        d = json.loads(Path(f).read_text())
        is_dark = d.get("verdict") == "dark"
        hypo = "hypothetical" in (d.get("product_h37rv") or "").lower()
        if which == "dark" and not is_dark:
            continue
        if which == "hypo" and (is_dark or not hypo):
            continue
        out.append(d)
    return sorted(out, key=lambda d: d["rv"])


def main() -> None:
    which = sys.argv[sys.argv.index("--set") + 1] if "--set" in sys.argv else "dark"
    genes = targets(which)
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"balayage PubMed VÉRIFIÉ ({which}) : {len(genes)} gènes\n")

    rows, truncated = [], []
    for i, d in enumerate(genes, 1):
        if i % CONTROL_EVERY == 1:
            c = esearch(f"{CONTROL_TERM}[tiab]", retmax=0)
            if not c or c[0] < CONTROL_MIN:
                raise SystemExit(f"ABORT: témoin '{CONTROL_TERM}' = {c} -> API cassée, on NE conclut PAS sur des zéros.")
            print(f"  [témoin OK: {CONTROL_TERM}={c[0]}]")
            time.sleep(SLEEP)

        al = aliases(d)
        union = " OR ".join(f"{v}[tiab]" for v in al.values())
        r = esearch(f"({union}) AND {CTX}")
        time.sleep(SLEEP)
        if r is None:
            rows.append({"rv": d["rv"], "verdict": d.get("verdict"), "n_cand": "NA", "n_verified": "NA",
                         "by_alias": "", "pmids": "", "aliases": ";".join(al.values())})
            print(f"  {d['rv']}: ÉCHEC esearch -> exclu (jamais compté comme 0)")
            continue
        n_cand, pmids = r

        by, confirmed = {}, []
        if pmids:
            if n_cand > MAX_VERIFY:
                truncated.append((d["rv"], n_cand))
            texts = efetch_texts(pmids[:MAX_VERIFY])
            time.sleep(SLEEP)
            pats = {lab: re.compile(rf"(?<![A-Za-z0-9_]){re.escape(tag)}(?![A-Za-z0-9_])", re.I)
                    for lab, tag in al.items()}
            for pmid, txt in texts.items():
                hit = [lab for lab, p in pats.items() if p.search(txt)]
                if hit:
                    confirmed.append(pmid)
                    for lab in hit:
                        by[lab] = by.get(lab, 0) + 1

        rows.append({"rv": d["rv"], "verdict": d.get("verdict"), "n_cand": n_cand,
                     "n_verified": len(confirmed),
                     "by_alias": ";".join(f"{k}:{v}" for k, v in sorted(by.items())),
                     "pmids": ",".join(confirmed[:12]), "aliases": ";".join(al.values())})
        if confirmed and "H37Rv" not in by:
            print(f"  *** {d['rv']} : {len(confirmed)} papier(s) VÉRIFIÉ(S) via ORTHOLOGUE seul -> {by}")
        if i % 50 == 0:
            print(f"  [{i}/{len(genes)}]")

    p = OUT / f"pubmed_verified_{which}.tsv"
    with p.open("w") as fh:
        fh.write("rv\tverdict\tn_candidates\tn_verified\tby_alias\tpmids\taliases\n")
        for r in rows:
            fh.write(f"{r['rv']}\t{r['verdict']}\t{r['n_cand']}\t{r['n_verified']}\t{r['by_alias']}\t{r['pmids']}\t{r['aliases']}\n")

    ok = [r for r in rows if r["n_verified"] != "NA"]
    lit = [r for r in ok if r["n_verified"] > 0]
    ortho = [r for r in lit if "H37Rv" not in (r["by_alias"] or "")]
    print(f"\n===== {which} : {len(ok)} interrogés =====")
    print(f"  AVEC littérature vérifiée : {len(lit)}")
    print(f"  SANS aucune mention       : {len(ok)-len(lit)}")
    print(f"  dont littérature UNIQUEMENT via ORTHOLOGUE (invisible au tag Rv) : {len(ortho)}")
    for r in sorted(ortho, key=lambda r: -r["n_verified"])[:25]:
        print(f"     {r['rv']:9} {r['n_verified']:3} papier(s)  [{r['by_alias']}]")
    if truncated:
        print(f"\n  ⚠ {len(truncated)} gène(s) avec plus de {MAX_VERIFY} candidats (vérif plafonnée, pas silencieuse) :")
        for rv, n in truncated[:10]:
            print(f"     {rv}: {n} candidats")
    print(f"\nTable : {p}")


if __name__ == "__main__":
    main()
