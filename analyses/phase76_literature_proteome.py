#!/usr/bin/env python3
"""phase76_literature_proteome.py -- P16.2a-quater : littérature du RESTE du protéome + « mieux étudié ailleurs ».

DEUX QUESTIONS (demande CG 2026-07-13) :
 A. Étendre la couche `literature` aux ~2795 gènes NON-hypothétiques (les « gènes connus »), pour que TOUT le
    protéome soit ancré dans la littérature (phase70+phase75 couvrent déjà les 1111 hypothétiques).
 B. **L'angle ORTHOLOGUE vaut aussi pour les gènes CONNUS** : un gène nommé chez Mtb peut avoir été étudié PLUS
    PROFONDÉMENT chez une autre espèce (la génétique mycobactérienne se fait chez *M. smegmatis*). On veut donc
    faire remonter les gènes dont le corpus ORTHOLOGUE est substantiel.

PIÈGE CAPITAL POUR LES GÈNES CONNUS : **la littérature ne les appelle JAMAIS par leur locus tag.** Mesuré :
`katG` = **3** papiers sous `Rv1908c[tiab]` mais **1857** sous `katG[tiab]` ; `rpoB` = 2 vs 5556. Il FAUT donc
ajouter le NOM DE GÈNE aux alias. Mais les noms courts sont explosifs (`eis` seul = 9723 papiers, c'est aussi un
acronyme courant ; `ald` = 10893) → le **filtre de contexte mycobactérien** les ramène à 207 et 40. C'est lui, plus
la vérification sur texte, qui rend la mesure honnête.

PIÈGE N°2 (découvert au test) : compter les seuls LOCUS TAGS d'orthologues SOUS-ESTIME « étudié ailleurs », car
chez *M. smegmatis* le gène `katG` s'appelle AUSSI `katG` — les papiers le nomment, ils n'écrivent pas `MSMEG_6384`.
Pour répondre honnêtement à la question B, il faut comparer les **CONTEXTES D'ESPÈCE**, pas les identifiants.

TROIS MÉTRIQUES, à ne pas confondre :
 - `n_lit`   : volume total (union Rv + nom de gène + tags d'orthologues) sous filtre de contexte mycobactérien.
               Vérifié sur texte quand le compte est petit (<= VERIFY_UNDER) ; au-delà le volume parle de lui-même.
 - `n_tb`    : le même jeu d'alias, mais en **contexte *M. tuberculosis*** (tuberculosis / H37Rv / bovis).
 - `n_other` : le même jeu d'alias en **contexte NON-TB** (*smegmatis* / *marinum* / *leprae* / *abscessus*),
               **vérifié sur le texte**. C'est LA métrique de la question B : un gène dont `n_other` rivalise avec
               (ou dépasse) `n_tb` est un gène **dont la biologie est explorée ailleurs que chez Mtb** — typiquement
               chez *M. smegmatis*, cheval de trait de la génétique mycobactérienne.
 - `n_ortho` : sous-ensemble de `n_other` où c'est le **LOCUS TAG** de l'orthologue qui apparaît (MSMEG_…, MMAR_…) :
               ces papiers sont totalement invisibles à une recherche par tag `Rv`.

Débit : NCBI autorise 3 req/s sans clé. On se limite à 2,5 req/s (jeton global partagé) avec 3 workers.

Sortie : résultats/phase76_literature/proteome_literature.tsv
Run: python analyses/phase76_literature_proteome.py   (lancer NU, jamais dans un pipe : cf. CLAUDE.md)
"""
from __future__ import annotations
import csv, glob, json, re, sys, threading, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

try:
    sys.stdout.reconfigure(line_buffering=True)   # type: ignore[attr-defined]  (stub TextIO incomplet ; vérifié OK à l'exécution)
except Exception:
    pass

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "résultats" / "phase76_literature"
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL, EMAIL = "annotation_mtbc", "guyeux@gmail.com"
CTX = ("(mycobacteri*[tiab] OR mycobacterium[tiab] OR tuberculosis[tiab] OR abscessus[tiab] "
       "OR smegmatis[tiab] OR marinum[tiab] OR leprae[tiab] OR bovis[tiab])")
CTX_TB = "(tuberculosis[tiab] OR H37Rv[tiab] OR bovis[tiab] OR BCG[tiab])"
CTX_OTHER = "(smegmatis[tiab] OR marinum[tiab] OR leprae[tiab] OR abscessus[tiab])"
SKIP_SPECIES = {"M. orygis"}          # locus RJtmp_ = tag interne jamais publié
VERIFY_UNDER = 25                     # au-delà, on ne vérifie pas le corpus principal (volume auto-parlant)
MAX_ORTHO_VERIFY = 40
WORKERS, RATE = 3, 2.5                # req/s global (limite NCBI sans clé : 3/s)

_lock = threading.Lock()
_last = [0.0]


def _throttle():
    with _lock:
        dt = time.monotonic() - _last[0]
        wait = (1.0 / RATE) - dt
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.monotonic()


def _get(url: str, retries: int = 4):
    for a in range(retries):
        _throttle()
        try:
            with urllib.request.urlopen(url, timeout=40) as r:
                return r.read()
        except Exception:
            time.sleep(1.0 * (a + 1))
    return None


def esearch(term: str, retmax: int = 40, retries: int = 3):
    """PIÈGE : NCBI peut renvoyer un HTTP **200** dont le corps est {"esearchresult": {"ERROR": …}} — sans `count`.
    Ce n'est PAS une exception réseau, donc les retries de `_get` ne le voient pas : le `d["count"]` lève alors
    un KeyError qui, dans un ThreadPoolExecutor, TUE TOUT LE RUN. (Vécu : 2 h de calcul perdues.)
    -> on retente sur corps d'erreur, puis on renvoie None (le gène sera marqué NA, jamais compté comme 0)."""
    q = urllib.parse.urlencode({"db": "pubmed", "retmode": "json", "retmax": retmax,
                                "term": term, "tool": TOOL, "email": EMAIL})
    for a in range(retries):
        b = _get(f"{EU}esearch.fcgi?{q}")
        if b is None:
            continue
        try:
            d = json.loads(b)["esearchresult"]
        except Exception:
            time.sleep(1.0 * (a + 1))
            continue
        if "count" in d:
            return int(d["count"]), d.get("idlist", [])
        time.sleep(2.0 * (a + 1))          # corps d'erreur NCBI (transitoire) -> on retente
    return None


def efetch_texts(pmids: list[str]) -> dict[str, str]:
    if not pmids:
        return {}
    q = urllib.parse.urlencode({"db": "pubmed", "retmode": "xml", "id": ",".join(pmids),
                                "tool": TOOL, "email": EMAIL})
    b = _get(f"{EU}efetch.fcgi?{q}")
    if b is None:
        return {}
    try:
        root = ET.fromstring(b)
    except ET.ParseError:
        return {}
    out = {}
    for art in root.iter("PubmedArticle"):
        pmid = art.findtext(".//PMID") or ""
        # itertext() OBLIGATOIRE : `.text` s'arrête au premier balisage (bug qui jetait 92 % du résumé)
        title = "".join(next(art.iter("ArticleTitle"), ET.Element("x")).itertext())
        abst = " ".join("".join(t.itertext()) for t in art.iter("AbstractText"))
        out[pmid] = f"{title} {abst}"
    return out


def pat(tag: str):
    return re.compile(rf"(?<![A-Za-z0-9_]){re.escape(tag)}(?![A-Za-z0-9_])", re.I)


def analyse(d: dict) -> dict:
    rv = d["rv"]
    name = (d.get("gene") or "").strip()
    ortho = {}
    for sp, v in (((d.get("orthologs") or {}).get("orthologs")) or {}).items():
        if sp in SKIP_SPECIES:
            continue
        loc = (v or {}).get("locus")
        if loc and not loc.startswith("RJtmp"):
            ortho[sp] = loc

    main = [rv] + ([name] if name and len(name) >= 3 else [])
    all_al = main + list(ortho.values())
    # un nom comme `erm(37)` (Rv1988) contient des parenthèses qui font MAL-PARSER la requête
    # (mesuré : `erm(37)[tiab]` -> 8 hits, `"erm(37)"[tiab]` -> 14) -> mettre entre guillemets si non alphanumérique.
    def q(a: str) -> str:
        return f'"{a}"[tiab]' if re.search(r"[^A-Za-z0-9_.-]", a) else f"{a}[tiab]"
    ALIAS = "(" + " OR ".join(q(a) for a in all_al) + ")"
    NA = {"rv": rv, "name": name, "n_lit": "NA", "n_tb": "NA", "n_other": "NA", "n_ortho": "NA",
          "other_species": "", "ortho_loci": ";".join(ortho.values()), "pmids": "", "other_pmids": ""}

    # --- 1. volume total (contexte mycobactérien) ---
    u = esearch(f"{ALIAS} AND {CTX}", retmax=40)
    if u is None:
        return NA
    n_lit, pmids = u
    verified_main = None
    if 0 < n_lit <= VERIFY_UNDER:          # petit compte -> vérifier (tue les faux positifs de nom court)
        texts = efetch_texts(pmids[:VERIFY_UNDER])
        pats = [pat(a) for a in all_al]
        verified_main = sum(1 for t in texts.values() if any(p.search(t) for p in pats))

    # --- 2. contexte M. tuberculosis ---
    t = esearch(f"{ALIAS} AND {CTX_TB}", retmax=0)
    n_tb = t[0] if t else 0

    # --- 3. contexte NON-TB (smegmatis/marinum/leprae/abscessus) ---
    # n_other = VOLUME BRUT (non plafonné, comparable à n_tb) ; la vérification porte sur un ÉCHANTILLON
    # de MAX_ORTHO_VERIFY papiers et donne un taux de précision + l'attribution par espèce.
    n_other, n_other_checked, n_other_ok, by_sp, n_ortho, opm = 0, 0, 0, {}, 0, []
    o = esearch(f"{ALIAS} AND {CTX_OTHER}", retmax=MAX_ORTHO_VERIFY)
    if o and o[0]:
        n_other = o[0]                      # volume brut, PAS plafonné
        texts = efetch_texts(o[1][:MAX_ORTHO_VERIFY])
        n_other_checked = len(texts)
        alias_pats = [pat(a) for a in all_al]
        loci_pats = {sp: pat(loc) for sp, loc in ortho.items()}
        sp_pats = {"M. smegmatis": pat("smegmatis"), "M. marinum": pat("marinum"),
                   "M. leprae": pat("leprae"), "M. abscessus": pat("abscessus")}
        for pmid, txt in texts.items():
            if not any(p.search(txt) for p in alias_pats):
                continue                                    # le gène n'est pas vraiment dans le texte
            n_other_ok += 1
            opm.append(pmid)
            for sp, p in sp_pats.items():
                if p.search(txt):
                    by_sp[sp] = by_sp.get(sp, 0) + 1
            if any(p.search(txt) for p in loci_pats.values()):
                n_ortho += 1                                # cité par le LOCUS TAG de l'orthologue
    prec = round(n_other_ok / n_other_checked, 2) if n_other_checked else None
    return {"rv": rv, "name": name,
            "n_lit": verified_main if verified_main is not None else n_lit,
            "n_lit_is_verified": verified_main is not None,
            "n_tb": n_tb,
            "n_other": n_other,                 # volume brut non-TB (comparable à n_tb)
            "n_other_checked": n_other_checked, # taille de l'échantillon vérifié
            "precision": prec,                  # fraction de l'échantillon où le gène est vraiment dans le texte
            "n_ortho": n_ortho,
            "other_species": ";".join(f"{k}:{v}" for k, v in sorted(by_sp.items())),
            "ortho_loci": ";".join(ortho.values()),
            "pmids": ",".join(pmids[:6]), "other_pmids": ",".join(opm[:6])}


def main() -> None:
    genes = []
    for f in glob.glob(str(GENES / "*.json")):
        d = json.loads(Path(f).read_text())
        hypo = "hypothetical" in (d.get("product_h37rv") or "").lower()
        if hypo or d.get("verdict") == "dark":
            continue                                   # déjà couverts par phase70 / phase75
        genes.append(d)
    genes.sort(key=lambda d: d["rv"])
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"reste du protéome (gènes connus) : {len(genes)} gènes")
    print(f"débit {RATE} req/s, {WORKERS} workers -> estimation ~{len(genes)*3/RATE/60:.0f} min\n")

    c = esearch("katG[tiab]", retmax=0)                # témoin positif d'amorçage
    if not c or c[0] < 500:
        raise SystemExit(f"ABORT: témoin katG = {c} -> API cassée, on ne conclut pas.")
    print(f"  [témoin OK: katG={c[0]}]")

    p = OUT / "proteome_literature.tsv"
    cols = ["rv", "name", "n_lit", "n_lit_is_verified", "n_tb", "n_other", "n_other_checked",
            "precision", "n_ortho", "other_species", "ortho_loci", "pmids", "other_pmids"]

    # CHECKPOINT + REPRISE : on écrit CHAQUE ligne dès qu'elle est prête. Un run de 2 h contre une API
    # externe DOIT être reprenable — sinon une seule erreur transitoire anéantit tout le calcul (vécu).
    done_rv = set()
    if p.exists():
        with p.open() as fh:
            done_rv = {l.split("\t")[0] for l in fh.readlines()[1:] if l.strip()}
        genes = [d for d in genes if d["rv"] not in done_rv]
        print(f"  REPRISE : {len(done_rv)} gènes déjà faits, {len(genes)} restants\n")
    fh_out = p.open("a" if done_rv else "w")
    if not done_rv:
        fh_out.write("\t".join(cols) + "\n")
        fh_out.flush()

    rows, done = [], [0]
    def work(d):
        try:
            r = analyse(d)
        except Exception as e:                 # ISOLATION : un gène qui échoue ne tue pas le run
            r = {"rv": d["rv"], "name": d.get("gene") or "", "n_lit": "NA", "n_tb": "NA",
                 "n_other": "NA", "n_ortho": "NA", "other_species": "", "ortho_loci": "",
                 "pmids": "", "other_pmids": ""}
            with _lock:
                print(f"  !! {d['rv']} : ÉCHEC ({type(e).__name__}: {e}) -> NA, jamais compté comme 0")
        with _lock:
            fh_out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
            fh_out.flush()                      # checkpoint immédiat
            done[0] += 1
            if done[0] % 150 == 0:
                print(f"  [{done[0]}/{len(genes)}]")
            no, tb = (int(r["n_other"]) if str(r.get("n_other")).isdigit() else 0,
                      int(r["n_tb"]) if str(r.get("n_tb")).isdigit() else 0)
            if no >= 3 and no >= tb:            # étudié AILLEURS autant ou plus que chez Mtb
                print(f"  *** {r['rv']} ({r['name'] or '-'}) : non-TB={no} vs TB={tb} -> [{r['other_species']}]")
        return r

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        rows = list(ex.map(work, genes))
    fh_out.close()

    # relire le TSV complet (inclut la reprise éventuelle)
    rows = list(csv.DictReader(p.open(), delimiter="\t"))
    ok = [r for r in rows if r["n_lit"] != "NA"]
    lit = [r for r in ok if int(r["n_lit"]) > 0]
    # QUESTION B : gènes dont la biologie est explorée AILLEURS que chez Mtb
    elsewhere = sorted([r for r in ok if int(r["n_other"]) >= 3 and int(r["n_other"]) >= int(r["n_tb"])],
                       key=lambda r: -int(r["n_other"]))
    invisible = sorted([r for r in ok if int(r["n_ortho"]) > 0], key=lambda r: -int(r["n_ortho"]))
    print(f"\n===== {len(ok)} gènes connus interrogés =====")
    print(f"  avec littérature : {len(lit)} | AUCUNE mention : {len(ok)-len(lit)}")
    print(f"\n>>> QUESTION B — {len(elsewhere)} gènes ÉTUDIÉS AUTANT OU PLUS HORS Mtb (non-TB >= TB, >=3 papiers). Top 30 :")
    for r in elsewhere[:30]:
        print(f"    {r['rv']:9} {str(r['name'] or '-'):9} non-TB={r['n_other']:3} vs TB={r['n_tb']:4}  [{r['other_species']}]")
    print(f"\n>>> {len(invisible)} gènes cités par le LOCUS TAG de leur orthologue (invisibles à une recherche Rv). Top 20 :")
    for r in invisible[:20]:
        print(f"    {r['rv']:9} {str(r['name'] or '-'):9} ortho-tag={r['n_ortho']:3}  [{r['other_species']}]")
    print(f"\nTable : {p}")


if __name__ == "__main__":
    main()
