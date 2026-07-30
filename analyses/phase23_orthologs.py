#!/usr/bin/env python3
"""phase23_orthologs.py -- orthologues nommés par espèce (P5.5).

Champ de parité Mycobrowser : par gène H37Rv, le locus tag orthologue chez les autres
mycobactéries de référence (M. bovis, M. leprae, M. marinum, M. smegmatis, M. orygis,
M. abscessus). Mycobrowser affiche ce panel sur la page live mais NE l'exporte PAS dans son
fichier plat (colonnes vides). On le CALCULE nous-mêmes par reciprocal-best-hit DIAMOND
(méthode transparente et reproductible), au lieu de scraper 3906 pages — et c'est plus honnête
que recopier (cf. garde-fou « Mycobrowser à vérifier, pas à recopier »).

Protéomes : release 5 de Mycobrowser (mêmes conventions de locus tag que le panel affiché),
data/mycobrowser/proteomes/*.faa (header `>LOCUS|name|product|len AA` → id = LOCUS).

Méthode : pour chaque espèce cible, DIAMOND blastp bidirectionnel (H37Rv->cible et cible->H37Rv),
meilleur hit par requête (bitscore max, e-value<=1e-5), orthologue = paire réciproque (RBH).
On stocke le locus, l'identité % et l'e-value pour transparence.

Sortie : résultats/phase23_orthologs/orthologs.json (keyed Rv). Fusion directe dans les fiches
(idempotent) + enregistrement dans phase4. DIAMOND bundlé eggNOG-mapper.
Run: python analyses/phase23_orthologs.py
"""
from __future__ import annotations
import json, glob, subprocess, re
from pathlib import Path

AA_KEEP = re.compile(r"[^ACDEFGHIKLMNPQRSTVWYXBZUO]")

ROOT = Path(__file__).resolve().parent.parent
PROT = ROOT / "data" / "mycobrowser" / "proteomes"
OUT = ROOT / "résultats" / "phase23_orthologs"
GENES = ROOT / "site" / "content" / "genes"
DIAMOND = ROOT.parent / "L8" / "eggnog-mapper-2.1.12" / "eggnogmapper" / "bin" / "diamond"

# nom d'affichage -> fichier protéome cible
SPECIES = {
    "M. bovis": "M_bovis",
    "M. leprae": "M_leprae",
    "M. marinum": "M_marinum",
    "M. smegmatis": "M_smegmatis",
    "M. orygis": "M_orygis",
    "M. abscessus": "M_abscessus",
}
EVALUE = "1e-5"
SOURCE = "Reciprocal-best-hit DIAMOND vs Mycobrowser reference proteomes (release 5)"
REFS = [
    {"authors": "Buchfink B, Reuter K, Drost HG", "year": 2021,
     "title": "Sensitive protein alignments at tree-of-life scale using DIAMOND",
     "journal": "Nature Methods", "doi": "10.1038/s41592-021-01101-x"},
]


def clean_fasta(src: Path, dst: Path):
    """Réécrit un FASTA avec header = locus tag seul (1er champ avant '|')."""
    with open(src, encoding="utf-8", errors="replace") as fh, open(dst, "w") as out:
        for line in fh:
            if line.startswith(">"):
                loc = line[1:].split("|")[0].strip()
                out.write(f">{loc}\n")
            else:
                # sanitation : diamond refuse tout caractère hors alphabet AA (ex. '+')
                out.write(AA_KEEP.sub("", line.strip().upper()) + "\n")


def makedb(faa: Path, db: Path):
    subprocess.run([str(DIAMOND), "makedb", "--in", str(faa), "-d", str(db), "--quiet"], check=True)


def blastp(query: Path, db: Path, out_tsv: Path):
    subprocess.run([str(DIAMOND), "blastp", "-q", str(query), "-d", str(db), "-o", str(out_tsv),
                    "-e", EVALUE, "-k", "1", "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue",
                    "--quiet", "--threads", "8"], check=True)


def best_hits(tsv: Path):
    """qseqid -> (sseqid, pident, evalue), meilleur (1er, -k 1 = top)."""
    out = {}
    if not tsv.exists():
        return out
    for line in tsv.read_text().splitlines():
        p = line.split("\t")
        if len(p) < 4:
            continue
        q, s, pid, ev = p[0], p[1], p[2], p[3]
        if q not in out:  # -k 1 => 1 ligne/req, mais garde-fou
            out[q] = (s, float(pid), float(ev))
    return out


def main():
    print("== phase23 : orthologues par espèce, RBH DIAMOND (P5.5) ==")
    OUT.mkdir(parents=True, exist_ok=True)
    work = OUT / "work"
    work.mkdir(exist_ok=True)

    # protéome H37Rv nettoyé + db
    h37_clean = work / "H37Rv.clean.faa"
    clean_fasta(PROT / "H37Rv.faa", h37_clean)
    h37_db = work / "H37Rv_db"
    makedb(h37_clean, h37_db)

    # ensemble des rv de nos fiches (restreint la sortie)
    rv_set = set()
    for f in glob.glob(str(GENES / "*.json")):
        rv_set.add(json.load(open(f))["rv"])

    ortho = {rv: {} for rv in rv_set}
    for disp, key in SPECIES.items():
        tgt_clean = work / f"{key}.clean.faa"
        clean_fasta(PROT / f"{key}.faa", tgt_clean)
        tgt_db = work / f"{key}_db"
        makedb(tgt_clean, tgt_db)

        fwd = work / f"H37Rv_vs_{key}.tsv"   # H37Rv -> cible
        rev = work / f"{key}_vs_H37Rv.tsv"   # cible -> H37Rv
        blastp(h37_clean, tgt_db, fwd)
        blastp(tgt_clean, h37_db, rev)
        bf = best_hits(fwd)   # rv -> (tgt_locus, pid, ev)
        br = best_hits(rev)   # tgt_locus -> (rv, pid, ev)

        n = 0
        for rv, (tgt, pid, ev) in bf.items():
            back = br.get(tgt)
            if back and back[0] == rv and rv in ortho:   # réciprocité
                ortho[rv][disp] = {"locus": tgt, "identity": round(pid, 1), "evalue": ev}
                n += 1
        print(f"  {disp:14s} : {n} orthologues RBH")

    # ne garder que les rv avec >=1 orthologue
    rec = {rv: {"orthologs": d, "source": SOURCE, "refs": REFS} for rv, d in ortho.items() if d}
    (OUT / "orthologs.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))
    print(f"Écrit {OUT/'orthologs.json'} ({len(rec)} gènes avec >=1 orthologue)")

    # fusion dans les fiches
    n_written = n_hit = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        d["orthologs"] = rec.get(d["rv"], {})
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n_written += 1
        if d["rv"] in rec:
            n_hit += 1
    print(f"Fusionné 'orthologs' dans {n_written} fiches ({n_hit} avec orthologue)")

    # bilan couverture par espèce
    from collections import Counter
    cov = Counter()
    for r in rec.values():
        for sp in r["orthologs"]:
            cov[sp] += 1
    print("Couverture par espèce :", dict(cov))


if __name__ == "__main__":
    main()
