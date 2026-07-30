#!/usr/bin/env python3
"""phase66_phylostratum_scan.py -- P10.3 : phylostratigraphie CROSS-GENRE (hors genre Mycobacterium).

Le cran au-dessus de phase57 (présence/absence DANS le genre). phase57 laisse les 3037 « genus-core »
indifférenciés : un gène qui s'arrête au genre Mycobacterium et un gène conservé jusqu'à E. coli sont
dans le même sac. Ce scan monte l'échelle taxonomique avec les 13 génomes de `bdd/hors_mycobacterium/`
(Mycobacteriaceae > Corynebacteriales > autres Actinomycetia > hors-phylum), et retient pour chaque gène
la PROFONDEUR la plus grande où tblastn le détecte encore.

Méthode identique à phase57 (KB) : tblastn protéine-vs-génome traduit 6 cadres, présence = pident>=30 &
qcov>=50, meilleur hit par requête. Réutilise le proteome `résultats/phase57_ntm/query_proteome.faa`.

GARDE-FOU null-first (KB [2026-07-11], Moyers & Zhang) : la détection tblastn décroît avec la longueur du
gène → un gène court paraît « jeune » par simple échec de détection. On ENREGISTRE la longueur par gène
(colonne len_aa) pour le contrôle de modèle nul fait en phase67 ; ici on ne fait que scanner.

Sortie : résultats/phase66_phylostratum/crossgenus_presence_per_gene.tsv (lecture SEULE des génomes).
Run: python analyses/phase66_phylostratum_scan.py   (~13 génomes × ~1 min)
"""
from __future__ import annotations
import glob, json, subprocess, sys, tempfile
from collections import defaultdict
from pathlib import Path

try:
    sys.stdout.reconfigure(line_buffering=True)   # type: ignore[attr-defined]  (stub TextIO incomplet ; vérifié OK à l'exécution)   # progression visible en arrière-plan (cf. KB)
except Exception:
    pass

ROOT = Path(__file__).resolve().parent.parent
HORS = ROOT.parent / "bdd" / "hors_mycobacterium"
OUT = ROOT / "résultats" / "phase66_phylostratum"
QUERY = ROOT / "résultats" / "phase57_ntm" / "query_proteome.faa"   # même proteome que P10.2

MIN_PIDENT = 30.0
MIN_QCOV = 50.0
THREADS = 8

# rang de PROFONDEUR taxonomique (plus grand = plus ancien/large), à partir du champ `clade` des meta.json.
# 1 = famille Mycobacteriaceae (genre sœur) ; 2 = ordre Corynebacteriales ; 3 = classe Actinomycetia
# (autres ordres) ; 4 = hors-phylum (Bacteria au sens large, témoins Proteobacteria/Firmicutes).
CLADE_RANK = {
    "Mycobacteriaceae(rapid)": 1,
    "Corynebacteriales": 2,
    "Bifidobacteriales": 3, "Propionibacteriales": 3, "Micrococcales": 3, "Streptomycetales": 3,
    "control_Proteobacteria": 4, "control_Firmicutes": 4,
}
RANK_DEEPEST_LABEL = {0: "Mycobacterium-or-shallower", 1: "Mycobacteriaceae",
                      2: "Corynebacteriales", 3: "Actinomycetia", 4: "Bacteria"}


def load_genomes() -> dict[str, tuple[Path, int, str]]:
    """label -> (genome.fna, rang, organism). Lit le clade dans ref/meta.json."""
    out = {}
    for meta in sorted(glob.glob(str(HORS / "*" / "ref" / "meta.json"))):
        d = json.load(open(meta))
        label = Path(meta).parents[1].name
        genome = Path(meta).parent / "genome.fna"
        clade = d.get("clade", "?")
        rank = CLADE_RANK.get(clade)
        if rank is None:
            print(f"  ATTENTION clade inconnu '{clade}' pour {label} — ignoré")
            continue
        if not genome.exists():
            print(f"  ATTENTION genome.fna absent pour {label} — ignoré")
            continue
        out[label] = (genome, rank, d.get("organism", label))
    return out


def query_lengths() -> dict[str, int]:
    """rv -> longueur (aa), depuis le proteome de requête."""
    lens, rv, n = {}, None, 0
    for line in open(QUERY):
        if line.startswith(">"):
            if rv is not None:
                lens[rv] = n
            rv = line[1:].split()[0]
            n = 0
        else:
            n += len(line.strip())
    if rv is not None:
        lens[rv] = n
    return lens


def tblastn_presence(genome: Path, tmp: Path) -> dict[str, tuple[float, float]]:
    """rv -> (pident, qcov) du meilleur hit (par bitscore) dans ce génome."""
    db = tmp / "db"
    subprocess.run(["makeblastdb", "-in", str(genome), "-dbtype", "nucl", "-out", str(db)],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out = tmp / "hits.tsv"
    subprocess.run(["tblastn", "-query", str(QUERY), "-db", str(db), "-out", str(out),
                    "-evalue", "1e-5", "-num_threads", str(THREADS), "-max_target_seqs", "1",
                    "-max_hsps", "1", "-outfmt", "6 qseqid pident qcovs bitscore"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    best: dict[str, tuple[float, float, float]] = {}
    for line in open(out):
        p = line.rstrip("\n").split("\t")
        if len(p) < 4:
            continue
        rv, pid, qcov, bit = p[0], float(p[1]), float(p[2]), float(p[3])
        if rv not in best or bit > best[rv][2]:
            best[rv] = (pid, qcov, bit)
    return {rv: (v[0], v[1]) for rv, v in best.items()}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    genomes = load_genomes()
    lens = query_lengths()
    print(f"{len(genomes)} génomes hors-genre, {len(lens)} protéines de requête.")
    # presence[rv][label] = (pident, qcov)
    presence: dict[str, dict[str, tuple[float, float]]] = defaultdict(dict)
    labels = sorted(genomes)
    for i, lab in enumerate(labels, 1):
        genome, rank, org = genomes[lab]
        with tempfile.TemporaryDirectory() as td:
            hits = tblastn_presence(genome, Path(td))
        k = 0
        for rv, (pid, qcov) in hits.items():
            if pid >= MIN_PIDENT and qcov >= MIN_QCOV:
                presence[rv][lab] = (pid, qcov)
                k += 1
        print(f"  [{i}/{len(labels)}] {lab} (rang {rank}, {org}): {k} orthologues détectés")

    all_rv = list(lens)
    with open(OUT / "crossgenus_presence_per_gene.tsv", "w") as fh:
        fh.write("rv\tlen_aa\txg_n_hits\txg_deepest_rank\txg_deepest_label\t"
                 "n_r1_family\tn_r2_coryneb\tn_r3_actino\tn_r4_bacteria\tmean_pident\tgenomes_hit\n")
        for rv in all_rv:
            pres = presence.get(rv, {})
            ranks = [genomes[lab][1] for lab in pres]
            deepest = max(ranks) if ranks else 0
            by = {1: 0, 2: 0, 3: 0, 4: 0}
            for r in ranks:
                by[r] += 1
            mean_pid = round(sum(v[0] for v in pres.values()) / len(pres), 1) if pres else 0.0
            fh.write(f"{rv}\t{lens[rv]}\t{len(pres)}\t{deepest}\t{RANK_DEEPEST_LABEL[deepest]}\t"
                     f"{by[1]}\t{by[2]}\t{by[3]}\t{by[4]}\t{mean_pid}\t{','.join(sorted(pres))}\n")

    from collections import Counter
    dc = Counter(max([genomes[lab][1] for lab in presence.get(rv, {})] or [0]) for rv in all_rv)
    print("\nprofondeur cross-genre atteinte (rang -> #gènes) :")
    for r in sorted(dc):
        print(f"  rang {r} ({RANK_DEEPEST_LABEL[r]}): {dc[r]}")
    # contrôles positifs attendus (housekeeping profonds) : doivent atteindre le rang 4
    for probe in ("Rv0667", "Rv0440", "Rv0006"):   # rpoB, groEL2, gyrA
        pres = presence.get(probe, {})
        deep = max([genomes[lab][1] for lab in pres] or [0])
        print(f"  [contrôle] {probe}: rang {deep} ({RANK_DEEPEST_LABEL[deep]}), {len(pres)}/{len(labels)} génomes")
    print(f"\nTable : {OUT/'crossgenus_presence_per_gene.tsv'}")


if __name__ == "__main__":
    main()
