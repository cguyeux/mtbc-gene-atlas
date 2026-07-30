#!/usr/bin/env python3
"""phase57_ntm_orthology.py -- P10.2 : présence/absence des gènes MTBC à travers le genre (NTM).

Répond à la question biologique que le dN/dS canettii ne peut pas (sous-puissant) : chaque gène MTBC est-il
un CŒUR ANCIEN (présent à travers les NTM) ou une INNOVATION SPÉCIFIQUE-MTBC (absente des NTM, y compris des
plus proches parents MTBAP = candidat facteur d'adaptation à l'hôte) ? Robuste car fondé sur la présence/absence
par tblastn (protéine MTBC vs génome NTM traduit 6 cadres), PAS sur un dN/dS ni un mapping H37Rv biaisé.

Garde-fous (KB) : tblastn protéine-vs-génome (divergent != absent) ; seuil de présence sur identité + COUVERTURE
de requête ; les plus proches parents (MTBAP) sont les discriminants du « spécifique-MTBC ». Un modèle nul (les
spécifiques-MTBC ne sont pas juste des ORF courts/low-complexity) sera vérifié à l'analyse.

[2026-07-30] TEST DE SENSIBILITÉ INTRINSÈQUE (correctif d'outil). Le seuil MIN_QCOV rate les homologues
PRÉSENTS MAIS DIVERGENTS : un gène peut n'aligner que sur un fragment de sa longueur et être classé
« MTBC-specific » à tort. Vécu sur Rv2438A, où un homologue de M. decipiens (79 % identité, 32 % couverture)
avait été manqué par qcov>=50 ; le faux « spécifique-MTBC » a fondé tout le cadrage d'un projet avant d'être
réfuté a posteriori par phase61. Le problème était que phase61 (test de sensibilité) est une phase SÉPARÉE,
lancée après coup et seulement sur une shortlist : la classification voyageait donc sans son test.
Correctif : les hits SOUS-SEUIL sont désormais conservés au lieu d'être jetés (ils sont déjà calculés, coût
nul) et attachés à chaque gène via les colonnes `sensitivity_flag` / `subthreshold_best`. Aucun gène ne peut
plus sortir « MTBC-specific » sans que son test de sensibilité voyage avec lui.
LIMITE assumée : la détection sous-seuil est bornée par l'e-value du tblastn de présence (1e-5). phase61
interroge plus permissivement (1e-3) et reste donc utile sur les gènes marqués `clean-at-evalue-1e-5`.
L'e-value de la passe principale est volontairement INCHANGÉE pour ne pas modifier `presence` (donc les
classifications déjà publiées) : ce correctif ajoute de l'information, il n'en réécrit aucune.

Sortie : résultats/phase57_ntm/ntm_presence_per_gene.tsv (+ matrice). ~55 s/génome × 58 ≈ 1 h.
Run: python analyses/phase57_ntm_orthology.py
"""
from __future__ import annotations
import glob, subprocess, sys, tempfile
from collections import defaultdict
from pathlib import Path

# progression VISIBLE même lancé en arrière-plan / piped (sinon stdout bloc-bufferisé, 0 octet jusqu'à la fin).
try:
    sys.stdout.reconfigure(line_buffering=True)   # type: ignore[attr-defined]  (stub TextIO incomplet ; vérifié OK à l'exécution)
except Exception:
    pass

ROOT = Path(__file__).resolve().parent.parent
HORS = ROOT.parent / "bdd" / "hors_mtbc"
OUT = ROOT / "résultats" / "phase57_ntm"
QUERY = OUT / "query_proteome.faa"

# seuils de présence d'un orthologue (permissifs en divergence, stricts en couverture)
MIN_PIDENT = 30.0
MIN_QCOV = 50.0
THREADS = 8

# groupes phylogénétiques (README bdd/hors_mtbc : Sapriel & Brosch, Khattak 2021)
MTBAP = {"M_decipiens", "M_lacus", "M_riyadhense", "M_shinjukuense"}  # clade le plus proche du MTBC
MKC = {"M_kansasii", "M_persicum", "M_pseudokansasii", "M_gastri", "M_attenuatum", "M_innocens",
       "M_basiliense", "M_neolactis", "M_novum"}  # complexe M. kansasii


def species_genomes() -> dict[str, Path]:
    """espèce -> un genome.fna (le premier trouvé)."""
    out = {}
    for f in sorted(glob.glob(str(HORS / "**" / "genome.fna"), recursive=True)):
        rel = Path(f).relative_to(HORS).parts
        if not rel:
            continue
        out.setdefault(rel[0], Path(f))
    return out


def tblastn_presence(genome: Path, tmp: Path) -> dict[str, tuple[float, float]]:
    """rv -> (pident, qcov) du meilleur hit, pour ce génome."""
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
    genomes = species_genomes()
    print(f"{len(genomes)} espèces NTM avec génome.")
    # presence[rv][species] = (pident, qcov)
    presence: dict[str, dict[str, tuple[float, float]]] = defaultdict(dict)
    # hits qui passent l'e-value de tblastn mais ÉCHOUENT le filtre pident/qcov : ne plus les jeter (cf. docstring).
    subthr: dict[str, dict[str, tuple[float, float]]] = defaultdict(dict)
    species = sorted(genomes)
    for i, sp in enumerate(species, 1):
        with tempfile.TemporaryDirectory() as td:
            hits = tblastn_presence(genomes[sp], Path(td))
        for rv, (pid, qcov) in hits.items():
            if pid >= MIN_PIDENT and qcov >= MIN_QCOV:
                presence[rv][sp] = (pid, qcov)
            else:
                subthr[rv][sp] = (pid, qcov)
        print(f"  [{i}/{len(species)}] {sp}: {sum(1 for rv in hits if hits[rv][0]>=MIN_PIDENT and hits[rv][1]>=MIN_QCOV)} orthologues")

    # tous les rv de la requête
    all_rv = [l[1:].strip() for l in open(QUERY) if l.startswith(">")]
    n_sp = len(species)
    with open(OUT / "ntm_presence_per_gene.tsv", "w") as fh:
        # colonnes de sensibilité AJOUTÉES EN FIN de ligne : les lecteurs par index (p[0..8]) restent valides.
        fh.write("rv\tn_ntm_present\tn_ntm_total\tfrac\tin_MTBAP\tin_MKC\tmean_pident\tclassification\t"
                 "species_present\tsensitivity_flag\tsubthreshold_best\n")
        n_flagged = 0
        for rv in all_rv:
            pres = presence.get(rv, {})
            npres = len(pres)
            frac = npres / n_sp
            in_mtbap = sum(1 for s in pres if s in MTBAP)
            in_mkc = sum(1 for s in pres if s in MKC)
            mean_pid = round(sum(v[0] for v in pres.values()) / npres, 1) if npres else 0.0
            if npres == 0:
                cls = "MTBC-specific"            # absent de tous les NTM = innovation candidate
            elif in_mtbap == 0 and frac < 0.15:
                cls = "MTBC-near-specific"        # absent des plus proches (MTBAP) + rare ailleurs
            elif frac >= 0.5:
                cls = "genus-core"
            else:
                cls = "restricted"
            # test de sensibilité attaché à la classification (cf. docstring) : un « absent » qui a des hits
            # sous-seuil est un PRÉSENT-MAIS-DIVERGENT candidat, pas une innovation.
            sub = subthr.get(rv, {})
            if npres == 0 and sub:
                # meilleur candidat divergent : on privilégie l'identité, en départageant par la couverture.
                best_sp = max(sub, key=lambda s: (sub[s][0], sub[s][1]))
                bp, bq = sub[best_sp]
                in_mtbap_sub = sum(1 for s in sub if s in MTBAP)
                flag = "SUBTHRESHOLD_HITS_MTBAP" if in_mtbap_sub else "SUBTHRESHOLD_HITS"
                sub_txt = f"{best_sp}:{bp:.1f}id/{bq:.0f}cov;n={len(sub)};mtbap={in_mtbap_sub}"
                n_flagged += 1
            elif npres == 0:
                flag, sub_txt = "clean-at-evalue-1e-5", "-"
            else:
                flag, sub_txt = "-", "-"
            fh.write(f"{rv}\t{npres}\t{n_sp}\t{frac:.3f}\t{in_mtbap}\t{in_mkc}\t{mean_pid}\t{cls}\t"
                     f"{','.join(sorted(pres))}\t{flag}\t{sub_txt}\n")
    from collections import Counter
    cls_counts = Counter()
    for rv in all_rv:
        pres = presence.get(rv, {})
        npres = len(pres); frac = npres / n_sp
        in_mtbap = sum(1 for s in pres if s in MTBAP)
        cls_counts["MTBC-specific" if npres == 0 else
                   "MTBC-near-specific" if (in_mtbap == 0 and frac < 0.15) else
                   "genus-core" if frac >= 0.5 else "restricted"] += 1
    print("classification :", dict(cls_counts))
    if n_flagged:
        print(f"\n  /!\\ SENSIBILITÉ : {n_flagged}/{cls_counts['MTBC-specific']} gènes « MTBC-specific » ont des hits "
              f"SOUS-SEUIL chez au moins un NTM (colonne sensitivity_flag).")
        print("      Ce sont des PRÉSENTS-MAIS-DIVERGENTS candidats, PAS des innovations MTBC : ne pas les cadrer")
        print("      comme spécifiques-MTBC sans avoir inspecté subthreshold_best (a fortiori si flag=...._MTBAP,")
        print("      c'est-à-dire hit chez un des 4 plus proches parents).")
    print(f"Table : {OUT/'ntm_presence_per_gene.tsv'}")


if __name__ == "__main__":
    main()
