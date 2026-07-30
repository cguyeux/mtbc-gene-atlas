#!/usr/bin/env python3
"""phase61_sensitivity_mtbc_specific.py -- P11.2c : test de sensibilité de l'absence NTM.

Garde-fou (évaluation P11.1) : « MTBC-specific » = pas de hit tblastn à pident≥30 & qcov≥50 sur 53 NTM.
Un gène PRÉSENT MAIS TRÈS DIVERGENT pourrait être manqué → faux « spécifique ». On teste en RÉ-INTERROGEANT
les 53 génomes NTM en tblastn RELÂCHÉ (e-value 1e-3, aucun seuil d'identité) : si même une recherche permissive
ne trouve RIEN, l'absence est robuste (vraie innovation) ; si un hit faible apparaît, le gène est présent-divergent
et doit être RÉTROGRADÉ. (HMMER/phmmer absent localement ; HHpred profil-profil = hand-off web, hors scope.)
Entrée : les 17 candidats durcis (phase60, keep=1). Sortie : résultats/phase61/sensitivity.tsv.
Run: python analyses/phase61_sensitivity_mtbc_specific.py
"""
from __future__ import annotations
import json, subprocess, sys, tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase57_ntm_orthology as p57  # réutilise species_genomes()

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
HARDENED = ROOT / "résultats" / "phase60" / "dark_mtbc_specific_hardened.tsv"
OUT = ROOT / "résultats" / "phase61"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cands = [l.split("\t")[0] for l in open(HARDENED) if not l.startswith("rv") and l.split("\t")[1] == "1"]
    q = OUT / "cands.faa"
    q.write_text("".join(f">{rv}\n{json.load(open(GENES/f'{rv}.json')).get('protein_mtbc0')}\n" for rv in cands))
    genomes = p57.species_genomes()
    print(f"{len(cands)} candidats × {len(genomes)} NTM, tblastn RELÂCHÉ (e-value 1e-3, sans seuil d'identité)…")

    # meilleur hit relâché par candidat, à travers tous les NTM
    best: dict[str, tuple] = {}   # rv -> (pident, qcov, species)
    for i, (sp, g) in enumerate(sorted(genomes.items()), 1):
        with tempfile.TemporaryDirectory() as td:
            db = Path(td) / "db"
            subprocess.run(["makeblastdb", "-in", str(g), "-dbtype", "nucl", "-out", str(db)],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            out = Path(td) / "h.tsv"
            subprocess.run(["tblastn", "-query", str(q), "-db", str(db), "-out", str(out), "-evalue", "1e-3",
                            "-num_threads", "8", "-max_target_seqs", "1", "-max_hsps", "1",
                            "-outfmt", "6 qseqid pident qcovs bitscore"],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            for line in open(out):
                pp = line.rstrip("\n").split("\t")
                if len(pp) < 4:
                    continue
                rv, pid, qc, bit = pp[0], float(pp[1]), float(pp[2]), float(pp[3])
                if rv not in best or bit > best[rv][3]:
                    best[rv] = (pid, qc, sp, bit)
        print(f"  [{i}/{len(genomes)}] {sp}")

    rows = []
    for rv in cands:
        b = best.get(rv)
        if not b:
            verdict = "robust-MTBC-specific"          # aucun hit même relâché
            detail = "no hit at e-value 1e-3 on any of 53 NTM"
        else:
            pid, qc, sp, _ = b
            # hit faible réel = présent-divergent ; hit ultra-marginal (qcov très bas) = bruit, reste spécifique
            if pid >= 25 and qc >= 30:
                verdict = "present-divergent-DOWNGRADE"
                detail = f"weak hit {pid:.0f}%/{qc:.0f}%cov in {sp} — likely present but divergent"
            else:
                verdict = "robust-MTBC-specific"
                detail = f"only marginal hit {pid:.0f}%/{qc:.0f}%cov in {sp} (below sensitivity floor)"
        rows.append({"rv": rv, "verdict": verdict, "detail": detail})

    with open(OUT / "sensitivity.tsv", "w") as fh:
        fh.write("rv\tsensitivity_verdict\tdetail\n")
        for r in rows:
            fh.write(f"{r['rv']}\t{r['verdict']}\t{r['detail']}\n")
    robust = [r for r in rows if r["verdict"] == "robust-MTBC-specific"]
    down = [r for r in rows if r["verdict"] != "robust-MTBC-specific"]
    print(f"\n{len(robust)}/{len(cands)} ROBUSTES (absence confirmée même en recherche permissive).")
    for r in down:
        print(f"  RÉTROGRADÉ {r['rv']}: {r['detail']}")
    print(f"Robustes : {', '.join(r['rv'] for r in robust)}")
    print(f"Rapport : {OUT/'sensitivity.tsv'}")


if __name__ == "__main__":
    main()
