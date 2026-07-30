#!/usr/bin/env python3
"""phase18a_run_deeptmhmm.py -- exécuter DeepTMHMM (référence) sur le protéome, en batches (P3.7).

DeepTMHMM (Hallgren et al. 2022) est le prédicteur de RÉFÉRENCE de topologie transmembranaire + peptide
signal. Il tourne en cloud biolib **ANONYMEMENT** (pas de licence ni de compte) via `pybiolib`. PIÈGE
découvert 2026-07-01 : l'upload biolib bascule sur un chemin *multipart* au-dessus de ~100-300 séquences,
et ce chemin est CASSÉ en anonyme (« invalid type: string "" , expected struct MultipartUploadStartResponse »).
=> il FAUT batcher à **<= 100 séquences/job** (100 validé OK, 300 KO). Ce script batche, est reprenable
(saute les batches déjà produits), réessaie (3x) et concatène les `.3line`.

DIAGNOSTIC (pourquoi ce runner et pas un naïf) — deux pièges biolib appris à la dure 2026-07-01 :
  1. ÉCHEC SILENCIEUX au-dessus de ~100-300 séq : `app.cli()` ne LÈVE PAS d'exception ; le job revient
     `status=failed` et `save_files()` ne crée simplement pas `predicted_topologies.3line`. => TOUJOURS
     vérifier l'EXISTENCE du fichier de sortie (pas seulement l'absence d'exception). D'où batch<=100.
  2. UPLOAD RELATIF au cwd : biolib uploade les fichiers passés en `--fasta` relativement au répertoire
     courant. Passer un chemin absolu casse le mapping côté cloud (`FileNotFoundError: 'hash/xxx.fasta'`).
     => on `os.chdir(work)` et on passe seulement `fa.name`.
  3. QUOTA anonyme : ~50 jobs/fenêtre puis tout échoue `status=failed` → run par vagues (reprenable) ou compte.

Prérequis : `pip install pybiolib` (pur Python, OK sous Python 3.14, contrairement à tmhmm.py qui ne compile pas).
Entrée : un FASTA du protéome (defaut : construit depuis site/content/genes/*.json, champ protein_mtbc0).
Sortie : résultats/phase18_localization/predicted_topologies.3line (concaténé) — consommé par phase18b.
Run : python analyses/phase18a_run_deeptmhmm.py [--batch 100] [--fasta proteome.fasta] [--workdir <dir>]
Note : long (protéome entier = ~40 batches). Lancer en tâche de fond.
"""
from __future__ import annotations
import argparse, glob, json, math, os, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUTCAT = ROOT / "résultats" / "phase18_localization" / "predicted_topologies.3line"


def build_fasta(path: Path):
    n = 0
    with open(path, "w") as fh:
        for f in sorted(glob.glob(str(GENES / "*.json"))):
            d = json.load(open(f)); s = d.get("protein_mtbc0") or ""
            if s:
                fh.write(f">{d['rv']}\n{s}\n"); n += 1
    return n


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=int, default=100, help="taille de batch (<=100 pour éviter le multipart cassé)")
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--workdir", default=str(ROOT / "résultats" / "phase18_localization" / "dt_batches"))
    a = ap.parse_args(argv)
    import biolib  # noqa: E402

    # Auth via token stocké (~/.config/biolib/token) si présent : lève le plafond multipart + le quota anonyme.
    tokf = Path.home() / ".config" / "biolib" / "token"
    if tokf.exists():
        biolib.set_api_token(tokf.read_text().strip())
        print("biolib: authentifié via ~/.config/biolib/token", flush=True)

    work = Path(a.workdir); work.mkdir(parents=True, exist_ok=True)
    fasta = Path(a.fasta) if a.fasta else work / "proteome.fasta"
    if not a.fasta:
        print("Construction FASTA:", build_fasta(fasta), "protéines", flush=True)

    recs, cur = [], None
    for line in open(fasta):
        if line.startswith(">"):
            cur = [line.strip(), ""]; recs.append(cur)
        else:
            cur[1] += line.strip()
    n = len(recs); nb = math.ceil(n / a.batch)
    print(f"{n} protéines, {nb} batches de {a.batch}", flush=True)

    app = biolib.load("DTU/DeepTMHMM")
    for b in range(nb):
        out3 = work / f"b{b:02d}.3line"
        if out3.exists() and out3.stat().st_size > 0:
            continue
        chunk = recs[b * a.batch:(b + 1) * a.batch]
        fa = work / f"b{b:02d}.fasta"
        fa.write_text("\n".join(f"{h}\n{s}" for h, s in chunk) + "\n")
        ok = False
        for attempt in range(3):
            cwd = os.getcwd()
            try:
                os.chdir(work)  # biolib uploade ET sauvegarde relativement au cwd
                job = app.cli(args=f"--fasta {fa.name}")
                # Passer un chemin RELATIF : un chemin ABSOLU se fait NICHER par biolib
                # (work/work/oNN/...), et le fichier attendu n'est alors pas là où on le cherche.
                job.save_files(f"o{b:02d}")
            except Exception as e:  # noqa: BLE001
                os.chdir(cwd)
                print(f"b{b} try{attempt}: {type(e).__name__} {str(e)[:80]}", flush=True)
                time.sleep(5)
                continue
            os.chdir(cwd)
            # Robuste au nesting : retrouver le 3-line où qu'il ait atterri sous work.
            hits = [p for p in glob.glob(str(work / "**" / "predicted_topologies.3line"), recursive=True)
                    if os.path.getsize(p) > 0]
            if hits:
                os.replace(hits[0], out3); ok = True; break
            print(f"b{b} try{attempt}: pas de predicted_topologies.3line trouvé", flush=True)
            time.sleep(5)
        print(f"b{b:02d}: {'OK' if ok else 'FAIL'} ({len(glob.glob(str(work/'b*.3line')))}/{nb})", flush=True)

    # Write the concatenated 3-line INTO the workdir (self-contained): a targeted run on a
    # custom --workdir must NOT clobber the whole-proteome predicted_topologies.3line. The
    # legacy global path is mirrored ONLY for the default (whole-proteome) workdir.
    parts = sorted(glob.glob(str(work / "b*.3line")))
    outcat = work / "predicted_topologies.3line"
    with open(outcat, "w") as fh:
        for p in parts:
            fh.write(open(p).read())
    default_work = ROOT / "résultats" / "phase18_localization" / "dt_batches"
    if work.resolve() == default_work.resolve():
        OUTCAT.parent.mkdir(parents=True, exist_ok=True)
        import shutil
        shutil.copy(outcat, OUTCAT)
    print(f"=== DONE {len(parts)}/{nb} batches -> {outcat} ===", flush=True)


if __name__ == "__main__":
    main()
