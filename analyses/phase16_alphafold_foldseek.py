#!/usr/bin/env python3
"""phase16_alphafold_foldseek.py -- couverture structurale génome-entière via AlphaFold DB (P3.2).

Les 350 modèles ESMFold (phase2c) ne couvraient que les gènes dark <=400 aa. AlphaFold DB
couvre TOUT le protéome H37Rv (y compris les longues protéines) avec des modèles de meilleure
qualité. On récupère le modèle AlphaFold de chaque gène (mapping via l'accession UniProt),
on le passe au MÊME Foldseek (vs tools/foldseek_db/pdb) et on stocke le résultat dans une
couche DISTINCTE (struct_af / plddt_af) : la couche ESMFold existante (cohérente avec le
manuscrit soumis) n'est PAS écrasée. Cette couche AF alimente ensuite P3.1 (M-CSA).

Route (cf. KB python-patterns 2026-07-01) : API `alphafold.ebi.ac.uk`, fichier
`AF-{acc}-F1-model_v6.pdb` (version courante = v6). pLDDT moyen = moyenne des B-factors CA.

Sorties : résultats/phase16_alphafold/{foldseek_af.json, plddt_af.json, missing.txt}
          + fusion struct_af/plddt_af dans site/content/genes/*.json
Modèles bruts : SCRATCH (non persistés dans le projet).
Run (fond) : python analyses/phase16_alphafold_foldseek.py
"""
from __future__ import annotations
import glob, json, os, subprocess, time, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "résultats" / "phase16_alphafold"
SCRATCH = Path(os.environ.get("AF_SCRATCH",
    "/tmp/claude-1000/-home-christophe-docs-codes-mtbc-annotation-mtbc/0af1e6b3-80ee-4e8d-9844-b83680f21762/scratchpad/af_models"))
FOLDSEEK = ROOT / "tools" / "foldseek" / "bin" / "foldseek"
FOLDSEEK_DB = ROOT / "tools" / "foldseek_db" / "pdb"
AF_URL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.pdb"
TOPN = 8
UA = "annotation_mtbc-research/1.0 (guyeux@univ-fcomte.fr)"


def rv_to_acc():
    m = {}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        acc = (d.get("uniprot") or {}).get("acc")
        if acc:
            m[d["rv"]] = acc
    return m


def download_models(mapping):
    SCRATCH.mkdir(parents=True, exist_ok=True)
    missing, ok = [], 0
    for i, (rv, acc) in enumerate(mapping.items(), 1):
        dest = SCRATCH / f"{rv}.pdb"
        if dest.exists() and dest.stat().st_size > 500:
            ok += 1; continue
        try:
            req = urllib.request.Request(AF_URL.format(acc=acc), headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
            if data[:6] in (b"HEADER", b"ATOM  ", b"REMARK", b"CRYST1") or b"\nATOM  " in data[:2000]:
                dest.write_bytes(data); ok += 1
            else:
                missing.append(f"{rv}\t{acc}\tbad_content")
        except Exception as e:
            missing.append(f"{rv}\t{acc}\t{type(e).__name__}")
        if i % 200 == 0:
            print(f"  [{i}/{len(mapping)}] téléchargés OK={ok} manquants={len(missing)}", flush=True)
        time.sleep(0.12)  # politesse EBI
    print(f"Téléchargement fini : OK={ok}, manquants={len(missing)}")
    return missing


def mean_plddt(pdb: Path):
    vals = []
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                pass
    if not vals:
        return None
    m = sum(vals) / len(vals)
    band = "very high" if m >= 90 else "confident" if m >= 70 else "low" if m >= 50 else "very low"
    return {"mean_plddt": round(m, 1), "band": band, "source": "AlphaFold DB v6"}


def run_foldseek_batch():
    tmp = SCRATCH.parent / "af_fs_tmp"
    qdb = SCRATCH.parent / "af_qdb"
    resdb = SCRATCH.parent / "af_res"
    m8 = OUT / "af.m8"
    subprocess.run([str(FOLDSEEK), "createdb", str(SCRATCH), str(qdb)], check=True,
                   capture_output=True, text=True)
    subprocess.run([str(FOLDSEEK), "search", str(qdb), str(FOLDSEEK_DB), str(resdb), str(tmp),
                    "-e", "10", "--max-seqs", "50", "-a"], check=True, capture_output=True, text=True)
    subprocess.run([str(FOLDSEEK), "convertalis", str(qdb), str(FOLDSEEK_DB), str(resdb), str(m8),
                    "--format-output", "query,target,evalue,prob,alntmscore,theader"],
                   check=True, capture_output=True, text=True)
    hits = {}
    for line in m8.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 6:
            continue
        q = f[0].split(".pdb")[0]  # 'Rv0001.pdb' ou 'Rv0001.pdb_A' → 'Rv0001'
        ev = float(f[2])
        hits.setdefault(q, []).append({"target": f[1], "evalue": ev, "prob": float(f[3]),
                                       "tmscore": float(f[4]), "description": f[5],
                                       "significant": ev < 0.01})
    for q in hits:
        hits[q].sort(key=lambda h: h["prob"], reverse=True)
        hits[q] = hits[q][:TOPN]
    return hits


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mapping = rv_to_acc()
    print(f"{len(mapping)} gènes avec accession UniProt → modèle AlphaFold", flush=True)
    missing = download_models(mapping)
    (OUT / "missing.txt").write_text("\n".join(missing))

    # pLDDT par modèle
    plddt = {}
    for rv in mapping:
        p = SCRATCH / f"{rv}.pdb"
        if p.exists() and p.stat().st_size > 500:
            mp = mean_plddt(p)
            if mp:
                plddt[rv] = mp
    (OUT / "plddt_af.json").write_text(json.dumps(plddt, ensure_ascii=False, indent=1))
    print(f"pLDDT calculé pour {len(plddt)} modèles", flush=True)

    print("Foldseek batch (createdb + search + convertalis)…", flush=True)
    hits = run_foldseek_batch()
    (OUT / "foldseek_af.json").write_text(json.dumps(hits, ensure_ascii=False, indent=1))
    n_sig = sum(1 for h in hits.values() if any(x["significant"] for x in h))
    print(f"Foldseek AF : {len(hits)} gènes avec hit(s), {n_sig} avec hit significatif (E<0.01)", flush=True)

    # gating pLDDT 70 + fusion dans les fiches (couche struct_af / plddt_af)
    n_written = n_gated = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f)); rv = d["rv"]
        pl = plddt.get(rv, {})
        gated = bool(pl and pl.get("mean_plddt", 0) >= 70)
        d["plddt_af"] = pl
        d["struct_af"] = {"hits": hits.get(rv, []), "plddt_gated": gated}
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n_written += 1
        if gated:
            n_gated += 1
    print(f"Fusionné struct_af/plddt_af dans {n_written} fiches ; {n_gated} avec modèle pLDDT>=70", flush=True)


if __name__ == "__main__":
    main()
