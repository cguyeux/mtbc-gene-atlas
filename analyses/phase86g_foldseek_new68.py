#!/usr/bin/env python3
"""phase86g — AlphaFold-model Foldseek structural layer for the 68 Mycobrowser genes (P16.5d-cont).

Correction: Foldseek was wrongly marked "blocked (no GPU/structures)". No GPU is needed — the
query structures come from the AlphaFold DB by UniProt accession (the 68 all have one, from
Mycobrowser), and the Foldseek target DB is the local PDB DB (tools/foldseek_db/pdb). Verified:
7/8 sampled genes have an AF model_v6.

Reuses phase16.{download_models, run_foldseek_batch, mean_plddt} on a DEDICATED 68-only model dir
(so the 3878 existing models are not re-searched) and merges struct_af / plddt_af into the 68
fiches only, with phase16's pLDDT<70 gating. A Foldseek hit is a fold match, NOT proof of
function; it never changes a verdict here (M-CSA active-site check would be needed). Idempotent.
"""
from __future__ import annotations
import importlib.util
import json, glob, os
from pathlib import Path

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CG = os.path.join(ROOT, "site", "content", "genes")
MODELDIR = os.path.join(ROOT, "résultats", "phase86g_foldseek", "af_models")
OUTDIR = os.path.join(ROOT, "résultats", "phase86g_foldseek")


def main():
    spec = importlib.util.spec_from_file_location(
        "p16", os.path.join(ROOT, "analyses", "phase16_alphafold_foldseek.py"))
    p16 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(p16)
    # redirect the model/output dirs to a dedicated 68-only space
    Path(MODELDIR).mkdir(parents=True, exist_ok=True)
    Path(OUTDIR).mkdir(parents=True, exist_ok=True)
    p16.SCRATCH = Path(MODELDIR)
    p16.OUT = Path(OUTDIR)

    # rv -> UniProt acc, restricted to the 68
    mapping = {}
    for f in glob.glob(os.path.join(CG, "Rv*.json")):
        d = json.load(open(f))
        if not str(d.get("curation_note", "")).startswith("Added P16.5d"):
            continue
        acc = (d.get("uniprot") or {}).get("acc")
        if acc:
            mapping[d["rv"]] = acc
    print(f"phase86g: {len(mapping)} gènes avec accession UniProt")

    missing = p16.download_models(mapping)
    hits = p16.run_foldseek_batch()

    written = n_gated = n_sig = 0
    for rv in mapping:
        model = Path(MODELDIR) / f"{rv}.pdb"
        pl = p16.mean_plddt(model) if model.exists() else {}
        gated = bool(pl and pl.get("mean_plddt", 0) >= 70)
        fn = os.path.join(CG, f"{rv}.json")
        d = json.load(open(fn))
        d["plddt_af"] = pl or {}
        d["struct_af"] = {"hits": hits.get(rv, []), "plddt_gated": gated}
        json.dump(d, open(fn, "w"), ensure_ascii=False, indent=2)
        written += 1
        n_gated += gated
        if gated and any(h["significant"] for h in hits.get(rv, [])):
            n_sig += 1

    print(f"phase86g: struct_af/plddt_af écrits sur {written}/{len(mapping)} gènes "
          f"(modèles manquants AFDB: {len(missing)})")
    print(f"  pLDDT>=70 (modèle fiable): {n_gated} | avec hit Foldseek significatif ET gated: {n_sig}")
    for rv in mapping:
        d = json.load(open(os.path.join(CG, f"{rv}.json")))
        sa = d.get("struct_af", {})
        if d.get("verdict") == "dark" and sa.get("plddt_gated"):
            sig = [h for h in sa.get("hits", []) if h["significant"]]
            if sig:
                print(f"    [dark] {rv} -> {sig[0]['description'][:55]} (tm={sig[0]['tmscore']:.2f})")


if __name__ == "__main__":
    main()
