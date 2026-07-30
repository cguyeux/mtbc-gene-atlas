#!/usr/bin/env python3
"""phase17_pfam_relaxed.py -- domaines Pfam TENTATIFS (sous le seuil de gathering) pour les gènes dark (P3.3).

Substitut tractable au HHpred (bloqué : pdb70 28,6 Go + UniRef30 ~50 Go ne tiennent pas sur le disque, et
un HHpred vs pdb70 serait redondant avec le Foldseek-vs-PDB génome-entier déjà fait en P3.2). Cible : les
~172 hypothétiques encore SANS aucun signal (ni Pfam --cut_ga, ni eggNOG, ni UniProt, ni structure, ni M-CSA).

Idée : re-scanner ces gènes contre Pfam-A SANS `--cut_ga` (que phase2b applique pour la précision) avec un
seuil E-value permissif. On récupère des correspondances de domaine SOUS le seuil de gathering : ce sont des
PISTES de fonction en basse confiance, clairement étiquetées comme telles (pas des assignations fermes).

Sortie : résultats/phase17_pfam_relaxed/pfam_relaxed.json {rv: {pfam_name, pfam_acc, i_evalue, ali_from, ali_to, description}}
         + fusion du champ `pfam_tentative` dans les fiches dark ciblées.
Réutilise hmmscan + Pfam-A de phase2b (env HMMSCAN_BIN / PFAM_DB). Stdlib.
Run: python analyses/phase17_pfam_relaxed.py [--write]
"""
from __future__ import annotations
import argparse, glob, json, os, subprocess, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "résultats" / "phase17_pfam_relaxed"
HMMSCAN = Path(os.environ.get("HMMSCAN_BIN", ROOT.parent / "L8" / "eggnog-mapper-2.1.12" / "eggnogmapper" / "bin" / "hmmscan"))
PFAM_DB = Path(os.environ.get("PFAM_DB", ROOT.parent / "projets_abandonnes" / "Mycobacterium_sp_novel" / "data" / "pfam" / "Pfam-A.hmm"))
KEEP_IEVAL = 0.1   # seuil « piste tentative » (bien au-dessus du bruit, sous le gathering strict)
ne = lambda x: bool(x) and str(x).strip() not in ("", "[]", "{}")


def residual_dark():
    """Hypothétiques sans AUCUN signal (Pfam strict, eggNOG, UniProt, structure, M-CSA)."""
    out = {}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        if "hypothetical" not in (d.get("product_h37rv") or "").lower():
            continue
        egg = d.get("eggnog") or {}
        e_handle = any(ne(egg.get(k)) for k in ("cog_cat", "ec", "kegg_ko"))
        up = ne((d.get("uniprot") or {}).get("function"))
        dom = bool(d.get("domains"))
        sh = any(h.get("significant") for h in (d.get("struct_hits") or []))
        sa = d.get("struct_af") or {}
        sa_hit = sa.get("plddt_gated") and any(h.get("significant") for h in sa.get("hits", []))
        mcsa = bool((d.get("mcsa") or {}).get("verdict"))
        if not (e_handle or up or dom or sh or sa_hit or mcsa):
            seq = d.get("protein_mtbc0") or ""
            if seq:
                out[d["rv"]] = seq
    return out


def run_hmmscan(fasta: Path, domtbl: Path):
    # SANS --cut_ga : seuil E-value permissif pour capter les hits sous le gathering.
    cmd = [str(HMMSCAN), "-E", "10", "--domE", "10", "--domtblout", str(domtbl), str(PFAM_DB), str(fasta)]
    subprocess.run(cmd, capture_output=True, text=True, check=True)


def parse_best(domtbl: Path):
    best = {}
    for line in domtbl.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 23:
            continue
        pfam_name, pfam_acc, rv = f[0], f[1], f[3]
        try:
            ieval = float(f[12]); afrom, ato = int(f[17]), int(f[18])
        except ValueError:
            continue
        desc = " ".join(f[22:])
        if ieval > KEEP_IEVAL:
            continue
        if rv not in best or ieval < best[rv]["i_evalue"]:
            best[rv] = {"pfam_name": pfam_name, "pfam_acc": pfam_acc, "i_evalue": ieval,
                        "ali_from": afrom, "ali_to": ato, "description": desc,
                        "note": "tentative: below Pfam gathering threshold (--cut_ga); low-confidence lead"}
    return best


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)

    dark = residual_dark()
    print(f"Gènes dark ciblés : {len(dark)}")
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as fh:
        for rv, seq in dark.items():
            fh.write(f">{rv}\n{seq}\n")
        faa = Path(fh.name)
    domtbl = OUT / "dark_relaxed.domtblout"
    print("hmmscan (sans --cut_ga)…")
    run_hmmscan(faa, domtbl)
    best = parse_best(domtbl)
    (OUT / "pfam_relaxed.json").write_text(json.dumps(best, ensure_ascii=False, indent=1))
    print(f"Pistes de domaine tentatives (i-Evalue < {KEEP_IEVAL}) : {len(best)}/{len(dark)}")
    for rv, h in sorted(best.items(), key=lambda kv: kv[1]["i_evalue"])[:15]:
        print(f"  {rv}  {h['pfam_name']} ({h['pfam_acc']})  iE={h['i_evalue']:.1e}  {h['description'][:48]}")

    if a.write:
        n = 0
        for f in glob.glob(str(GENES / "*.json")):
            d = json.load(open(f))
            if d["rv"] in best:
                d["pfam_tentative"] = best[d["rv"]]
                Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
                n += 1
        print(f"\n(--write) champ pfam_tentative écrit dans {n} fiches")
    faa.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
