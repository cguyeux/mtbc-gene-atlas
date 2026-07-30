#!/usr/bin/env python3
"""phase86e — intra-MTBC conservation (pN/pS + pseudogene) for the 68 Mycobrowser genes (P16.5d-cont).

The 68 genes added in phase86b are absent from the RefSeq CDS fasta that phase2f reads,
so phase2f skips them. This reuses phase2f's *validated* pass-2 logic (codon-aware pN/pS,
disruption load, coordinate auto-detection) WITHOUT re-scanning 145k strains: the heavy
pass-1 cache (résultats/phase2f_conservation/snp_counts.tsv, indel_counts.tsv, meta.json)
already exists. We build an AUGMENTED CDS fasta = the original + the 68 (nucleotide CDS
extracted from the H37Rv genome at their H37Rv coordinates), monkey-patch phase2f.CDS_FA to
it, and run phase2f.main() (which reuses the cache and merges `conservation` into every fiche,
now including the 68). Idempotent; re-affirms the 3906's conservation unchanged.
"""
from __future__ import annotations
import importlib.util
import json, glob, os, re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CG = os.path.join(ROOT, "site", "content", "genes")
GENOME_FA = os.path.join(ROOT, "..", "investigate_phylo", "resources", "NC_000962.3.fasta")
CDS_ORIG = os.path.join(ROOT, "..", "investigate_phylo", "resources", "NC_000962.3_CDS.fasta")
AUG = os.path.join(ROOT, "résultats", "phase86c_eggnog", "cds_augmented_new68.fasta")

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def load_genome():
    return "".join(l.strip() for l in open(GENOME_FA) if not l.startswith(">"))


def main():
    genome = load_genome()
    # collect the 68 P16.5d genes with H37Rv coords
    new = []
    for f in glob.glob(os.path.join(CG, "Rv*.json")):
        d = json.load(open(f))
        if not str(d.get("curation_note", "")).startswith("Added P16.5d"):
            continue
        s, e, st = d.get("start_mtbc0"), d.get("end_mtbc0"), d.get("strand")
        if s and e:
            new.append((d["rv"], int(s), int(e), st or "+"))

    # build augmented CDS fasta = original + 68 (nucleotide CDS from the genome)
    with open(AUG, "w") as out:
        out.write(open(CDS_ORIG).read())
        for rv, s, e, st in new:
            sub = genome[s - 1:e]
            if st == "-":
                sub = sub.translate(COMP)[::-1]
            loc = f"complement({s}..{e})" if st == "-" else f"{s}..{e}"
            out.write(f">lcl|NC_000962.3_cds_{rv}_p165d [gene={rv}] [locus_tag={rv}] [location={loc}] [gbkey=CDS]\n{sub}\n")

    # load phase2f as a module, point it at the augmented fasta, run (reuses pass-1 cache)
    spec = importlib.util.spec_from_file_location(
        "p2f", os.path.join(ROOT, "analyses", "phase2f_conservation.py"))
    p2f = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(p2f)
    from pathlib import Path
    if not (p2f.OUTDIR / "snp_counts.tsv").exists():
        raise SystemExit("cache pass1 absent — lancer phase2f complet d'abord")
    p2f.CDS_FA = Path(AUG)
    print(f"phase86e: {len(new)} gènes ajoutés au CDS fasta augmenté ; phase2f pass2 (cache réutilisé)…")
    p2f.main()

    # report on the 68
    got = pn = pseudo = 0
    for rv, *_ in new:
        d = json.load(open(os.path.join(CG, f"{rv}.json")))
        c = d.get("conservation")
        if c and c.get("snp_sites") is not None:
            got += 1
            if c.get("selection", "").startswith("strong"):
                pn += 1
            if c.get("pseudogene_flag"):
                pseudo += 1
    print(f"phase86e: conservation calculée pour {got}/{len(new)} gènes")
    print(f"  sous sélection purifiante forte (vrais gènes): {pn} | candidats pseudogène: {pseudo}")


if __name__ == "__main__":
    main()
