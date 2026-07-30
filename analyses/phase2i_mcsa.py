#!/usr/bin/env python3
"""phase2i_mcsa.py -- couche M-CSA : « même repli » -> « enzyme active ? » sur le protéome.

Pour chaque fiche gène ayant un hit Foldseek significatif (phase2c, champ struct_hits) dont la cible
PDB est une enzyme cataloguée M-CSA, aligne la protéine requête sur la séquence de référence M-CSA et
vérifie la conservation des RÉSIDUS CATALYTIQUES (skill active-site-check). Écrit un champ `mcsa` sur la
fiche (mode --write) : {mcsa_id, ec, n_catalytic, n_present, n_identical, verdict, residues}.

Réutilise le client M-CSA testé du skill active-site-check ; Biopython pour l'alignement.
Séquences requête : H37Rv CDS fasta (header [locus_tag=...]).

NOTE DE PORTÉE (mesurée le 2026-06-09) : seules 86/3906 fiches ont un hit Foldseek significatif (Foldseek
n'a tourné que sur ~350 modèles dark), et l'intersection avec les ~1002 PDB de RÉFÉRENCE M-CSA = 1 gène
(Rv3205c, déjà curé comme faux positif). Cette couche n'a donc d'intérêt proteome-wide que si (a) la
couverture Foldseek est étendue (AlphaFold DB pour tout le protéome), et/ou (b) on bascule sur le route
HOMOLOGUES de M-CSA (homologues_residues.json, ~15541 PDB, >200 Mo). En l'état : tool prêt, gain ~nul.

Usage :  python analyses/phase2i_mcsa.py            # rapport seul
         python analyses/phase2i_mcsa.py --write    # + enrichit les fiches matchées
"""
from __future__ import annotations
import argparse, glob, json, os, sys, urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ATLAS = HERE.parent
SKILL = Path.home() / "docs/codes/claude_plugins/bio_pathogens/skills/active-site-check"
CDS = ATLAS.parent / "investigate_phylo/resources/NC_000962.3_CDS.fasta"
GENES = ATLAS / "site/content/genes"
sys.path.insert(0, str(SKILL))
import active_site_check as asc  # noqa: E402

try:
    from Bio.Align import PairwiseAligner, substitution_matrices
except ImportError:
    sys.exit("Biopython requis (pip install biopython).")

UNIPROT_CACHE = Path(os.environ.get("MCSA_CACHE", Path.home() / ".cache" / "mcsa")) / "uniprot"


def pdb_of(target: str) -> str:
    return target.split("-")[0].lower()


def load_cds_seqs() -> dict:
    seqs, tag, buf = {}, None, []
    for line in open(CDS):
        if line.startswith(">"):
            if tag:
                seqs[tag] = "".join(buf)
            buf = []
            tag = None
            for tok in line.split("["):
                if tok.startswith("locus_tag="):
                    tag = tok[len("locus_tag="):].rstrip("] \n")
        else:
            buf.append(line.strip())
    if tag:
        seqs[tag] = "".join(buf)
    return seqs


def uniprot_seq(acc: str) -> str:
    UNIPROT_CACHE.mkdir(parents=True, exist_ok=True)
    f = UNIPROT_CACHE / f"{acc}.fasta"
    if not f.exists():
        data = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=60).read().decode()
        f.write_text(data)
    return "".join(l for l in f.read_text().splitlines() if not l.startswith(">"))


def aligner():
    a = PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score, a.extend_gap_score = -11, -1
    return a


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true", help="écrire le champ mcsa dans les fiches matchées")
    a = ap.parse_args(argv)

    print("Chargement M-CSA (cache)...")
    entries = asc.load_all_entries()
    pdb2entry = {}
    for e in entries:
        for r in asc.catalytic_residues(e):
            if r.get("pdb_id"):
                pdb2entry.setdefault(r["pdb_id"].lower(), e)
    print(f"  PDB de référence M-CSA : {len(pdb2entry)}")

    seqs = load_cds_seqs()
    al = aligner()
    n_fiche = n_match = n_active = n_fold = 0
    rows_out = []
    for f in sorted(glob.glob(str(GENES / "*.json"))):
        d = json.load(open(f))
        hits = [h for h in (d.get("struct_hits") or []) if h.get("significant")]
        # inclure les hits AlphaFold génome-entier (P3.2), seulement si le modèle est confiant (pLDDT>=70)
        af = d.get("struct_af") or {}
        if af.get("plddt_gated"):
            hits = hits + [h for h in af.get("hits", []) if h.get("significant")]
        if not hits:
            continue
        n_fiche += 1
        # premier hit dont la cible PDB est une enzyme M-CSA de référence
        match = next((pdb2entry[pdb_of(h["target"])] for h in hits if pdb_of(h["target"]) in pdb2entry), None)
        if not match:
            continue
        rv = d["rv"]
        qseq = d.get("protein_mtbc0") or d.get("protein_mtbc") or seqs.get(rv) or ""
        if not qseq:
            print(f"  {rv}: pas de séquence requête -> skip")
            continue
        cat = asc.catalytic_residues(match)
        # certaines entrées M-CSA listent plusieurs accessions ("Q06128, Q06129") -> prendre la première
        up = (match.get("reference_uniprot_id") or "").split(",")[0].strip()
        try:
            ref = uniprot_seq(up)
        except Exception as ex:  # noqa: BLE001
            print(f"  {rv}: UniProt {up} indisponible ({ex}) -> skip")
            continue
        # sanitiser : BLOSUM62 rejette les résidus hors alphabet (U sélénocystéine, O, *, ...)
        clean = lambda s: "".join(c if c in "ACDEFGHIKLMNPQRSTVWYBZX*" else "X" for c in s.upper())
        try:
            aln = al.align(clean(qseq), clean(ref))[0]
        except Exception as ex:  # noqa: BLE001
            print(f"  {rv}: alignement échoué ({ex}) -> skip")
            continue
        rows, n, present, ident = asc.map_active_site(str(aln[0]), str(aln[1]), 1, cat)
        verdict = asc._verdict(n, present, ident)
        n_match += 1
        n_active += verdict.startswith("ACTIVE")
        n_fold += verdict.startswith("FOLD")
        rows_out.append((rv, match.get("mcsa_id"), asc._ec_of(match), n, present, ident, verdict))
        if a.write:
            d["mcsa"] = {"mcsa_id": match.get("mcsa_id"), "ec": asc._ec_of(match),
                         "n_catalytic": n, "n_present": present, "n_identical": ident,
                         "verdict": verdict, "residues": rows}
            Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False) + "\n")

    print(f"\nFiches à hit significatif : {n_fiche} | matchées M-CSA (réf PDB) : {n_match} "
          f"| ACTIVE : {n_active} | FOLD-ONLY : {n_fold}")
    for rv, mid, ec, n, present, ident, v in rows_out:
        print(f"  {rv}  M-CSA {mid} EC {ec}  cat {ident}/{n} id ({present}/{n} aln)  -> {v}")
    if a.write:
        print("\n(--write) Champs mcsa écrits. Relancer l'ingestion : cd site && python -m backend.ingest")
    else:
        print("\n(rapport seul ; --write pour enrichir les fiches)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
