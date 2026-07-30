#!/usr/bin/env python3
"""phase60_harden_mtbc_specific.py -- P11.2 : durcir la shortlist des dark spécifiques-MTBC.

Corrige les faiblesses trouvées par l'évaluation adversariale (Workflow P11.1) avant toute valorisation :
  (a) TROU MOBILE : l'exclusion des IS/phage par regex sur le NOM rate les dark (sans nom). Remplacé par un
      test agnostique au nom = MULTIPLICITÉ dans le génome (tblastn du candidat vs H37Rv ; ≥2 copies fortes =
      élément répété/mobile → exclu).
  (b) exclure les PSEUDOGÈNES (pseudogene_flag, ou pN/pS diversifiant ≥1.5).
  (c) exclure les DÉJÀ-CARACTÉRISÉS repérés par l'évaluation (Rv1048c, Rv3126c, Rv0378) — liste à compléter
      par une passe littérature/HHpred (P11.2c, hand-off).
Entrée : phase59 confirmed-insertion + dark. Sortie : résultats/phase60/dark_mtbc_specific_hardened.tsv.
Run: python analyses/phase60_harden_mtbc_specific.py
"""
from __future__ import annotations
import json, subprocess, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
CONF = ROOT / "résultats" / "phase59_synteny" / "mtbc_specific_confirmed.tsv"
H37RV = ROOT.parent / "investigate_phylo" / "resources" / "NC_000962.3.fasta"
OUT = ROOT / "résultats" / "phase60"
CHARACTERISED = {"Rv1048c", "Rv3126c", "Rv0378"}   # repérés par l'évaluation P11.1 (à compléter)


def self_multiplicity(cands: dict[str, str], tmp: Path) -> dict[str, int]:
    """nb de loci du génome H37Rv touchés par chaque protéine candidate (>1 = répété/mobile)."""
    q = tmp / "cands.faa"
    q.write_text("".join(f">{rv}\n{seq}\n" for rv, seq in cands.items()))
    db = tmp / "h37rv"
    subprocess.run(["makeblastdb", "-in", str(H37RV), "-dbtype", "nucl", "-out", str(db)],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out = tmp / "self.tsv"
    subprocess.run(["tblastn", "-query", str(q), "-db", str(db), "-out", str(out), "-evalue", "1e-5",
                    "-num_threads", "8", "-outfmt", "6 qseqid sstart send pident qcovs"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # compter les loci DISTINCTS (fenêtres génomiques non chevauchantes) à hit fort par requête
    hits: dict[str, list] = {}
    for line in open(out):
        p = line.rstrip("\n").split("\t")
        if len(p) < 5:
            continue
        rv = p[0]; s = min(int(p[1]), int(p[2])); pid = float(p[3]); qc = float(p[4])
        if pid >= 40 and qc >= 40:
            hits.setdefault(rv, []).append(s)
    n_loci = {}
    for rv, starts in hits.items():
        starts.sort()
        loci = 1
        for i in range(1, len(starts)):
            if starts[i] - starts[i - 1] > 5000:   # >5 kb = locus distinct
                loci += 1
        n_loci[rv] = loci
    return n_loci


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    conf = []
    for l in open(CONF):
        if l.startswith("rv"):
            continue
        p = l.rstrip("\n").split("\t")
        if p[3] == "1" and p[2] == "confirmed-insertion":   # dark + confirmed-insertion
            conf.append(p[0])
    cands = {}
    meta = {}
    for rv in conf:
        d = json.load(open(GENES / f"{rv}.json"))
        cands[rv] = d.get("protein_mtbc0") or ""
        cons = d.get("conservation") or {}
        meta[rv] = {"len": len(cands[rv]), "pnps": cons.get("pN_pS"),
                    "pseudo": bool(cons.get("pseudogene_flag")),
                    "ess": bool((d.get("essentiality") or {}).get("essential")),
                    "vi": (d.get("vulnerability") or {}).get("vi"),
                    "vi_lo": (d.get("vulnerability") or {}).get("ci_low"),
                    "vi_hi": (d.get("vulnerability") or {}).get("ci_high")}
    with tempfile.TemporaryDirectory() as td:
        nloci = self_multiplicity(cands, Path(td))

    rows = []
    for rv in conf:
        m = meta[rv]
        loci = nloci.get(rv, 1)
        reasons = []
        if m["pseudo"] or (m["pnps"] is not None and m["pnps"] >= 1.5):
            reasons.append("pseudogene/diversifying")
        if loci > 1:
            reasons.append(f"repeated x{loci} (mobile?)")
        if rv in CHARACTERISED:
            reasons.append("already characterised")
        # vulnérabilité robuste = IC entièrement < 0
        robust_vuln = (m["vi_hi"] is not None and m["vi_hi"] < 0)
        rows.append({"rv": rv, "len": m["len"], "ess": m["ess"], "vi": m["vi"],
                     "robust_vuln": robust_vuln, "n_loci": loci, "pnps": m["pnps"],
                     "excluded": "; ".join(reasons), "keep": not reasons})
    kept = [r for r in rows if r["keep"]]
    kept.sort(key=lambda r: (not (r["ess"] or r["robust_vuln"]), -r["len"]))
    rows.sort(key=lambda r: (not r["keep"], -r["len"]))
    with open(OUT / "dark_mtbc_specific_hardened.tsv", "w") as fh:
        fh.write("rv\tkeep\tlen\tessential\trobust_vuln\tvi\tn_genome_loci\tpN_pS\texcluded_reason\n")
        for r in rows:
            fh.write(f"{r['rv']}\t{int(r['keep'])}\t{r['len']}\t{int(r['ess'])}\t{int(r['robust_vuln'])}\t"
                     f"{r['vi'] if r['vi'] is not None else ''}\t{r['n_loci']}\t{r['pnps'] if r['pnps'] is not None else ''}\t{r['excluded']}\n")
    print(f"26 dark confirmed-insertion → {len(kept)} RETENUS après durcissement, {len(rows)-len(kept)} exclus.")
    print("Exclusions :")
    for r in [x for x in rows if not x["keep"]]:
        print(f"  {r['rv']}: {r['excluded']}")
    print("\nShortlist DURCIE (priorité essentiel/vulnérable-robuste) :")
    for r in kept[:12]:
        tag = "ESS" if r["ess"] else ("VULN" if r["robust_vuln"] else "")
        print(f"  {r['rv']} ({r['len']}aa) {tag} VI={r['vi']}")
    print(f"\nRapport : {OUT/'dark_mtbc_specific_hardened.tsv'}")


if __name__ == "__main__":
    main()
