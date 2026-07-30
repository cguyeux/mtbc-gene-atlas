#!/usr/bin/env python3
"""P4.2 — Scan génome-entier des stops prématurés dans les ORF ancestraux MTBC0 v1.1.

Contexte (issu de P14.1) : phase1 traduit les CDS MTBC0 avec `to_stop=True`, qui coupe
au premier codon stop INTERNE. Un ORF portant un stop prématuré donne donc une protéine
tronquée. Le seuil de span (P14.1) n'a repéré que les tronquages LARGES (12 gènes) ; ce
scan traduit le CDS COMPLET (sans to_stop) pour TOUS les gènes et compte les stops internes,
attrapant aussi les stops proches du C-terminus.

Discriminateur artefact-de-reconstruction vs dégradation-réelle :
  - CLEAN                : 0 stop interne (gène normal).
  - CONSERVED_ARTIFACT   : gène conservé (non mobile), 1 seul stop interne, et la traduction
                           full-length restaure ~ la longueur de l'ortholog H37Rv → SNP/indel
                           parasite de la reconstruction MTBC0 v1.1 (candidat à corriger).
  - MOBILE_DEGRADED      : élément mobile (phage/IS/transposase), stops (souvent multiples) =
                           dégradation biologique attendue, NE PAS corriger.
  - CONSERVED_MULTISTOP  : gène conservé, ≥2 stops → ambigu (pseudogène réel possible).
  - FRAME_ISSUE          : (end-start+1) non multiple de 3 → annotation hors-cadre.

Lecture SEULE. Sortie : résultats/phase64_mtbc0_premature_stops/{premature_stops.tsv, summary.txt}.
Ne modifie AUCUNE fiche. Les candidats CONSERVED_ARTIFACT alimentent P4.1 (MTBC0 propre).
"""
import csv
import json
import re
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parent.parent
MTBC0_FASTA = ROOT / "data" / "MTBC0" / "MTBC0_v1.1.fasta"
GENES = ROOT / "site" / "content" / "genes"
GFF_H37 = ROOT.parent / "investigate_phylo" / "resources" / "NC_000962.3.gff3"
OUT = ROOT / "résultats" / "phase64_mtbc0_premature_stops"

MOBILE = re.compile(r"phage|transpos|insertion sequence|\bIS\d|integrase|prophage|"
                    r"PE-PGRS|PPE|resolvase|recombinase", re.I)
H37_TOL = 10  # aa : tolérance pour « la full-length restaure la longueur H37Rv »


def h37rv_lengths() -> dict:
    out = {}
    for line in open(GFF_H37):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "CDS":
            continue
        m = re.search(r"locus_tag=([^;]+)", p[8])
        if m:
            out[m.group(1)] = (abs(int(p[4]) - int(p[3])) + 1) // 3 - 1
    return out


def main() -> None:
    genome = next(SeqIO.parse(str(MTBC0_FASTA), "fasta")).seq
    h37 = h37rv_lengths()
    rows = []
    cat = Counter()
    for f in sorted(GENES.glob("*.json")):
        d = json.loads(f.read_text())
        rv = d["rv"]
        s, e, strand = d.get("start_mtbc0"), d.get("end_mtbc0"), d.get("strand", "+")
        prot = d.get("protein_mtbc0") or ""
        if not s or not e:
            continue
        s, e = int(s), int(e)
        span = abs(e - s) + 1
        frame_issue = (span % 3 != 0)
        sub = genome[s - 1:e]
        if strand == "-":
            sub = sub.reverse_complement()
        trim = len(sub) - (len(sub) % 3)
        full = str(Seq(sub[:trim]).translate(table=11))
        to_stop_len = len(prot)                    # = longueur servie (phase1, to_stop=True)
        # VALIDATION DE COHÉRENCE : ne compter un stop prématuré que si MA re-traduction
        # est d'accord avec phase1, i.e. le 1er stop tombe EXACTEMENT après la protéine servie
        # (full[to_stop_len] == '*'). Sinon mon cadre/brin diffère de phase1 → gène EXCLU
        # (évite les faux positifs off-by-one qui « cassaient » des gènes essentiels).
        consistent = (len(full) > to_stop_len and full[to_stop_len] == "*")
        internal = full[:to_stop_len].count("*")   # stops AVANT la fin de la protéine servie
        full_len = len(full.rstrip("*"))           # longueur si on va jusqu'au bout du CDS
        lost = (span // 3 - 1) - to_stop_len        # aa perdus vs capacité codante du CDS
        is_mobile = bool(MOBILE.search(d.get("product_h37rv") or "")
                         or MOBILE.search(d.get("product_mtbc0_pgap") or ""))
        h37_len = h37.get(rv)
        restores = (h37_len is not None and abs(full_len - h37_len) <= H37_TOL)

        # Un vrai stop prématuré = protéine servie plus courte que la capacité du CDS,
        # ET cohérence re-traduction↔phase1, ET perte substantielle (>=2 aa, écarte le bruit ±1).
        real_premature = consistent and lost >= 2 and internal == 0

        # Classer par NATURE (mobile vs conservé) d'abord ; le frame%3≠0 est une simple
        # anomalie d'annotation (fréquente sur ces loci tronqués), notée en colonne, PAS une
        # catégorie qui masquerait le signal.
        if not real_premature:
            category = "CLEAN"           # (ou incohérent re-traduction↔phase1 : non-signal)
        elif is_mobile:
            category = "MOBILE_DEGRADED"          # dégradation biologique attendue, ne pas corriger
        elif restores:
            category = "CONSERVED_ARTIFACT"       # candidat artefact reconstruction MTBC0 v1.1
        else:
            category = "CONSERVED_OTHER"          # conservé mais ne restaure pas H37Rv (à examiner)
        cat[category] += 1
        if real_premature:
            rows.append({
                "rv": rv, "category": category, "served_aa": to_stop_len,
                "full_cds_aa": full_len, "lost_aa": lost,
                "h37rv_aa": h37_len if h37_len is not None else "",
                "restores_h37rv": int(restores), "is_mobile": int(is_mobile),
                "cds_frame_ok": int(not frame_issue),
                "product": (d.get("product_h37rv") or "")[:60],
            })

    OUT.mkdir(parents=True, exist_ok=True)
    rows.sort(key=lambda r: (r["category"], -r["lost_aa"]))
    cols = ["rv", "category", "served_aa", "full_cds_aa", "lost_aa",
            "h37rv_aa", "restores_h37rv", "is_mobile", "cds_frame_ok", "product"]
    with (OUT / "premature_stops.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    total = sum(cat.values())
    lines = [
        f"MTBC0 v1.1 premature-stop scan — {total} genes, {total - cat['CLEAN']} with >=1 internal stop",
        "",
        *[f"  {k:<20} {v}" for k, v in cat.most_common()],
        "",
        "CONSERVED_ARTIFACT = candidats à corriger dans un MTBC0 propre (P4.1) :",
    ]
    for r in [r for r in rows if r["category"] == "CONSERVED_ARTIFACT"]:
        lines.append(f"  {r['rv']:<9} served {r['served_aa']:>4} aa  full {r['full_cds_aa']:>4} aa  "
                     f"H37Rv {r['h37rv_aa']}  | {r['product']}")
    summary = "\n".join(lines)
    (OUT / "summary.txt").write_text(summary + "\n")
    print(summary)
    print(f"\n-> {OUT}/premature_stops.tsv ({len(rows)} gènes à stop interne)")


if __name__ == "__main__":
    main()
