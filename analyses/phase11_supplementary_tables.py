#!/usr/bin/env python3
"""phase11_supplementary_tables.py -- matérialise les tables supplémentaires S1-S4.

Ferme la remarque du relecteur « supplementary vide » : produit les livrables FAIR
attendus d'un article-ressource, dérivés de la source de vérité site/content/genes/*.json.

  S1  table_s1_atlas.tsv                  atlas complet (3906 gènes, toutes couches)
  S2  table_s2_ahead_of_uniprot.tsv       gènes 'ahead of UniProt' (handle eggNOG, UniProt sans fonction)
  S3  table_s3_pseudogene_candidates.tsv  candidats pseudogène (charge de disruption clonale 1-50 %)
  S4  table_s4_ta_repertoire.tsv          répertoire toxine-antitoxine (cas d'étude), depuis data/ta_case_study/
  (S5 table_s5_ec_discordances.tsv est produite par phase12_concordance.py)

Écrit aussi deux fragments longtable LaTeX (s3_longtable.tex, s4_longtable.tex) à \\input
dans supplementary.tex, pour que les chiffres typographiés restent synchronisés avec les données.

Sortie : article/supplementary_materials/. Stdlib pure.
Run: python analyses/phase11_supplementary_tables.py
"""
from __future__ import annotations
import csv, glob, json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
TA = ROOT / "data" / "ta_case_study" / "ta_candidates.tsv"
OUT = ROOT / "article" / "supplementary_materials"
OUT.mkdir(parents=True, exist_ok=True)


def ne(x):
    return bool(x) and str(x).strip() not in ("", "[]", "{}")


def j(x, sep=";"):
    if isinstance(x, (list, tuple)):
        return sep.join(str(e) for e in x if e and str(e) != "-")
    return str(x or "")


def is_hyp(d):
    return "hypothetical" in (d.get("product_h37rv") or "").lower()


def texesc(s):
    s = str(s or "")
    for a, b in [("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("_", r"\_"),
                 ("#", r"\#"), ("$", r"\$"), ("{", r"\{"), ("}", r"\}"), ("~", r"\textasciitilde{}"),
                 ("^", r"\textasciicircum{}")]:
        s = s.replace(a, b)
    return s


def trunc(s, n):
    s = str(s or "")
    return s if len(s) <= n else s[: n - 1].rstrip() + "…"


def load():
    return [json.load(open(f)) for f in sorted(glob.glob(str(GENES / "*.json")))]


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    print(f"  {path.name}: {len(rows)} rows")


def main():
    g = load()
    print(f"Loaded {len(g)} fiches")

    # ── S1 : atlas complet ────────────────────────────────────────────────────
    h1 = ["locus", "gene", "len_aa", "product_h37rv", "function_revised", "verdict",
          "confidence", "cog_cat", "ec_eggnog", "kegg_ko", "uniprot_acc", "uniprot_reviewed",
          "pN_pS", "selection", "pseudogene_flag", "string_anchor"]
    r1 = []
    for d in g:
        e = d.get("eggnog") or {}; u = d.get("uniprot") or {}
        c = d.get("conservation") or {}; s = (d.get("string") or {}).get("anchor") or {}
        r1.append([d["rv"], d.get("gene") or "", d.get("len_aa") or "", d.get("product_h37rv") or "",
                   d.get("function_revised") or "", d.get("verdict") or "", d.get("confidence") or "",
                   e.get("cog_cat") or "", j(e.get("ec")), j(e.get("kegg_ko")),
                   u.get("acc") or "", u.get("reviewed") if u else "",
                   c.get("pN_pS", ""), c.get("selection") or "", c.get("pseudogene_flag", ""),
                   s.get("gene") or s.get("rv") or ""])
    write_tsv(OUT / "table_s1_atlas.tsv", h1, r1)

    # ── S2 : ahead of UniProt ─────────────────────────────────────────────────
    def egg_handle(d):
        e = d.get("eggnog") or {}
        return any(ne(e.get(k)) for k in ("cog_cat", "ec", "kegg_ko"))
    def up_unchar(d):
        return not ne((d.get("uniprot") or {}).get("function"))
    h2 = ["locus", "gene", "product_h37rv", "legacy_hypothetical", "cog_cat", "ec_eggnog", "kegg_ko", "uniprot_acc"]
    r2 = []
    for d in g:
        if egg_handle(d) and up_unchar(d):
            e = d.get("eggnog") or {}
            r2.append([d["rv"], d.get("gene") or "", d.get("product_h37rv") or "", is_hyp(d),
                       e.get("cog_cat") or "", j(e.get("ec")), j(e.get("kegg_ko")),
                       (d.get("uniprot") or {}).get("acc") or ""])
    write_tsv(OUT / "table_s2_ahead_of_uniprot.tsv", h2, r2)

    # ── S3 : candidats pseudogène ─────────────────────────────────────────────
    h3 = ["locus", "gene", "product_h37rv", "legacy_hypothetical", "disrupt_mode",
          "disrupt_frac_pct", "disrupt_strains", "n_strains", "pN_pS"]
    r3 = []
    for d in g:
        c = d.get("conservation") or {}
        if c.get("pseudogene_flag") is True:
            r3.append([d["rv"], d.get("gene") or "", d.get("product_h37rv") or "", is_hyp(d),
                       c.get("disrupt_mode") or "", round(100 * (c.get("disrupt_frac") or 0), 2),
                       c.get("disrupt_strains") or "", c.get("n_strains") or "", c.get("pN_pS", "")])
    r3.sort(key=lambda x: (-float(x[5] or 0)))
    write_tsv(OUT / "table_s3_pseudogene_candidates.tsv", h3, r3)

    # ── S4 : répertoire TA (copie normalisée) ─────────────────────────────────
    ta_rows = list(csv.DictReader(open(TA), delimiter="\t"))
    h4 = ["locus", "len_aa", "refseq_product", "verdict", "function_revised"]
    r4 = [[t["rv"], t["len_aa"], t["product"], t["verdict"], t["function_revised"]] for t in ta_rows]
    write_tsv(OUT / "table_s4_ta_repertoire.tsv", h4, r4)

    # ── fragment longtable S3 ─────────────────────────────────────────────────
    hyp3 = [r for r in r3 if r[3]]  # legacy hypothetical uniquement (focus manuscrit)
    with open(OUT / "s3_longtable.tex", "w") as fh:
        fh.write("% auto-généré par phase11 -- NE PAS éditer à la main\n")
        fh.write("\\begin{longtable}{@{}llp{6.8cm}rr@{}}\n")
        fh.write("\\caption{Pseudogene candidates among the historically hypothetical genes "
                 f"({len(hyp3)} loci), flagged by a clonal nonsense/frameshift disruption load "
                 "(fewer than four distinct disrupting alleles) within the conservative 1--50\\% "
                 "frequency window across 145{,}209 genomes. Disruption frequency (fraction of "
                 "genomes carrying a disrupting allele) and supporting strain count are shown, sorted "
                 "by frequency. The full list of all "
                 f"{len(r3)} pseudogene candidates (all products) is in \\texttt{{table\\_s3\\_pseudogene\\_candidates.tsv}}.}}"
                 "\\label{tab:s3}\\\\\n")
        fh.write("\\toprule\nLocus & Gene & RefSeq product & Freq.\\ (\\%) & Strains \\\\\n\\midrule\n\\endfirsthead\n")
        fh.write("\\toprule\nLocus & Gene & RefSeq product & Freq.\\ (\\%) & Strains \\\\\n\\midrule\n\\endhead\n")
        for rv, gene, prod, _, mode, frac, strains, nstr, pnps in hyp3:
            fh.write(f"{texesc(rv)} & {texesc(gene)} & {texesc(trunc(prod,56))} & "
                     f"{frac} & {strains} \\\\\n")
        fh.write("\\bottomrule\n\\end{longtable}\n")
    print(f"  s3_longtable.tex: {len(hyp3)} hypothetical pseudogene rows")

    # ── fragment longtable S4 ─────────────────────────────────────────────────
    with open(OUT / "s4_longtable.tex", "w") as fh:
        fh.write("% auto-généré par phase11 -- NE PAS éditer à la main\n")
        fh.write("\\begin{longtable}{@{}llp{4.3cm}lp{4.6cm}@{}}\n")
        fh.write(f"\\caption{{Toxin--antitoxin repertoire surfaced by the atlas ({len(r4)} candidate loci). "
                 "Verdict is the atlas outcome (requalified, family\\_assigned, or removed as non-TA); the "
                 "revised function is abbreviated. Full text in \\texttt{table\\_s4\\_ta\\_repertoire.tsv}.}"
                 "\\label{tab:s4}\\\\\n")
        fh.write("\\toprule\nLocus & Len & RefSeq product & Verdict & Revised function \\\\\n\\midrule\n\\endfirsthead\n")
        fh.write("\\toprule\nLocus & Len & RefSeq product & Verdict & Revised function \\\\\n\\midrule\n\\endhead\n")
        for rv, ln, prod, verdict, fn in r4:
            fh.write(f"{texesc(rv)} & {ln} & {texesc(trunc(prod,34))} & {texesc(verdict.replace('_',' '))} & "
                     f"{texesc(trunc(fn,58))} \\\\\n")
        fh.write("\\bottomrule\n\\end{longtable}\n")
    print(f"  s4_longtable.tex: {len(r4)} TA rows")

    # récapitulatif chiffres
    print("\nCounts: S1", len(r1), "| S2", len(r2), "(hyp", sum(1 for r in r2 if r[3]),
          ") | S3", len(r3), "(hyp", len(hyp3), ") | S4", len(r4))


if __name__ == "__main__":
    main()
