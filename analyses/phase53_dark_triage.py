#!/usr/bin/env python3
"""phase53_dark_triage.py -- P8 étape 0 : triage de résolvabilité des gènes dark.

READ-ONLY (n'écrit qu'un TSV de worklist dans résultats/, ne touche aucune fiche). Profile tous les
verdict=="dark", compte chaque signal d'évidence, et partitionne par catégorie de résolvabilité pour
ordonner la reprise par batch (HIGH d'abord). Handles FORTS (dérivent une fonction) :
  - opéron : ≥1 voisin co-transcrit NOMMÉ et de fonction connue (fiche non-dark, gene_name/eggnog name)
  - STRING : partenaire NOMMÉ non-hypothétique à combined_no_tm ≥700
  - iModulon : appartenance à un iModulon dont le regulator est NOMMÉ (contexte de condition)
  - phénotype : mutant_phenotypes significatif OU phenotype_lead
Garde-fou d'orthogonalité : le canal STRING `neighborhood` == synténie == l'opéron (même fait) → non
compté comme 2e handle ; opéron et STRING(no-tm, non-neighborhood-only) comptent séparément.
Sortie : résultats/phase53_triage/dark_worklist.tsv + décomptes.
Run: python analyses/phase53_dark_triage.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "résultats" / "phase53_triage"


def load_all():
    G = {}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        G[d["rv"]] = d
    return G


def named_known(G, rv):
    """le gène rv (voisin) est-il NOMMÉ et de fonction connue (non-dark) ?"""
    d = G.get(rv) or {}
    if d.get("verdict") == "dark":
        return None
    nm = d.get("gene_name") or (d.get("eggnog") or {}).get("preferred_name")
    return nm or None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    G = load_all()
    dark = [d for d in G.values() if d.get("verdict") == "dark"]

    rows = []
    counts = {k: 0 for k in ("operon_named", "string_named_strong", "string_exp", "imod_reg",
                             "phenotype", "vuln", "vuln_hi", "essential", "purifying_strong",
                             "purifying", "ptm", "struct", "proteomics", "micro65", "micro100",
                             "mobile", "pseudogene")}
    loc_counts = {}
    for d in dark:
        rv = d["rv"]
        seq = d.get("protein_mtbc0") or ""
        L = len(seq) or (d.get("len_aa") or 0)

        # --- handles forts ---
        op = (d.get("genomic_context") or {}).get("operon") or []
        op_named = [named_known(G, m.get("locus")) for m in op if m.get("locus") != rv]
        op_named = [x for x in op_named if x]
        h_operon = bool(op_named)

        s = d.get("string") or {}
        parts = s.get("partners") or []
        h_string = any((not p.get("hypothetical")) and p.get("gene") and (p.get("combined_no_tm") or 0) >= 700
                       for p in parts)
        h_string_exp = any((p.get("channels") or {}).get("experimental", 0) >= 400 for p in parts)

        # --- PAIRE ORTHOGONALE (leçon batch 1 P8.1) : opéron nommé + partenaire STRING EXTERNE à l'opéron
        # via un canal non-neighborhood (cooccurence/fusion/coexpression/experimental/database ≥400).
        # Une fusion/cooccurrence ENTRE voisins d'opéron = co-hérédité synténique = même fait, PAS orthogonal.
        NON_NB = ("cooccurence", "cooccurrence", "fusion", "coexpression", "experimental", "database")
        op_loci = {m.get("locus") for m in op}
        ortho_partner = None
        for p in parts:
            if p.get("hypothetical") or not p.get("gene") or p.get("rv") in op_loci or p.get("rv") == rv:
                continue
            ch = p.get("channels") or {}
            if max((ch.get(k, 0) for k in NON_NB), default=0) >= 400:
                ortho_partner = ortho_partner or f"{p['gene']}({max((k for k in NON_NB if ch.get(k,0)>=400), key=lambda k: ch.get(k,0))})"
        ortho_pair = bool(op_named) and bool(ortho_partner)

        exprs = (d.get("expression") or {}).get("imodulons") or []
        h_imod = any(m.get("regulator") for m in exprs)

        mp = d.get("mutant_phenotypes") or {}
        h_pheno = bool(mp.get("n_significant")) or bool(d.get("phenotype_lead"))

        strong = sum([h_operon, h_string, h_imod, h_pheno])

        # --- signaux secondaires / classifiants ---
        cons = d.get("conservation") or {}
        sel = (cons.get("selection") or "")
        is_pseudo = bool(cons.get("pseudogene_flag") or cons.get("pseudogene_candidate"))
        rd = d.get("rd") or {}
        is_mobile = bool(rd.get("rds") or rd.get("in_rd") or (isinstance(rd, dict) and rd.get("deleted")))
        vi = (d.get("vulnerability") or {}).get("vi")
        loc = (d.get("localization") or {}).get("type") or "?"
        pdb = d.get("pdb") or {}; plaf = d.get("plddt_af") or {}; sa = d.get("struct_af") or {}
        has_struct = bool((isinstance(pdb, dict) and pdb.get("n_structures")) or
                          (isinstance(plaf, dict) and (plaf.get("mean_plddt") or 0) >= 70) or
                          (isinstance(sa, dict) and sa.get("hits")) or d.get("pfam_tentative"))
        has_ptm = bool((d.get("ptm") or {}).get("n_sites"))
        ms = bool((d.get("proteomics") or {}).get("detected"))
        secondary = any([vi is not None and vi <= -6, has_struct, has_ptm, bool(d.get("regulation")),
                         (d.get("expression") or {}).get("n_imodulons")])

        # --- catégorie (précédence : pseudo > mobile > HIGH > MEDIUM > micro > weak > orphan) ---
        if is_pseudo:
            cat = "PSEUDOGENE"
        elif is_mobile:
            cat = "MOBILE"
        elif strong >= 2:
            cat = "HIGH"
        elif strong == 1:
            cat = "MEDIUM"
        elif L and L < 65:
            cat = "MICRO-ORF"
        elif secondary:
            cat = "WEAK"
        else:
            cat = "ORPHAN"

        # --- décomptes ---
        counts["operon_named"] += h_operon
        counts["string_named_strong"] += h_string
        counts["string_exp"] += h_string_exp
        counts["imod_reg"] += h_imod
        counts["phenotype"] += h_pheno
        counts["vuln"] += vi is not None
        counts["vuln_hi"] += (vi is not None and vi <= -6)
        counts["essential"] += bool((d.get("essentiality") or {}).get("essential"))
        counts["purifying_strong"] += sel.startswith("strong purifying")
        counts["purifying"] += ("purifying" in sel)
        counts["ptm"] += has_ptm
        counts["struct"] += has_struct
        counts["proteomics"] += ms
        counts["micro65"] += bool(L and L < 65)
        counts["micro100"] += bool(L and L < 100)
        counts["mobile"] += is_mobile
        counts["pseudogene"] += is_pseudo
        loc_counts[loc] = loc_counts.get(loc, 0) + 1

        handles = ";".join(h for h, ok in [("operon", h_operon), ("string", h_string),
                                           ("imod", h_imod), ("pheno", h_pheno)] if ok)
        sig = []
        if op_named:
            sig.append(f"op:{','.join(op_named[:2])}")
        if h_string:
            best = max((p for p in parts if p.get("gene") and not p.get("hypothetical")),
                       key=lambda p: p.get("combined_no_tm", 0), default=None)
            if best:
                sig.append(f"str:{best['gene']}({best.get('combined_no_tm')})")
        if h_imod:
            sig.append("iM:" + ",".join(m["regulator"] for m in exprs if m.get("regulator"))[:24])
        if vi is not None:
            sig.append(f"VI={vi:.1f}")
        rows.append({"rv": rv, "cat": cat, "strong": strong, "len": L, "sel": sel[:16],
                     "loc": loc, "handles": handles, "signals": " | ".join(sig),
                     "ortho": ortho_pair, "ortho_partner": ortho_partner or ""})

    # tri : PAIRE ORTHOGONALE d'abord (vrais candidats-flip), puis HIGH, par nb de handles forts décroissant
    order = {"HIGH": 0, "MEDIUM": 1, "WEAK": 2, "MICRO-ORF": 3, "MOBILE": 4, "PSEUDOGENE": 5, "ORPHAN": 6}
    rows.sort(key=lambda r: (not r["ortho"], order[r["cat"]], -r["strong"], r["rv"]))

    with open(OUT / "dark_worklist.tsv", "w") as fh:
        fh.write("rv\tcategory\tortho_pair\tn_strong_handles\tlen_aa\tselection\tlocalization\thandles\tortho_partner\tsignals\n")
        for r in rows:
            fh.write(f"{r['rv']}\t{r['cat']}\t{int(r['ortho'])}\t{r['strong']}\t{r['len']}\t{r['sel']}\t{r['loc']}\t"
                     f"{r['handles']}\t{r['ortho_partner']}\t{r['signals']}\n")

    from collections import Counter
    cat_counts = Counter(r["cat"] for r in rows)
    print(f"=== TRIAGE DES {len(dark)} DARK ===")
    print("Catégories :", dict(cat_counts), "| somme =", sum(cat_counts.values()))
    print("\nSignaux (parmi les dark) :")
    for k, v in counts.items():
        print(f"  {k:<20} {v}")
    print("\nLocalisation :", loc_counts)
    ortho = [r for r in rows if r["ortho"]]
    print(f"\n=== CANDIDATS-FLIP prioritaires : {len(ortho)} PAIRES ORTHOGONALES "
          f"(opéron nommé + partenaire STRING EXTERNE non-neighborhood) ===")
    for r in ortho:
        print(f"  {r['rv']:<9} [{r['cat']}] len={r['len']} {r['loc']} :: op:{r['handles']} + ext-STRING:{r['ortho_partner']} :: {r['signals']}")
    print(f"\nWorklist (ortho_pair d'abord) : {OUT/'dark_worklist.tsv'}")


if __name__ == "__main__":
    main()
