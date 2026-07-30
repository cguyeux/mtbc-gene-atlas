#!/usr/bin/env python3
"""phase56_outgroup_conservation.py -- P10 : couche `outgroup` (conservation M. canettii).

Ingère la table dN/dS par gène M. canettii vs MTBC (produite par
`../Canettii/export_canettii_per_gene_dnds.py`, réutilisant la machinerie Nei-Gojobori de
`analyse_dnds_categories.py` sur le consensus canettii ≥50% des 153 génomes). Apporte à l'atlas
la seule couche de sélection à DIVERGENCE PROFONDE (tout le reste est intra-MTBC quasi-clonal) +
un flag de substitution disruptive fixée dans le clade canettii (pseudogénisation de lignée).

GARDE-FOU (anti-survente) : porter le flag `power` (ok si ≥8 substitutions canettii, sinon low).
Les gènes courts/très conservés ont peu de substitutions → dN/dS non fiable, NE PAS sur-interpréter.
Écrit EN PLACE la couche `outgroup` sur chaque fiche présente dans la table.
Run: python analyses/phase56_outgroup_conservation.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
TSV = ROOT.parent / "Canettii" / "canettii_vs_mtbc" / "canettii_per_gene_dnds.tsv"


def interpret(n, disr, dnds, power):
    bits = []
    if power == "ok":
        if isinstance(dnds, float) and dnds < 0.5:
            bits.append(f"under purifying selection vs M. canettii (deep divergence; dN/dS={dnds}) — a real, "
                        "constrained gene predating the MTBC clonal expansion")
        elif isinstance(dnds, float):
            bits.append(f"elevated dN/dS vs M. canettii ({dnds}) — relaxed or positive selection at deep divergence")
        else:
            bits.append("substitutions present but dS=0 (dN/dS undefined)")
    else:
        bits.append(f"low power ({n} canettii-consensus substitution(s)); present in M. canettii but dN/dS not reliable")
    if disr:
        bits.append(f"carries {disr} M. canettii-clade-fixed disruptive substitution(s) (candidate lineage-specific "
                    "pseudogenisation — cross-check the intra-MTBC pseudogene layer)")
    return "; ".join(bits)


def main():
    tbl = {}
    for l in open(TSV):
        if l.startswith("rv"):
            continue
        p = l.rstrip("\n").split("\t")
        dnds = p[7]
        try:
            dnds = float(dnds)
        except ValueError:
            dnds = None if dnds in ("None", "") else dnds  # keep "Inf"
        tbl[p[0]] = {"n_substitutions": int(p[1]), "n_syn": int(p[2]), "n_nonsyn": int(p[3]),
                     "n_disruptive": int(p[4]), "dN_dS": dnds, "power": p[8]}
    n = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        rec = tbl.get(d["rv"])
        if not rec:
            continue
        rec = dict(rec)
        rec["interpretation"] = interpret(rec["n_substitutions"], rec["n_disruptive"], rec["dN_dS"], rec["power"])
        rec["source"] = "M. canettii (immediate outgroup, 153 genomes) vs H37Rv, Nei-Gojobori dN/dS over the "
        rec["source"] += ">=50% canettii-consensus substitutions (annotation_mtbc P10)"
        d["outgroup"] = {"canettii": rec}
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
    powered = sum(1 for r in tbl.values() if r["power"] == "ok")
    print(f"couche `outgroup` écrite sur {n} fiches ({powered} avec puissance ok, "
          f"{sum(1 for r in tbl.values() if r['n_disruptive'])} avec substitution disruptive canettii).")


if __name__ == "__main__":
    main()
