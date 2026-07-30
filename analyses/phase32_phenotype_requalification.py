#!/usr/bin/env python3
"""phase32_phenotype_requalification.py -- requalification des hypothétiques par PHÉNOTYPE (P5.2b).

Née de P5.2 : 272 « hypothetical » ont un phénotype conditionnel (MtbTnDB). Ici on ISOLE et
PRIORISE ceux à signal fort (défaut in vivo net, sensibilisation à un antibiotique, stress),
et on SYNTHÉTISE une hypothèse de fonction citable en INTÉGRANT les couches déjà en ligne
(conservation, essentialité, STRING, régulon, opéron, localisation, structure).

GARDE-FOU (leçons KB) : le lead est une HYPOTHÈSE de contexte fonctionnel dérivée du phénotype
+ corroboration, PAS une fonction prouvée. Cadré « lead / à valider ». Pour la conservation,
un gène quasi-invariant (peu de sites segrégeants sur 145k souches) EST contraint quel que soit
le label pN/pS (dominé par le bruit à faible n) — on utilise snp_sites, pas le ratio.

Sortie : résultats/phase32_phenotype_requalification/{leads.json, priorities.md}. Fusion
`phenotype_lead` dans les fiches hypothétiques concernées + enregistrement phase4. Stdlib.
Run: python analyses/phase32_phenotype_requalification.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "phase32_phenotype_requalification"
GENES = ROOT / "site" / "content" / "genes"

CELLWALL_DRUGS = {"meropenem", "vancomycin", "ethambutol", "rifampicin", "isoniazid", "ethionamide"}
DORMANCY_TF = {"devR", "dosR", "Rv3133c"}
VIRULENCE_TF = {"phoP", "phoR", "Rv0757"}


def constrained(cons):
    if not cons:
        return False
    if cons.get("selection") == "purifying":
        return True
    ss = cons.get("snp_sites")
    return ss is not None and ss <= 10   # quasi-invariant sur 145k = contraint


def build_lead(d):
    mp = d.get("mutant_phenotypes") or {}
    phs = mp.get("phenotypes") or []
    if not phs:
        return None
    cats = {p["category"] for p in phs}
    loc = (d.get("localization") or {}).get("type")
    anchor = (d.get("string") or {}).get("anchor") or {}
    regs = [(r.get("tf_gene") or r.get("tf")) for r in (d.get("regulation") or {}).get("regulated_by", [])]
    operon = [m.get("gene") or m.get("locus") for m in (d.get("genomic_context") or {}).get("operon", [])]
    operon = [x for x in operon if x != d["rv"]]

    # direction fonctionnelle depuis le phénotype
    bits = []
    drugs = sorted({p["condition"].split("under ")[-1] for p in phs if p["category"] == "drug exposure"})
    if drugs:
        cw = [x for x in drugs if x.lower() in CELLWALL_DRUGS]
        bits.append(f"disruption sensitises to {', '.join(drugs)}" +
                    (" (cell-envelope / intrinsic drug tolerance)" if cw else ""))
    if "in vivo" in cats or "macrophage" in cats:
        bits.append("required for fitness in vivo (virulence / persistence factor)")
    stresses = sorted({p["condition"].split("under ")[-1] for p in phs
                       if p["category"] == "stress"})
    if stresses:
        bits.append(f"required under {', '.join(stresses)}")
    if any(p["category"] == "carbon source" for p in phs):
        bits.append("altered fitness on cholesterol (lipid catabolism)")

    # corroboration (indices que c'est un vrai gène fonctionnel)
    corr = []
    if constrained(d.get("conservation")):
        corr.append("conserved / under constraint intra-MTBC")
    if (d.get("essentiality") or {}).get("essential"):
        corr.append("Tn-seq growth-defect/essential")
    if loc in ("SP", "SP+TM"):
        bits.append("predicted secreted (signal peptide)")
    elif loc == "TM":
        bits.append("predicted membrane protein")
    if anchor.get("rv"):
        corr.append(f"STRING-coupled to {anchor.get('gene') or anchor.get('rv')}"
                    + (f" ({anchor.get('product')})" if anchor.get("product") else ""))
    reg_notable = [t for t in regs if t in DORMANCY_TF] and "DosR dormancy regulon" or \
                  ([t for t in regs if t in VIRULENCE_TF] and "PhoP virulence regulon" or "")
    if reg_notable:
        corr.append(f"in the {reg_notable}")
    if len(operon) >= 1:
        named = [x for x in operon if not x.startswith("Rv")]
        if named:
            corr.append(f"co-transcribed with {', '.join(named[:3])}")
    if (d.get("struct_af") or {}).get("hits") or (d.get("mcsa") or {}).get("verdict"):
        corr.append("structural lead available")

    if not bits:
        return None
    lead = "; ".join(bits) + "."

    # score de priorité
    score = 0.0
    invivo_lfcs = [p["log2fc"] for p in phs if p["category"] in ("in vivo", "macrophage")]
    if invivo_lfcs:
        score += min(3.0, abs(min(invivo_lfcs)))          # défaut in vivo le plus fort
    score += min(2, sum(1 for p in phs if p["category"] == "drug exposure"))
    score += min(1, sum(1 for p in phs if p["category"] == "stress"))
    score += min(4, len(corr))                             # corroboration
    return {
        "lead": lead,
        "priority": round(score, 1),
        "categories": sorted(cats),
        "corroboration": corr,
        "n_significant": mp.get("n_significant"),
        "n_in_vivo": mp.get("n_in_vivo"),
    }


def main():
    print("== phase32 : requalification des hypothétiques par phénotype (P5.2b) ==")
    leads = {}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        if "hypothetical" not in (d.get("product_h37rv") or "").lower():
            continue
        lead = build_lead(d)
        if lead:
            leads[d["rv"]] = lead

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "leads.json").write_text(json.dumps(leads, ensure_ascii=False, indent=1))
    print(f"Écrit {OUT/'leads.json'} ({len(leads)} hypothétiques avec lead phénotype)")

    ranked = sorted(leads.items(), key=lambda kv: -kv[1]["priority"])
    md = ["# Requalification des hypothétiques par phénotype — priorisation (P5.2b)", "",
          f"{len(leads)} hypothétiques avec un lead phénotype-driven. Top 30 par priorité.",
          "Lead = HYPOTHÈSE (phénotype + corroboration), à valider.", "",
          "| Rang | Rv | prio | catégories | lead | corroboration |",
          "|---|---|---|---|---|---|"]
    for i, (rv, l) in enumerate(ranked[:30], 1):
        md.append(f"| {i} | {rv} | {l['priority']} | {','.join(l['categories'])} | {l['lead']} | "
                  f"{'; '.join(l['corroboration'])} |")
    (OUT / "priorities.md").write_text("\n".join(md) + "\n")
    print("Top 8 :")
    for rv, l in ranked[:8]:
        print(f"  {rv} (prio {l['priority']}) : {l['lead'][:70]}")

    # fusion dans les fiches
    n = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        if d["rv"] in leads:
            d["phenotype_lead"] = leads[d["rv"]]
            Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
            n += 1
    print(f"Fusionné 'phenotype_lead' dans {n} fiches")


if __name__ == "__main__":
    main()
