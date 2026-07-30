#!/usr/bin/env python3
"""dark_dossier.py <Rv> [Rv2 ...] -- P8 étape 1 : dump lisible de TOUTES les couches d'un gène.

READ-ONLY. Agrège, pour un locus, l'évidence pertinente à une décision de curation intégrative :
identité, « vrai gène ? » (conservation via snp_sites, essentialité, MS), vulnérabilité, localisation,
contexte d'opéron (voisins NOMMÉS + fonction), régulon, iModulons (co-expression conditionnelle),
réseau STRING (partenaires nommés + CANAUX pour distinguer experimental/database de neighborhood),
phénotypes de mutants, PTM, structure (PDB/AF/Foldseek/pLDDT/pfam_tentative), HHpred, legacy Mycobrowser,
leads déjà posés. C'est l'ENTRÉE de la synthèse intégrative — il n'écrit rien.
Run: python analyses/dark_dossier.py Rv2680
"""
from __future__ import annotations
import json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"


def g(rv):
    p = GENES / f"{rv}.json"
    return json.load(open(p)) if p.exists() else None


def neigh_fn(rv):
    d = g(rv) or {}
    nm = d.get("gene_name") or (d.get("eggnog") or {}).get("preferred_name") or ""
    desc = (d.get("eggnog") or {}).get("description") or (d.get("uniprot") or {}).get("protein_name") or ""
    fr = (d.get("function_revised") or "")
    vd = d.get("verdict")
    detail = nm or (desc if "hypothetical" not in desc.lower() else "") or (fr[:50] if not fr.startswith("Conserved hypothetical") else "hypothetical")
    return f"{rv}({nm or '—'}, {vd}): {detail[:60]}"


def dossier(rv):
    d = g(rv)
    if not d:
        return f"[{rv}] fiche introuvable\n"
    out = [f"================ {rv} ================"]
    seq = d.get("protein_mtbc0") or ""
    out.append(f"len={len(seq)}aa | verdict={d.get('verdict')}/{d.get('confidence')} | product={d.get('product_h37rv')}")
    out.append(f"function_revised: {(d.get('function_revised') or '')[:200]}")

    # --- vrai gène ? ---
    c = d.get("conservation") or {}
    out.append(f"\n[VRAI GÈNE ?] conservation: selection={c.get('selection')} pN/pS={c.get('pN_pS')} "
               f"snp_sites={c.get('snp_sites')} syn/mis/non/fs={c.get('syn')}/{c.get('missense')}/{c.get('nonsense')}/{c.get('frameshift')} "
               f"pseudo={c.get('pseudogene_flag') or c.get('pseudogene_candidate')} (n={c.get('n_strains')})")
    ess = d.get("essentiality") or {}
    out.append(f"  essentiality: {ess.get('dejesus2017')} essential={ess.get('essential')} | "
               f"proteomics MS: detected={(d.get('proteomics') or {}).get('detected')} n_datasets={(d.get('proteomics') or {}).get('n_datasets')}")
    v = d.get("vulnerability") or {}
    out.append(f"  CRISPRi vulnerability VI={v.get('vi')} (95%CI {v.get('ci_low')}..{v.get('ci_high')})  [plus négatif=plus vulnérable]")

    # --- localisation / structure ---
    loc = d.get("localization") or {}
    out.append(f"\n[LOCALISATION] {loc.get('type')} TM={loc.get('tm_helices')} SP={loc.get('signal_peptide')} lipo={loc.get('lipoprotein')} ({loc.get('prediction')})")
    pl = d.get("plddt") or {}; plaf = d.get("plddt_af") or {}
    pdb = d.get("pdb") or {}; sa = d.get("struct_af") or {}; pf = d.get("pfam_tentative")
    out.append(f"[STRUCTURE] pLDDT(ESM)={pl.get('mean_plddt')} pLDDT(AF)={plaf.get('mean_plddt')} "
               f"PDB={pdb.get('n_structures') if isinstance(pdb,dict) else None} "
               f"Foldseek_AFDB_hits={len(sa.get('hits',[])) if isinstance(sa,dict) else 0} pfam_tentative={bool(pf)}")
    if isinstance(sa, dict) and sa.get("hits"):
        for h in sa["hits"][:2]:
            out.append(f"    AFDB: {h.get('description','')[:55]} tm={h.get('tmscore')} sig={h.get('significant')}")
    mcsa = d.get("mcsa")
    if mcsa:
        out.append(f"[M-CSA] {json.dumps(mcsa, ensure_ascii=False)[:120]}")
    hh = d.get("hhpred")
    if hh:
        out.append(f"[HHpred] no_hit={hh.get('no_confident_hit')} note={(hh.get('note') or '')[:120]} top={hh.get('top_hits')}")

    # --- contexte : opéron ---
    gc = d.get("genomic_context") or {}
    op = gc.get("operon") or []
    out.append(f"\n[OPÉRON] size={gc.get('operon_size')} strand={gc.get('strand')}")
    for m in op:
        loci = m.get("locus")
        tag = " <-- CE GÈNE" if loci == rv else ""
        out.append(f"    {neigh_fn(loci) if loci != rv else loci}{tag}")

    # --- régulon / iModulon ---
    reg = d.get("regulation") or {}
    if reg.get("regulated_by"):
        out.append(f"\n[RÉGULON] régulé par: {', '.join((r.get('tf_gene') or r.get('tf')) + '(' + r.get('effect','?') + ')' for r in reg['regulated_by'])}")
    ex = d.get("expression") or {}
    if ex.get("imodulons"):
        out.append(f"[iMODULON co-expression] " + "; ".join(
            f"{m['name']}" + (f"/reg:{m['regulator']}" if m.get('regulator') else "") for m in ex["imodulons"]))

    # --- STRING (avec canaux) ---
    s = d.get("string") or {}
    if s.get("partners"):
        out.append(f"\n[STRING] {s.get('n_partners')} partenaires ; ancre contexte={(s.get('anchor') or {}).get('gene')}")
        named = [p for p in s["partners"] if p.get("gene") and not p.get("hypothetical")]
        for p in sorted(named, key=lambda p: -(p.get("combined_no_tm") or 0))[:5]:
            ch = p.get("channels") or {}
            chans = " ".join(f"{k}:{val}" for k, val in ch.items() if val and val >= 200)
            out.append(f"    {p['gene']}({p['rv']}) no_tm={p.get('combined_no_tm')} [{chans}] :: {(p.get('product') or '')[:40]}")

    # --- phénotype / PTM / mycobrowser / leads ---
    mp = d.get("mutant_phenotypes") or {}
    if mp.get("n_significant"):
        phs = mp.get("phenotypes") or []
        out.append(f"\n[PHÉNOTYPE] {mp['n_significant']} significatifs, {mp.get('n_in_vivo')} in vivo ; ex: " +
                   "; ".join(f"{p.get('condition','')[:40]}(lfc {p.get('log2fc')})" for p in phs[:3]))
    ptm = d.get("ptm") or {}
    if ptm.get("n_sites"):
        out.append(f"[PTM] {ptm['n_sites']} sites ({ptm.get('n_phospho')} phospho)")
    for lead in ("phenotype_lead", "context_lead", "struct_cluster"):
        if d.get(lead):
            out.append(f"[LEAD {lead}] {json.dumps(d[lead], ensure_ascii=False)[:150]}")
    myc = d.get("mycobrowser") or {}
    if myc:
        out.append(f"[MYCOBROWSER legacy] {json.dumps(myc, ensure_ascii=False)[:150]}")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: dark_dossier.py Rv0001 [Rv0002 ...]")
    for rv in sys.argv[1:]:
        print(dossier(rv))
