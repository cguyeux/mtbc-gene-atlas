#!/usr/bin/env python3
"""Changelog par gène : ce que l'atlas a APPORTÉ vs Mycobrowser (et via quelle évidence).

Une ligne par gène (3906). Colonnes : l'annotation de départ (H37Rv + Mycobrowser legacy),
le résultat atlas (verdict + fonction révisée + EC), et surtout le DELTA — ce qui est nouveau :
catégorie de nouveauté, « ahead of Mycobrowser », statut EC, et l'ÉVIDENCE fonctionnelle qui a
permis l'annotation (quelles couches). But : qu'un lecteur voie d'un coup d'œil ce que la
production de l'atlas a changé/ajouté par rapport à la base de référence.

Lecture SEULE. Sortie : résultats/phase65_changelog/atlas_changelog.csv (+ résumé à l'écran).
"""
import csv
import glob
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "résultats" / "phase65_changelog"

WS = re.compile(r"\s+")


def clean(s) -> str:
    return WS.sub(" ", str(s or "")).strip()


def functional_basis(d: dict) -> list[str]:
    """Couches apportant un HANDLE fonctionnel (ce que l'atlas met en face du gène)."""
    b = []
    egg = d.get("eggnog") or {}
    if egg.get("og") or egg.get("ec") or egg.get("kegg_ko") or egg.get("description"):
        b.append("eggNOG")
    up = d.get("uniprot") or {}
    if up.get("reviewed"):
        b.append("SwissProt")
    elif up.get("function"):
        b.append("UniProt")
    if d.get("domains"):
        b.append("Pfam")
    hh = d.get("hhpred") or {}
    if hh and not hh.get("no_confident_hit"):
        b.append("HHpred")
    if d.get("struct_af") or d.get("struct_hits"):
        b.append("Foldseek")
    if d.get("pdb"):
        b.append("PDB")
    mc = d.get("mcsa") or {}
    if mc and (mc.get("active_site") or mc.get("has_active_site")):
        b.append("M-CSA")
    st = d.get("string") or {}
    if st.get("n_partners"):
        b.append("STRING")
    if d.get("context_lead"):
        b.append("operon-context")
    if d.get("phenotype_lead"):
        b.append("phenotype")
    if d.get("integrative_lead"):
        b.append("integrative")
    if d.get("struct_cluster"):
        b.append("struct-cluster")
    return b


def support_evidence(d: dict) -> list[str]:
    """Signaux « gène réel / important » (pas la fonction, mais la nouveauté atlas)."""
    s = []
    ess = d.get("essentiality") or {}
    if ess.get("essential"):
        s.append("essential")
    vul = d.get("vulnerability") or {}
    if isinstance(vul.get("vi"), (int, float)) and vul["vi"] <= -6:
        s.append("vulnerable")
    con = d.get("conservation") or {}
    if "purif" in (con.get("selection") or ""):
        s.append("purifying-selection")
    if con.get("pseudogene_flag"):
        s.append("pseudogene-candidate")
    pr = d.get("proteomics") or {}
    if pr.get("detected") or pr.get("n_datasets"):
        s.append("MS-detected")
    return s


def annotation_route(d: dict) -> str:
    """Par quel CANAL l'atlas est arrivé à l'annotation (littérature vs homologie/structure)."""
    if d.get("verdict") == "dark":
        return "unknown"
    up = d.get("uniprot") or {}
    if up.get("reviewed"):
        return "literature (SwissProt-curated)"
    if d.get("curation_note") or d.get("integrative_lead") or d.get("phenotype_lead"):
        return "curation/synthesis"
    hh = d.get("hhpred") or {}
    if hh and not hh.get("no_confident_hit"):
        return "remote-homology (HHpred)"
    egg = d.get("eggnog") or {}
    if egg.get("og") or egg.get("kegg_ko") or egg.get("ec"):
        return "orthology (eggNOG)"
    if d.get("domains"):
        return "domain (Pfam)"
    if d.get("pdb") or d.get("struct_af") or d.get("struct_hits"):
        return "structure (PDB/Foldseek)"
    if d.get("context_lead"):
        return "genomic-context"
    return "other"


def novelty(d: dict, myc: dict, was_unknown: bool) -> str:
    """Le titre du changement vs Mycobrowser (une catégorie)."""
    verdict = d.get("verdict")
    ec_status = myc.get("ec_status")
    if verdict == "dark":
        return "still_unknown"
    if was_unknown:
        return "new_function" if verdict == "requalified" else "new_family_or_domain"
    if ec_status == "reclassified":
        return "ec_reclassified"
    if ec_status == "differ":
        return "ec_disagreement"
    if myc.get("atlas_ahead"):
        return "handle_added"
    return "confirmed"


def atlas_ec(d: dict, myc: dict) -> str:
    egg = d.get("eggnog") or {}
    up = d.get("uniprot") or {}
    for src in (myc.get("atlas_ec"), egg.get("ec"), up.get("ec")):
        if src:
            return ";".join(src)
    return ""


def main() -> None:
    rows = []
    nov = Counter()
    for f in sorted(glob.glob(str(GENES / "*.json"))):
        d = json.loads(Path(f).read_text())
        myc = d.get("mycobrowser") or {}
        was_unknown = bool(myc.get("conserved_hypothetical"))
        n = novelty(d, myc, was_unknown)
        nov[n] += 1
        rows.append({
            "rv": d.get("rv"),
            "gene_name": d.get("gene") or "",
            "mtbc0": d.get("mtbc0") or "",
            "h37rv_product": clean(d.get("product_h37rv")),
            "mycobrowser_product": clean(myc.get("legacy_product")),
            "mycobrowser_function": clean(myc.get("legacy_function")),
            "mycobrowser_ec": ";".join(myc.get("legacy_ec") or []),
            "was_unknown_in_mycobrowser": int(was_unknown),
            "atlas_verdict": d.get("verdict") or "",
            "atlas_confidence": d.get("confidence") or "",
            "atlas_function": clean(d.get("function_revised")),
            "atlas_ec": atlas_ec(d, myc),
            "novelty": n,
            "ahead_of_mycobrowser": int(bool(myc.get("atlas_ahead"))),
            "ahead_via": ";".join(myc.get("ahead_via") or []),
            "ec_vs_mycobrowser": myc.get("ec_status") or "",
            "annotation_route": annotation_route(d),
            "swissprot_curated": int(bool((d.get("uniprot") or {}).get("reviewed"))),
            "functional_basis": ";".join(functional_basis(d)),
            "support_evidence": ";".join(support_evidence(d)),
            "n_references": len(d.get("references") or []),
        })

    OUT.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0].keys())
    # Double sortie : résultats/ (archive) ET site/content/ (servi par le site, baked dans l'image).
    served = ROOT / "site" / "content" / "atlas_changelog.csv"
    for csv_path in (OUT / "atlas_changelog.csv", served):
        with csv_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)

    print(f"{len(rows)} gènes -> {OUT/'atlas_changelog.csv'}  +  {served}\n")
    print("Répartition 'annotation_route' (canal d'annotation) :")
    for k, v in Counter(r["annotation_route"] for r in rows).most_common():
        print(f"  {k:<32} {v}")
    print()
    print("Répartition 'novelty' (ce qui a changé vs Mycobrowser) :")
    for k, v in nov.most_common():
        print(f"  {k:<24} {v}")
    ahead = sum(r["ahead_of_mycobrowser"] for r in rows)
    ecrec = sum(1 for r in rows if r["ec_vs_mycobrowser"] == "reclassified")
    unknown_resolved = sum(1 for r in rows if r["was_unknown_in_mycobrowser"] and r["atlas_verdict"] != "dark")
    print(f"\n  ahead_of_mycobrowser        {ahead}")
    print(f"  EC reclassifiés             {ecrec}")
    print(f"  hypothétiques Mycobrowser requalifiés (verdict != dark) : {unknown_resolved} / "
          f"{sum(r['was_unknown_in_mycobrowser'] for r in rows)}")


if __name__ == "__main__":
    main()
