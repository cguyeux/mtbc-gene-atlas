#!/usr/bin/env python3
"""phase9_propagate_dark_enzymes_curation.py

Back-propagate the curated verdicts of the companion study
(dark_enzymes: "Structure-guided functional hypotheses for uncharacterised
enzymes of Mycobacterium tuberculosis") into the annotation_mtbc *site* gene
fiches (site/content/genes/<rv>.json), the source of truth ingested into
site/data/db.sqlite by `python -m backend.ingest`.

Mirror of phase8_propagate_ta_curation.py. The auto-curation layer
(PGAP+Pfam+Foldseek) did not capture the article's manual findings: the 3D
active-site mapping, the HHpred fold confirmations, the AlphaFold3 cross-
validation (incl. the Rv3577 metal co-folding that confirmed the binuclear
site), the fold-paralogue safeguards, and the requalification of the Rv3196
structural false positive. This script flips those six fiches to the
hand-review layer (auto=false), updating verdict / confidence /
function_revised / evidence and attaching CrossRef-verified references.

Idempotent: re-running overwrites the same fields with the same content.
Run:  python analyses/phase9_propagate_dark_enzymes_curation.py
Then: cd site && python -m backend.ingest
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

GENES = Path(__file__).resolve().parent.parent / "site" / "content" / "genes"

# ── References (verbatim from dark_enzymes/article/references.bib, DOIs verified) ──
COMPANION = {
    "authors": "Guyeux C",
    "year": 2026,
    "title": "Structure-guided functional hypotheses for uncharacterised enzymes of Mycobacterium tuberculosis",
    "journal": "in preparation",
    "doi": "10.5281/zenodo.20571950",
}
REFS = {
    "quirk2026": {
        "authors": "Quirk NF, Gregory KN, Morita YS, Tan S",
        "year": 2026,
        "title": "Rv3839-Rv3840 links the endogenous heme biosynthesis pathway with Mycobacterium tuberculosis adaptation to nitric oxide and iron limitation stress",
        "journal": "bioRxiv (preprint)",
        "doi": "10.64898/2026.02.17.706279",
    },
}


def R(*keys, companion=True):
    refs = [REFS[k] for k in keys]
    if companion:
        refs.append(COMPANION)
    return refs


CUR: dict[str, dict] = {}

# ── Rv1118c: permuted NlpC/P60 amidase, competent Cys234-His97 dyad ──
CUR["Rv1118c"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "Circularly permuted NlpC/P60 (YaeF/YiiX-family) cysteine amidase with a structurally competent "
        "Cys234-His97 catalytic dyad (Glu115 a candidate third member). RefSeq leaves this locus "
        "'hypothetical protein'; here it is re-annotated by structure-guided active-site mapping. The "
        "catalytic histidine partner of the nucleophile Cys234 is His97 (Sgamma-Ndelta1 = 3.6 A on the "
        "ESMFold model, 2.9 A on AlphaFold3), invisible to sequence proximity because the fold is "
        "circularly permuted - the three sequence-proximal histidines lie 7-23 A away. HHpred matches the "
        "permuted Peptidase_C92 / YaeF-YiiX family (>=99.7%) and lipid-acting permuted members "
        "(LRAT, H-REV107), and the 86%-hydrophobic pocket points to an N-acyl-amino-acid / lipoprotein "
        "amide substrate rather than a peptidoglycan muropeptide. Distinct from the five canonical "
        "peptidoglycan-hydrolase NlpC/P60 enzymes of H37Rv (Rv0024, RipA, RipB, RipD, Rv2190c). Catalytic "
        "codons essentially invariant across ~250,724 MTBC genomes. A structural prediction, not a "
        "biochemical assay."),
    evidence=[
        "RefSeq: hypothetical protein",
        "Active-site mapping: catalytic dyad Cys234-His97 (Sgamma-Ndelta1 3.6 A ESMFold / 2.9 A AlphaFold3); circularly permuted topology",
        "HHpred: permuted Peptidase_C92 / YaeF-YiiX family >=99.7%; lipid-acting permuted members (LRAT, H-REV107)",
        "Hydrophobic pocket (86%) -> N-acyl-amino-acid / lipoprotein amide substrate, not peptidoglycan",
        "Distinct from the 5 canonical peptidoglycan-hydrolase NlpC/P60 enzymes of H37Rv",
        "Catalytic codons invariant across ~250,724 MTBC genomes (purifying selection)",
        "Curated against the companion dark-enzymes re-annotation (Guyeux 2026)"],
    references=R(companion=True),
)

# ── Rv3577: binuclear MBL-fold metallo-hydrolase, metal site confirmed by AF3 co-folding ──
CUR["Rv3577"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "Binuclear metallo-beta-lactamase (MBL)-fold metallo-hydrolase of the UPF0173/UlaG family "
        "(InterPro IPR050114), substrate unassigned. RefSeq leaves it 'hypothetical protein'. The MBL "
        "HxHxDH motif (His74-His76-Asp78-His79) was confirmed as a genuine two-metal centre by co-folding "
        "the chain with two metal ions on AlphaFold Server (independent Zn and Fe jobs, top model "
        "iPTM 0.98): the ions bind a single bridged binuclear site 3.2-3.4 A apart, partitioned into a "
        "three-histidine metal (His74/His76/His137) and an aspartate-plus-two-histidine metal "
        "(Asp78/His79/His235), with identical geometry for Zn or Fe. HHpred top hits (>=99.8%) are "
        "UPF0173/UlaG metal-dependent hydrolases (UlaG 2wyl; COG2220). The fold-paralogue safeguard "
        "withholds the RNase Z (held by Rv2407) and glyoxalase II (Rv0634c/Rv2581c) labels: only the "
        "fold-level family is claimed. The six metal ligands are effectively invariant across ~250,724 "
        "genomes (most frequent non-synonymous ligand variant 0.0032%). A structural prediction, not a "
        "biochemical assay."),
    evidence=[
        "RefSeq: hypothetical protein; eggNOG COG2220 (Zn-dependent beta-lactamase-fold hydrolases)",
        "MBL motif HxHxDH (His74/His76/Asp78/His79); Foldseek 3bv6 / 2wyl (prob 1.00)",
        "AlphaFold Server metal co-folding (Zn + Fe): bridged binuclear site 3.2-3.4 A; ligands His74/His76/His137 + Asp78/His79/His235; iPTM 0.98",
        "HHpred >=99.8%: UPF0173/UlaG metal-dependent hydrolases (UlaG 2wyl, COG2220)",
        "Fold-paralogue safeguard: RNase Z (Rv2407) and glyoxalase II (Rv0634c/Rv2581c) labels withheld",
        "Six metal-ligand codons effectively invariant across ~250,724 genomes (most frequent NS 0.0032%)",
        "Curated against the companion dark-enzymes re-annotation (Guyeux 2026)"],
    references=R(companion=True),
)

# ── Rv2492: TS-superfamily (hydroxymethyl)transferase, deoxyuridylate-type substrate ──
CUR["Rv2492"] = dict(
    verdict="requalified", confidence="medium",
    function_revised=(
        "Thymidylate-synthase (TS)-superfamily (hydroxymethyl)transferase acting on a deoxyuridylate-type "
        "nucleotide, a distinct non-essential locus separate from the essential thymidylate synthases "
        "ThyA/ThyX. RefSeq leaves it 'hypothetical protein'. The model superposes on MilA (a CMP "
        "hydroxymethyltransferase, PDB 5b6d) with a complete catalytic constellation at the nucleotide "
        "(nucleophile Cys151, phosphate-binding Arg175, base-reading Asp178, ribose-binding His222/Tyr224); "
        "HHpred places it 100% in the thymidylate-synthase / dCMP-hydroxymethylase superfamily. The "
        "specificity residues reassign the substrate away from the template's cytidylate: Asn186 (vs MilA "
        "Asp186) and Ser176 (deoxy-type) point to a (deoxy)uridylate (dUMP-type) substrate with one-carbon "
        "chemistry on the pyrimidine C5. Vertically conserved through the human and animal MTBC and "
        "M. canettii (an earlier anti-phage-defence hypothesis is rejected). A structural prediction; "
        "substrate narrowed, not fixed."),
    evidence=[
        "RefSeq: hypothetical protein",
        "Superposition on MilA+CMP (5b6d); HHpred 100% thymidylate-synthase / dCMP-hydroxymethylase superfamily",
        "Competent catalytic core: Cys151, Arg175, Asp178, His222, Tyr224 at H-bonding distance of the nucleotide",
        "Specificity residues Asn186 + Ser176 -> (deoxy)uridylate (dUMP-type), not the template CMP",
        "Distinct, non-essential locus (25% id); not redundant with ThyA/ThyX",
        "Curated against the companion dark-enzymes re-annotation (Guyeux 2026)"],
    references=R(companion=True),
)

# ── Rv2558: cofactor-independent ABM monooxygenase (weaker candidate) ──
CUR["Rv2558"] = dict(
    verdict="family_assigned", confidence="low",
    function_revised=(
        "Cofactor-independent ABM (antibiotic-biosynthesis-monooxygenase) of the quinone-monooxygenase "
        "family (ActVA-Orf6/SnoaB/YgiN), substrate unknown (weaker candidate). RefSeq leaves it "
        "'hypothetical protein'. Tandem-ABM protein (paralogue Rv2557, 69% identity); only the C-terminal "
        "lobe forms a competent pocket, centred on His124 and lined by Phe168/Tyr154 - the signature of a "
        "cofactorless ABM site where an aromatic substrate replaces flavin. A heme-oxygenase assignment is "
        "excluded structurally (geometry not heme-compatible; IsdG/MhuD catalytic determinants absent) and "
        "by the fold-paralogue safeguard (heme degradation in M. tuberculosis is already carried by "
        "MhuD/Rv3592). HHpred top hits (>=99.9%) are cofactorless ActVA-Orf6/YgiN-family monooxygenases; "
        "the IsdG/MhuD heme-oxygenase subfamily is absent. The C-terminal site residues are near-invariant "
        "across ~250,724 genomes. Fold and sub-type are well supported; substrate and activity remain "
        "undemonstrated."),
    evidence=[
        "RefSeq: hypothetical protein; tandem ABM, paralogue Rv2557 (69% id)",
        "Competent C-terminal pocket: His124 with Phe168/Tyr154 (cofactorless ABM signature)",
        "HHpred >=99.9%: ActVA-Orf6/YgiN cofactorless monooxygenases; IsdG/MhuD heme-oxygenase subfamily absent",
        "Heme-oxygenase excluded (structure + fold-paralogue safeguard: MhuD/Rv3592 holds that function)",
        "C-terminal site residues near-invariant across ~250,724 genomes; substrate undemonstrated",
        "Curated against the companion dark-enzymes re-annotation (Guyeux 2026)"],
    references=R(companion=True),
)

# ── Rv3196: structural false positive for YcaO, re-routed to ThiF/MoeB/E1 adenylyltransferase ──
CUR["Rv3196"] = dict(
    verdict="family_assigned", confidence="low",
    function_revised=(
        "Probable ThiF/MoeB/E1-family adenylyltransferase of undetermined substrate; NOT a YcaO "
        "heterocyclase. RefSeq leaves it 'hypothetical protein'. A strong Foldseek hit to the cyanobactin "
        "heterocyclase TruD (4bs9, query-normalised TM-score 0.72) is a normalised-TM-score artefact: the "
        "alignment covers only the N-terminal scaffold (TruD residues 3-321), stops before the YcaO "
        "catalytic domain, and the target-normalised TM-score is 0.33 (~38% coverage); none of the YcaO "
        "catalytic residues lies in the aligned region. HHpred re-routes it (>=99.9%) to standalone "
        "adenylyltransferases (MccB, PaaA) and the ThiF/MoeB/E1-like adenylation superfamily, within which "
        "the CxxC pair (C208/C211) matches the MoeB-type zinc motif. A RiPP role is independently "
        "improbable (no precursor-peptide cassette, maturation protease, or TfuA in H37Rv). A worked "
        "example of the normalised-TM-score trap."),
    evidence=[
        "RefSeq: hypothetical protein",
        "Foldseek hit to TruD (4bs9) is a normalised-TM-score artefact: qtmscore 0.72 vs ttmscore 0.33 (~38% coverage); YcaO catalytic domain NOT in the aligned region",
        "HHpred re-routes >=99.9% to ThiF/MoeB/E1 adenylyltransferase superfamily (MccB, PaaA); CxxC (C208/C211) = MoeB-type Zn motif",
        "RiPP / heterocyclase role excluded (no precursor cassette / protease / TfuA / second YcaO in H37Rv)",
        "Scaffold-of-fold caution; conserved in >99.98% of ~250,724 genomes",
        "Curated against the companion dark-enzymes re-annotation (Guyeux 2026)"],
    references=R(companion=True),
)

# ── Rv3839: heme-homeostasis regulator, independently characterised (Quirk 2026) ──
CUR["Rv3839"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "Heme-homeostasis regulator: the principal driver of the Rv3839-Rv3840 operon that down-regulates "
        "endogenous heme biosynthesis under nitric-oxide and iron-limitation stress (Quirk 2026). It "
        "carries a DUF2470/DRI (Domain Related to Iron) regulatory module - the GluTR-binding-protein "
        "domain, also found in HugZ - and lies immediately upstream of bacterioferritin bfrB (Rv3841), in "
        "an operon with the transcriptional regulator Rv3840. The fold-paralogue safeguard is passed: it "
        "is not a copy of the genuine glutamyl-tRNA reductase HemA/Rv0509 but shares only the regulatory "
        "DUF2470 domain. RefSeq left this locus uncharacterised; the structure-guided hypothesis "
        "(heme/iron homeostasis) was confirmed experimentally by Quirk 2026 while the companion study was "
        "in progress - a RefSeq-behind case and an independent validation of the approach."),
    evidence=[
        "RefSeq: hypothetical protein (DUF2470); operon Rv3839-Rv3840 immediately upstream of bfrB (Rv3841)",
        "DUF2470/DRI + GluTR-binding-protein regulatory module (also in HugZ)",
        "Fold-paralogue safeguard: not the genuine glutamyl-tRNA reductase HemA/Rv0509",
        "Literature: Rv3839-Rv3840 down-regulates heme biosynthesis under NO / iron-limitation stress, Rv3839 principal driver (Quirk 2026)",
        "RefSeq-behind case; independent validation of the structure-guided hypothesis",
        "Curated against the companion dark-enzymes re-annotation (Guyeux 2026)"],
    references=R("quirk2026", companion=True),
)


def apply_fiche(rv: str, patch: dict) -> str:
    fp = GENES / f"{rv}.json"
    if not fp.exists():
        return f"  MISSING {rv}.json"
    d = json.loads(fp.read_text())
    before = json.dumps(d, sort_keys=True)
    d["auto"] = False
    d["needs_review"] = False
    for k in ("verdict", "confidence", "function_revised", "evidence", "references"):
        if k in patch:
            d[k] = patch[k]
    after = json.dumps(d, sort_keys=True)
    if before == after:
        return f"  = {rv} (unchanged)"
    fp.write_text(json.dumps(d, indent=2, ensure_ascii=False) + "\n")
    return f"  + {rv} -> {d['verdict']}/{d['confidence']}"


def main() -> int:
    if not GENES.is_dir():
        print(f"ERROR: gene fiches dir not found: {GENES}", file=sys.stderr)
        return 1
    print(f"Gene fiches: {GENES}")
    print(f"\nDark-enzymes curations ({len(CUR)} fiches):")
    for rv, patch in CUR.items():
        print(apply_fiche(rv, patch))
    print(f"\nDone. {len(CUR)} fiches curated (now auto=false). "
          f"Re-ingest with: cd site && python -m backend.ingest")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
