#!/usr/bin/env python3
"""phase8_propagate_ta_curation.py

Back-propagate the curated verdicts of the companion study
(TA_persistence: "Toxin-antitoxin systems and persistence in the MTBC,
a structure-guided re-annotation") into the annotation_mtbc *site* gene
fiches (site/content/genes/<rv>.json), the source of truth ingested into
site/data/db.sqlite by `python -m backend.ingest`.

The auto-curation layer (PGAP+Pfam+Foldseek) did not capture the article's
manual findings: the HHpred resolution of two dark ORFs (Rv2664, Rv2308a),
the synteny+AlphaFold3 reconstruction of the two operons, the 16
RefSeq-behind literature identities, and the paralogue-of-fold / non-TA
corrections. This script flips those fiches to the hand-review layer
(auto=false), updating verdict / confidence / function_revised / evidence
and attaching CrossRef-verified references.

Idempotent: re-running overwrites the same fields with the same content.
Run:  python analyses/phase8_propagate_ta_curation.py
Then: cd site && python -m backend.ingest
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

GENES = Path(__file__).resolve().parent.parent / "site" / "content" / "genes"

# ── References (verbatim from TA_persistence/article/references.bib, DOIs CrossRef-verified) ──
COMPANION = {
    "authors": "Guyeux C",
    "year": 2026,
    "title": "Toxin-antitoxin systems and persistence in the Mycobacterium tuberculosis complex: a structure-guided re-annotation",
    "journal": "in preparation",
    "doi": "",
}
REFS = {
    "cai2020": {"authors": "Cai Y, Usher B, Gutierrez C, et al.", "year": 2020, "title": "A nucleotidyltransferase toxin inhibits growth of Mycobacterium tuberculosis through inactivation of tRNA acceptor stems", "journal": "Science Advances", "doi": "10.1126/sciadv.abb6651"},
    "kim2016": {"authors": "Kim Y, Choi E, Hwang J", "year": 2016, "title": "Functional studies of five toxin-antitoxin modules in Mycobacterium tuberculosis H37Rv", "journal": "Frontiers in Microbiology", "doi": "10.3389/fmicb.2016.02071"},
    "akarsu2019": {"authors": "Akarsu H, Bordes P, Mansour M, Bigot D-P, Genevaux P, Falquet L", "year": 2019, "title": "TASmania: a bacterial toxin-antitoxin systems database", "journal": "PLoS Computational Biology", "doi": "10.1371/journal.pcbi.1006946"},
    "kang2023": {"authors": "Kang S-M", "year": 2023, "title": "Mycobacterium tuberculosis Rv0229c shows ribonuclease activity and reveals its corresponding role as toxin VapC51", "journal": "Antibiotics", "doi": "10.3390/antibiotics12050840"},
    "beck2020": {"authors": "Beck IN, Usher B, Hampton HG, Fineran PC, Blower TR", "year": 2020, "title": "Antitoxin autoregulation of Mycobacterium tuberculosis toxin-antitoxin expression through negative cooperativity arising from multiple inverted repeat sequences", "journal": "Biochemical Journal", "doi": "10.1042/BCJ20200368"},
    "yu2020": {"authors": "Yu X, et al.", "year": 2020, "title": "Characterization of a toxin-antitoxin system in Mycobacterium tuberculosis suggests neutralization by phosphorylation as the antitoxicity mechanism", "journal": "Communications Biology", "doi": "10.1038/s42003-020-0941-1"},
    "tandon2019pezat": {"authors": "Tandon H, Sharma A, Sandhya S, Srinivasan N, Singh R", "year": 2019, "title": "Mycobacterium tuberculosis Rv0366c-Rv0367c encodes a non-canonical PezAT-like toxin-antitoxin pair", "journal": "Scientific Reports", "doi": "10.1038/s41598-018-37473-y"},
    "tandon2019ficta": {"authors": "Tandon H, Sharma A, Wadhwa S, Varadarajan R, Singh R, Srinivasan N, Sandhya S", "year": 2019, "title": "Bioinformatic and mutational studies of related toxin-antitoxin pairs in Mycobacterium tuberculosis predict and identify key functional residues", "journal": "Journal of Biological Chemistry", "doi": "10.1074/jbc.RA118.006814"},
    "agarwal2020": {"authors": "Agarwal S, Sharma A, Bouzeyen R, et al.", "year": 2020, "title": "VapBC22 toxin-antitoxin system from Mycobacterium tuberculosis is required for pathogenesis and modulation of host immune response", "journal": "Science Advances", "doi": "10.1126/sciadv.aba6944"},
    "chi2018": {"authors": "Chi X, Chang Y, Li M, et al.", "year": 2018, "title": "Biochemical characterization of mt-PemIK, a novel toxin-antitoxin system in Mycobacterium tuberculosis", "journal": "FEBS Letters", "doi": "10.1002/1873-3468.13280"},
    "sun2015": {"authors": "Sun J, et al.", "year": 2015, "title": "The tuberculosis necrotizing toxin kills macrophages by hydrolyzing NAD", "journal": "Nature Structural and Molecular Biology", "doi": "10.1038/nsmb.3064"},
}


def R(*keys, companion=True):
    refs = [REFS[k] for k in keys]
    if companion:
        refs.append(COMPANION)
    return refs


# ── Tier 1: the two reconstructed operons (HHpred + synteny + AlphaFold3) ──
CUR: dict[str, dict] = {}

CUR["Rv2663"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "BrnT/RelE-superfamily ribonuclease toxin and the toxin half of a HigBA-like module "
        "Rv2663-Rv2664 (toxin-first operon). Left 'hypothetical protein' by RefSeq, it is re-annotated "
        "here on three concordant lines of evidence: an ESMFold model matching a type II TA-complex toxin "
        "fold (Foldseek 7vd7, prob 1.00), MTBC-restricted operon synteny (Rv2664 11 bp downstream, 100% "
        "identity in M. bovis; no detectable M. canettii orthologue, i.e. a young MTBC-specific module), "
        "and an AlphaFold3 heterodimer with the antitoxin Rv2664 predicting a credible interface "
        "(inter-chain ipTM ~0.62). The complex is an in-silico prediction, not a biochemical assay."),
    evidence=[
        "RefSeq/MTBC0 PGAP product: hypothetical protein",
        "Foldseek: 7vd7 type II TA-complex toxin fold (prob 1.00, E=4e-3, TM=0.66)",
        "Synteny (tblastn): 100% identity in M. bovis, partner Rv2664 11 bp downstream; MTBC-restricted (no M. canettii orthologue)",
        "AlphaFold3 Rv2663-Rv2664 heterodimer: inter-chain ipTM 0.61-0.63 (credible interface)",
        "Curated against the companion TA re-annotation (Guyeux 2026)"],
    references=R(companion=True),
)
CUR["Rv2664"] = dict(
    verdict="requalified", confidence="medium",
    function_revised=(
        "HTH/HigA-like antitoxin and the antitoxin half of the HigBA-like module Rv2663-Rv2664. RefSeq "
        "leaves it 'hypothetical protein' and single-structure search (Foldseek/ESMFold) failed to assign "
        "it (best hits non-significant). Profile-profile search (HHpred) confidently matches HigA-family "
        "ribbon-helix-helix/HTH antitoxins (probability ~99%), resolving a fold that structure search "
        "missed. Supported by operon synteny (immediately downstream of the toxin Rv2663, 100% identity in "
        "M. bovis) and by the AlphaFold3 heterodimer with Rv2663 (inter-chain ipTM ~0.62)."),
    evidence=[
        "RefSeq/MTBC0 PGAP product: hypothetical protein",
        "Foldseek/ESMFold: no significant hit (dark to single-structure search)",
        "HHpred: HigA-family HTH antitoxin, probability ~99% (resolves the fold)",
        "Synteny: 11 bp downstream of toxin Rv2663; 100% identity in M. bovis",
        "AlphaFold3 Rv2663-Rv2664: inter-chain ipTM ~0.62",
        "Curated against the companion TA re-annotation (Guyeux 2026)"],
    references=R(companion=True),
)
CUR["Rv2308"] = dict(
    verdict="family_assigned", confidence="medium",
    function_revised=(
        "Putative antitoxin (DUF433 + HTH) and a certain paralogue of the validated antitoxin Rv2018; "
        "proposed antitoxin half of an ancient Rv2308-Rv2308a operon. The operon is conserved beyond the "
        "MTBC to the M. canettii outgroup, with the antitoxin (96% identity) more conserved than its "
        "faster-evolving putative toxin partner Rv2308a (45%), the canonical asymmetry of a genuine TA "
        "pair. The toxin-antitoxin pairing nonetheless remains a hypothesis: AlphaFold3 predicts no "
        "Rv2308-Rv2308a complex (inter-chain ipTM ~0.11)."),
    evidence=[
        "RefSeq/MTBC0 PGAP product: DUF433 domain-containing protein",
        "Paralogue of the validated antitoxin Rv2018 (DUF433 + HTH)",
        "Synteny: conserved to M. canettii (96% id), ancient operon predating the MTBC",
        "AlphaFold3: no Rv2308-Rv2308a interface (ipTM ~0.11), pairing remains hypothetical",
        "Curated against the companion TA re-annotation (Guyeux 2026)"],
    references=R(companion=True),
)
CUR["Rv2308a"] = dict(
    verdict="family_assigned", confidence="low",
    function_revised=(
        "Putative VapC/PIN-like ribonuclease toxin, proposed partner of the antitoxin Rv2308 in an ancient "
        "operon. Dark to Foldseek/ESMFold; HHpred matches ribonuclease/PIN-like toxin profiles "
        "(probability ~99%), i.e. distant family-level homology. Conserved to the M. canettii outgroup "
        "(45% identity, full length, 44 bp downstream of Rv2308) and present at 100% identity in M. bovis "
        "although unannotated there, so not an H37Rv annotation artefact. The pairing is tentative: "
        "AlphaFold3 yields neither a confident monomer fold (pTM ~0.20) nor a Rv2308-Rv2308a complex "
        "(ipTM ~0.11); presented as an unconfirmed hypothesis."),
    evidence=[
        "Not annotated as a gene by RefSeq/PGAP (small ORF downstream of Rv2308)",
        "Foldseek/ESMFold: no significant hit (dark)",
        "HHpred: ribonuclease/PIN-like toxin profile, probability ~99% (distant family homology)",
        "Synteny: 45% id to M. canettii (full length); 100% id in M. bovis though unannotated (not an H37Rv artefact)",
        "AlphaFold3: poor monomer fold (pTM ~0.20), no complex (ipTM ~0.11), pairing unconfirmed",
        "Curated against the companion TA re-annotation (Guyeux 2026)"],
    references=R(companion=True),
)

# ── Tier 2: the 16 RefSeq-behind (literature identity left generic/hypothetical by RefSeq) ──
CUR["Rv0836c"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "MenT2 (MenAT2 system), a DUF1814 nucleotidyltransferase toxin of the MenAT family. RefSeq leaves "
        "this locus 'hypothetical protein'; the toxin identity is established in the literature (Cai 2020). "
        "Paralogue-of-fold caution: the ESMFold model matches a MenAT1 complex (Foldseek 8an5) by paralogy "
        "- the locus is MenT2, not MenT1; the assignment rests on genomic context, not the fold name."),
    evidence=["RefSeq: hypothetical protein", "Literature: MenT2 toxin of the MenAT2 module (Cai 2020)",
              "Paralogue-of-fold caution: Foldseek hit is MenAT1 (8an5); locus is MenT2 by genomic context"],
    references=R("cai2020"),
)
CUR["Rv3749c"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "VapC50, a PIN-domain ribonuclease toxin of a VapBC type II module (Pfam VapC50_C PF26343). RefSeq "
        "leaves it 'hypothetical protein'; toxin activity was characterised by Kim 2016."),
    evidence=["RefSeq: hypothetical protein", "Pfam: VapC50_C (PF26343)", "Literature: VapC50 ribonuclease toxin (Kim 2016)"],
    references=R("kim2016"),
)
CUR["Rv0268c"] = dict(
    verdict="family_assigned", confidence="medium",
    function_revised=(
        "Predicted Phd/YefM-type antitoxin (Pfam PhdYeFM_antitox PF02604). RefSeq leaves it 'hypothetical "
        "protein'; the antitoxin architecture is predicted (TASmania, Akarsu 2019) and awaits experimental "
        "confirmation of its cognate toxin."),
    evidence=["RefSeq: hypothetical protein", "Pfam: PhdYeFM_antitox (PF02604)", "Prediction: Phd/YefM antitoxin (TASmania, Akarsu 2019)"],
    references=R("akarsu2019"),
)
CUR["Rv0634A"] = dict(
    verdict="family_assigned", confidence="medium",
    function_revised=(
        "Predicted VapB-type antitoxin (Pfam VapB_antitoxin PF09957). RefSeq leaves it 'hypothetical "
        "protein'; antitoxin architecture predicted (TASmania, Akarsu 2019)."),
    evidence=["RefSeq: hypothetical protein", "Pfam: VapB_antitoxin (PF09957)", "Prediction: VapB antitoxin (TASmania, Akarsu 2019)"],
    references=R("akarsu2019"),
)
CUR["Rv0229c"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "VapC51, a PIN-domain ribonuclease toxin of a VapBC module (Pfam PIN PF01850). RefSeq annotates "
        "only 'PIN domain nuclease'; ribonuclease activity and the VapC51 identity were shown by Kang 2023."),
    evidence=["RefSeq: PIN domain nuclease", "Pfam: PIN (PF01850)", "Literature: VapC51 ribonuclease toxin (Kang 2023)"],
    references=R("kang2023"),
)
CUR["Rv0837c"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "MenA2, the type IV antitoxin of the MenAT2 module (AbiEi family, Pfam AbiEi_2). RefSeq names only "
        "'type IV TA AbiEi family antitoxin'; the MenA2 identity and autoregulation are established "
        "(Cai 2020; Beck 2020)."),
    evidence=["RefSeq: type IV TA AbiEi family antitoxin", "Pfam: AbiEi_2", "Literature: MenA2 antitoxin (Cai 2020; Beck 2020)"],
    references=R("cai2020", "beck2020"),
)
CUR["Rv1044"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "MenA3 / TakA, a kinase antitoxin of the MenAT3 module that neutralises its cognate toxin by "
        "phosphorylation. RefSeq names only 'AbiEi family antitoxin domain protein'; the kinase mechanism "
        "was shown by Yu 2020."),
    evidence=["RefSeq: AbiEi family antitoxin domain protein", "Literature: MenA3/TakA kinase antitoxin (Yu 2020)"],
    references=R("yu2020"),
)
CUR["Rv1045"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "MenT3 (TglT), a DUF1814 nucleotidyltransferase toxin of the MenAT3 module that blocks translation "
        "by tRNA inactivation. RefSeq names only the generic 'AbiEii/AbiGii toxin family'; the MenT3 toxin "
        "identity is established (Cai 2020)."),
    evidence=["RefSeq: AbiEii/AbiGii toxin family", "Literature: MenT3/TglT tRNA-inactivating toxin (Cai 2020)"],
    references=R("cai2020"),
)
CUR["Rv2826c"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "MenT4, a DUF1814 nucleotidyltransferase toxin of the MenAT4 module. RefSeq names only the generic "
        "'AbiEii/AbiGii toxin family'; the MenT4 identity is established (Cai 2020)."),
    evidence=["RefSeq: AbiEii/AbiGii toxin family", "Literature: MenT4 toxin (Cai 2020)"],
    references=R("cai2020"),
)
CUR["Rv2827c"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "MenA4, the type IV antitoxin of the MenAT4 module (AbiEi family). RefSeq names only 'type IV TA "
        "AbiEi family antitoxin'; the MenA4 identity is established (Beck 2020)."),
    evidence=["RefSeq: type IV TA AbiEi family antitoxin", "Literature: MenA4 antitoxin (Beck 2020)"],
    references=R("beck2020"),
)
CUR["Rv0366c"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "PezT/zeta-family toxin (Pfam Zeta_toxin PF06414) of a non-canonical PezAT-like module with the "
        "antitoxin Rv0367c. RefSeq names only 'zeta toxin family protein'; the PezAT-like pair was "
        "characterised by Tandon 2019."),
    evidence=["RefSeq: zeta toxin family protein", "Pfam: Zeta_toxin (PF06414)", "Literature: PezT-like toxin, pair with Rv0367c (Tandon 2019)"],
    references=R("tandon2019pezat"),
)
CUR["Rv3642c"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "VbhA/FicA-family antitoxin (Pfam VbhA PF18495) of a FicTA module. RefSeq names only 'VbhA family "
        "protein'; the FicTA pairing was characterised by Tandon 2019."),
    evidence=["RefSeq: VbhA family antitoxin", "Pfam: VbhA (PF18495)", "Literature: FicTA antitoxin (Tandon 2019)"],
    references=R("tandon2019ficta"),
)
CUR["Rv2375"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "VapB22, the antitoxin of the VapBC22 type II module (Pfam VapB PF23719). RefSeq names only 'VapB "
        "family antitoxin'; VapBC22 is required for pathogenesis (Agarwal 2020)."),
    evidence=["RefSeq: VapB family antitoxin", "Pfam: VapB (PF23719)", "Literature: VapB22, VapBC22 required for pathogenesis (Agarwal 2020)"],
    references=R("agarwal2020"),
)
CUR["Rv3098A"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "mt-PemK (mt-PemIK / MazF-mt11), an endoribonuclease toxin (Pfam PemK_toxin PF02452). RefSeq names "
        "only 'PemK-like protein'; biochemical characterisation by Chi 2018. This is the genuine mt-PemK "
        "locus, distinct from the orphan PemK-fold gene Rv2405."),
    evidence=["RefSeq: PemK-like protein", "Pfam: PemK_toxin (PF02452)", "Literature: mt-PemK endoribonuclease toxin (Chi 2018)",
              "Distinct from the orphan PemK-fold Rv2405"],
    references=R("chi2018"),
)
CUR["Rv3188"] = dict(
    verdict="family_assigned", confidence="medium",
    function_revised=(
        "Predicted Xre/MbcA/ParS-type antitoxin, proposed antitoxin of a RES-Xre pair Rv3188-Rv3189. "
        "RefSeq names only the toxin-binding-domain protein; the pairing is a structural prediction "
        "awaiting experimental confirmation."),
    evidence=["RefSeq: Xre/MbcA/ParS toxin-binding-domain protein", "Prediction: antitoxin of a RES-Xre pair Rv3188-Rv3189 (structural)"],
    references=R(companion=True),
)
CUR["Rv3902c"] = dict(
    verdict="requalified", confidence="high",
    function_revised=(
        "IFT, the immunity factor for the tuberculosis necrotizing toxin (TNT), the antitoxin/immunity "
        "partner of the CpnT-TNT system (Pfam Imm61 PF15598). The TNT-IFT pair was characterised by "
        "Sun 2015. (RefSeq already names this locus; included for completeness of the CpnT-TNT module.)"),
    evidence=["RefSeq: necrotizing toxin immunity factor IFT", "Pfam: Imm61 (PF15598)", "Literature: immunity factor for TNT, CpnT-TNT system (Sun 2015)"],
    references=R("sun2015"),
)

# ── Tier 3: paralogue-of-fold caution + the 5 candidates evaluated and excluded as non-TA ──
CUR["Rv2405"] = dict(
    verdict="family_assigned", confidence="low",
    function_revised=(
        "Orphan PemK/MazF-fold protein (Pfam PemK_toxin PF02452). Surfaced as a TA candidate by its "
        "endoRNase fold, but it is NOT the characterised mt-PemK toxin (that is Rv3098A) and it has no "
        "adjacent cognate antitoxin. Reported as a paralogue-of-fold caution: a PemK-fold does not by "
        "itself establish an active TA toxin. Pairing and activity unconfirmed."),
    evidence=["Pfam: PemK_toxin (PF02452)", "Paralogue-of-fold caution: not mt-PemK (Rv3098A); no adjacent antitoxin",
              "TA candidate evaluated and left unconfirmed (Guyeux 2026)"],
    references=R(companion=True),
)

# The 5 non-TA: keep the existing (correct, non-TA) function, only flag as reviewed-and-excluded,
# except Rv2049c whose dark fiche is enriched into an explicit non-TA regulator call.
NON_TA_NOTE = "TA candidate evaluated and excluded: this is not a toxin-antitoxin gene (Guyeux 2026 companion)."
NON_TA_KEEP = ["Rv1662", "Rv3105c", "Rv1717", "Rv2619c"]
CUR["Rv2049c"] = dict(
    verdict="family_assigned", confidence="low",
    function_revised=(
        "Putative FrmR/CsoR-family transcriptional regulator (ESMFold model matches a formaldehyde-"
        "responsive regulator, Foldseek 5lcy). Surfaced as a TA candidate but re-annotated here as a "
        "non-TA regulator: no cognate toxin/antitoxin partner and no TA-family domain. Excluded from the "
        "TA set."),
    evidence=["Foldseek: 5lcy formaldehyde-responsive regulator (FrmR-like)", NON_TA_NOTE],
    references=R(companion=True),
)


def apply_fiche(rv: str, patch: dict, keep_function: bool = False) -> str:
    fp = GENES / f"{rv}.json"
    if not fp.exists():
        return f"  MISSING {rv}.json"
    d = json.loads(fp.read_text())
    before = json.dumps(d, sort_keys=True)
    d["auto"] = False
    d["needs_review"] = False
    for k in ("verdict", "confidence", "evidence", "references"):
        if k in patch:
            d[k] = patch[k]
    if not keep_function and "function_revised" in patch:
        d["function_revised"] = patch["function_revised"]
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
    print(f"\nTier 1+2+3 explicit curations ({len(CUR)} fiches):")
    for rv, patch in CUR.items():
        print(apply_fiche(rv, patch))
    print(f"\nNon-TA exclusion note (keep existing function, flag reviewed) ({len(NON_TA_KEEP)} fiches):")
    for rv in NON_TA_KEEP:
        fp = GENES / f"{rv}.json"
        d = json.loads(fp.read_text())
        ev = list(d.get("evidence") or [])
        if NON_TA_NOTE not in ev:
            ev.append(NON_TA_NOTE)
        patch = dict(evidence=ev, verdict=d.get("verdict") or "requalified",
                     confidence=d.get("confidence") or "medium")
        print(apply_fiche(rv, patch, keep_function=True))
    total = len(CUR) + len(NON_TA_KEEP)
    print(f"\nDone. {total} fiches curated (now auto=false). "
          f"Re-ingest with: cd site && python -m backend.ingest")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
