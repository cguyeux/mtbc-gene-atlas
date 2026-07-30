#!/usr/bin/env python3
"""phase34_hhpred_curation.py -- requalification des dark par HHpred profil-profil (P7.1).

Applique, gène par gène, l'interprétation des résultats HHpred (path web, exécuté par CG car
l'API est cassée en headless) sur le résidu dark. Discipline (KB) : Probability>95 = confiant ;
nommer le REPLI/domaine, PAS la sous-famille si le hit top n'est pas un membre canonique ; la
convergence de plusieurs hits vers un même repli compte ; hits à E-value élevée = bruit.
Hand-curation traçable (`auto:false` + réf méthode HHpred + accessions Pfam/PDB).

Couche extensible : on ajoute une entrée à HITS à chaque gène rapporté par CG.
Sortie : résultats/phase34_hhpred_curation/hhpred_curation.json + fusion dans les fiches +
enregistrement phase4 (via la même voie que phase33, curation littérature qui gagne).
Run: python analyses/phase34_hhpred_curation.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "phase34_hhpred_curation"
GENES = ROOT / "site" / "content" / "genes"

HHPRED_REF = {"authors": "Zimmermann L, Stephens A, Nam SZ, et al.", "year": 2018,
              "title": "A Completely Reimplemented MPI Bioinformatics Toolkit with a New HHpred Server at its Core",
              "journal": "Journal of Molecular Biology", "doi": "10.1016/j.jmb.2017.12.007"}

# interprétation experte des résultats HHpred rapportés par CG.
# Entrée AVEC verdict/function_revised = requalification ; entrée "no_confident_hit" = HHpred
# tenté mais sans hit confiant (reste dark, on enregistre l'essai pour ne pas le refaire).
HITS = {
    "Rv0610c": {
        "verdict": "family_assigned", "confidence": "medium",
        "function_revised": "N-terminal region matches the SseB N-terminal domain family (Pfam PF07179, HHpred probability ~96%, E=0.09), with converging hits to DUF3234 (PF11572) and the structural-genomics protein TTHA0547 (PDB 2Z0R). A domain/fold is assigned by profile-profile search (HHpred), but both families are poorly characterised, so the precise function in M. tuberculosis remains undefined. (The lower VbhA/Fic-antitoxin hits are high-E-value noise.)",
        "top_hits": ["PF07179 SseB N-term (96%, E=0.09)", "PF11572 DUF3234 (96%, E=0.15)", "2Z0R TTHA0547 unknown fn (96%, E=0.11)"],
    },
    "Rv0698": {
        "no_confident_hit": True,
        "top_hits": ["4NSX U3 snoRNA-assoc (44%, E=30)", "PF29823 Rv2386A Mtb-hypothetical (35%, E=65)", "1FI8 ecotin (20%, E=300)"],
        "note": "HHpred run gave no confident hit (best 44% probability, E=30) — remains a genuine unknown.",
    },
    "Rv3178a": {
        "verdict": "requalified", "confidence": "medium",
        "function_revised": "Minimal nucleotidyltransferase of the polymerase-beta-like superfamily (HHpred COG1669 95.7%, E=0.16, with strongly converging hits to putative/minimal nucleotidyltransferases and to MntA-type nucleotidyltransferase antitoxins). Confident fold assignment. Genomic context does NOT support a toxin-antitoxin antitoxin role (Rv3178a is a singleton flanked by a nitroreductase and an AAA ATPase at 239/287 bp, no adjacent toxin), so most likely a standalone nucleotidyltransferase; the precise substrate/acceptor remains undefined.",
        "top_hits": ["COG1669 predicted nucleotidyltransferase (95.7%, E=0.16)", "SCOP/1WOT minimal nucleotidyltransferase (88-92%)", "6M6U MntA antitoxin NTase fold (86%, E=15)"],
    },
    "Rv2809": {
        "verdict": "family_assigned", "confidence": "low",
        "function_revised": "Tentatively resembles a peptide-chain-release-factor fold: several converging HHpred hits (~56-62%) to the Actinobacteria/Chloroflexi VLRF1 release factor (Pfam PF18859) and to the eRF1 middle domain / eRF1 domain 2. The top hit is an Actinobacteria-specific family, taxonomically apt for M. tuberculosis. LOW confidence: all probabilities are below the confident threshold and E-values are high, so this is a lead (possible translation-termination-related fold) to validate, not an assignment.",
        "top_hits": ["PF18859 acVLRF1 Actinobacteria release factor (62%, E=67)", "eRF1 middle domain (57-59%)", "PF03464 eRF1 domain 2 (39%)"],
    },
    "Rv2644c": {
        "no_confident_hit": True,
        "top_hits": ["3RLG sphingomyelinase-D / PLC-like TIM-barrel (66%, E=37)", "scattered DUF/DHFR/NAC hits, no convergence"],
        "note": "No confident hit (best 66%, E=37, no convergence). Remains dark.",
    },
    "Rv2808": {
        "verdict": "requalified", "confidence": "medium",
        "function_revised": "Belongs to the SdhE/Sdh5/YgfY family (Pfam DUF5669 89.6% full-length, plus strongly converging YgfY-like / SdhE / Sdh5 hits at 83-88%): the FAD-assembly / flavinylation-factor fold of succinate dehydrogenase. The same fold is also used by CptB-type toxin-antitoxin antitoxins. Genomic context is ambiguous (singleton, not adjacent to sdh genes nor to a toxin; flanked by an integrase-domain hypothetical Rv2807 and Rv2809), so the family/fold is confidently assigned but the specific role (SDH flavinylation vs CptB-like antitoxin) is undetermined.",
        "top_hits": ["PF18918 DUF5669 (89.6%, E=6.7, full-length)", "SdhE FAD assembly factor / YgfY fold (86-88%)", "PF03937 Sdh5 flavinator (83%)"],
    },
    "Rv2662": {
        "no_confident_hit": True,
        "top_hits": ["9LVU de novo synthetic (67%, E=38)", "loose T6SS TssA/ImpA/VasJ series (42-59%, E=50-86)", "no convergence"],
        "note": "No confident hit (best 67%, E=38). Note: the cognate antitoxin of the RelS-RelI toxin Rv2663 is Rv2664 (HTH antitoxin), not Rv2662, so the absence of a TA fold here is consistent. Remains dark.",
    },
    "Rv0997a": {
        "verdict": "requalified", "confidence": "high",
        "function_revised": "N-terminal (KOW-like / SH3-barrel) domain of a translation elongation factor P (EF-P) / eIF-5A family protein (Pfam EFP_N, HHpred 96.5%, E=0.15, with a large confident cascade of EF-P / eIF-5A hits at 90-96%). Confident fold assignment: an EF-P-type ribosome-associated translation factor domain (EF-P assists translation of polyproline motifs). Likely a standalone EF-P N-terminal-domain protein / paralogue.",
        "top_hits": ["EF-P N-terminal domain SH3-barrel (96.5%, E=0.15)", "COG0231 Efp EF-P/eIF-5A (94%)", "PF08207 EFP_N KOW-like (93%)"],
    },
    "Rv2706c": {
        "no_confident_hit": True,
        "top_hits": ["GatA amidotransferase / amidase-signature fold, many converging hits but all 28-42% (E=66-150)", "1M22 peptide amidase (42%, E=66)"],
        "note": "Converging but weak lean toward the amidase-signature (GatA / peptide-amidase, SCOP c.117.1.1) fold; all probabilities <42% (E>=53) — below the confident threshold, not assertable (the amidase fold is large and promiscuous). Remains dark.",
    },
    "Rv3190A": {
        "verdict": "family_assigned", "confidence": "low",
        "function_revised": "Tentatively a small DNA-binding module of the ribbon-helix-helix / HTH type: several converging HHpred hits (~72-75%) to ParB C-terminal, the CcdA post-segregation antitoxin, a plasmid ribbon-helix-helix DNA-binding protein and ECF sigma factors. LOW confidence (no hit >95%; the top LOV-domain hit at 91% is a short 18-column match, likely spurious). Flanked by an IS30-family transposase (Rv3191c), so possibly an IS-associated DNA-binding peptide; NOT asserted as a toxin-antitoxin antitoxin (no adjacent toxin). A lead to validate.",
        "top_hits": ["ParB C-terminal / CcdA antitoxin / plasmid RHH DNA-binding (72-75%, E=20-39)", "ECF sigma factor HTH (60%)", "(2YON LOV 91% but 18-col short, discounted)"],
    },
    "Rv3202a": {
        "verdict": "family_assigned", "confidence": "low",
        "function_revised": "Belongs to a conserved Mycobacterium-specific protein family (Pfam PF29473, the 'Rv3202A family', HHpred 99.6%). This is a self-named family (built from Rv3202A and its homologues) with NO known function, so the assignment confirms it is a genuine conserved family (not a spurious ORF) but the molecular function remains unknown (DUF-like).",
        "top_hits": ["PF29473 Rv3202A family, self-named, unknown function (99.6%)"],
    },
    "Rv2386a": {
        "verdict": "family_assigned", "confidence": "medium",
        "function_revised": "Sm-like fold (SCOP b.38.1.6) small protein of the DUF903 / YgdI-YgdR family (a large converging cascade of HHpred hits at 90-94%; plus the self-named Mtb family PF29823 at 99.95%). Confident fold and family assignment. The DUF903 family is generically annotated as putative lipoproteins, but our own DeepTMHMM prediction sees Rv2386a as globular (no signal peptide / lipobox), so the lipoprotein label is NOT supported here. The molecular function of this family is not established.",
        "top_hits": ["PF29823 Rv2386A self-family (99.95%)", "b.38.1.6 Sm-like fold / DUF903 YgdI-YgdR (90-94%, E~1-4, 13+ converging hits)"],
    },
    "Rv2307A": {
        "no_confident_hit": True,
        "top_hits": ["eukaryotic homeodomain Ultrabithorax (47%, E=35)", "scattered low hits"],
        "note": "No confident hit (best 47%, E=35, a eukaryotic developmental homeodomain — implausible). Remains dark.",
    },
    "Rv0115a": {
        "no_confident_hit": True,
        "top_hits": ["RIF1 telomere regulator eukaryotic (71%, E=20)", "scattered low hits, no convergence"],
        "note": "No confident hit (best 71%, E=20, a eukaryotic telomere protein — implausible; no convergence). Remains dark.",
    },
    "Rv2256a": {
        "no_confident_hit": True,
        "top_hits": ["all hits <29% probability, E>=130 (YuaX, DUF5805, GntR)"],
        "note": "No confident hit (best 29%, E=190) — pure noise. Remains dark.",
    },
    "Rv0699": {
        "no_confident_hit": True,
        "top_hits": ["PF08557 sphingolipid desaturase (51%, E=42)", "scattered low hits, no convergence"],
        "note": "No confident hit (best 51%, E=42, isolated, no convergence). Remains dark.",
    },
    "Rv0397A": {
        "no_confident_hit": True,
        "top_hits": ["6MJP LptC LPS-export (81%, E=29)", "scattered DUF/peptidase/photosystem hits, no convergence"],
        "note": "No confident hit (best 81%, E=29, no convergence; the LptC LPS-transport top hit is biologically implausible as M. tuberculosis has no LPS). Remains dark.",
    },
    "Rv0150c": {
        "no_confident_hit": True,
        "top_hits": ["PF31149 TA_Actino tail-anchor membrane (33%, E=120)"],
        "note": "No confident hit (single hit at 33%, E=120) — remains dark.",
    },
    "Rv0395": {
        "no_confident_hit": True,
        "top_hits": ["PF15539 CAF1-p150 eukaryotic (73%, E=13)", "PF04613 LpxD acyltransferase (61%, E=41)", "no fold convergence, all E>>1"],
        "note": "No confident hit (best 73%, E=13, a eukaryotic chromatin factor — implausible; no convergence). Remains dark.",
    },
    "Rv2513": {
        "no_confident_hit": True,
        "top_hits": ["COG0048 RpsL/S12 ribosomal (40%, E=66)", "KOG1750 organellar S12 (25%)"],
        "note": "No confident hit (best 40%, E=66) — remains dark.",
    },
    "Rv2307B": {
        "no_confident_hit": True,
        "top_hits": ["PF15284 PAGK phage virulence factor (83%, E=3.3)", "3LAY zinc-resistance surface protein (75%, E=17)", "scattered surface/secreted hits, no fold convergence"],
        "note": "No confident hit (best 83%, E=3.3, no convergence). Scattered surface/secreted-protein hits are consistent with its predicted secretion (SP) but name no specific function — remains dark.",
    },
    "Rv1374c": {
        "no_confident_hit": True,
        "top_hits": ["PF01275 Myelin_PLP (44%, E=180)", "KOG4800 neuronal membrane glycoprotein (42%, E=66)", "2NC8 Mtb lipoprotein LppM (25%, E=110)"],
        "note": "No confident hit (best 44%, E=180). A weak, non-significant lean toward a membrane proteolipid fold (Myelin-PLP / LppM), not assertable — remains dark.",
    },
    "Rv0612": {
        "verdict": "family_assigned", "confidence": "medium",
        "function_revised": "Belongs to the DUF6155 family (Pfam PF19652, HHpred 96.6%, near full-length), with a converging hit to DUF6880 (PF21810, 96%) and to the nuclear proteasome-tether protein Cut8 (Pfam PF08559 / PDB 3Q5W, ~94%, E~4). A domain is assigned; the Cut8 match is a suggestive lead toward a proteasome-associated role (plausible given the Mtb Pup-proteasome), but at moderate E-value it is not asserted, and the DUF families are uncharacterised, so the precise function remains undefined. (The many ribosomal-S20 / TPR hits below ~85% with E>>1 are noise.)",
        "top_hits": ["PF19652 DUF6155 (96.6%, E=0.53)", "PF21810 DUF6880 (96%, E=0.41)", "PF08559/3Q5W Cut8 proteasome tether (94%, E~4)"],
    },
    "Rv2512a": {
        "no_confident_hit": True,
        "top_hits": ["2L0Z_A arenavirus envelope glycoprotein G2 zinc-binding domain (35.6%, E=0.6)", "PF22635 I-TevI_ZnF zinc-finger endonuclease domain (24.7%, E=0.4)", "6H8F TssA / 4LMQ SDF-1 / PF24044 DUF7353 / knottin (all 20-22%, E>=0.9, scattered)"],
        "note": "No confident hit (best 35.6%, an arenavirus envelope glycoprotein — implausible for a mycobacterial protein; E=0.6). The second, biologically apter hit (I-TevI zinc-finger, 24.7%) coincides with a cysteine-rich stretch in the sequence (C29-x-C31-...-C36 with flanking His27/His50), which is compatible with a possible metal/zinc-binding module; but at <25% probability, below the confident threshold and without fold convergence, this is a recorded lead, not an assignment. Remains dark.",
    },
    "Rv1116": {
        "no_confident_hit": True,
        "top_hits": ["PF19946 DUF6408 domain of unknown function (22.1%, E=1.3, ~res 18-42)"],
        "note": "Single sub-threshold hit to DUF6408 (22.1%, E=1.3), a domain of unknown function by definition, over ~24 of 61 residues. No convergence, probability well below the confident threshold. Remains dark.",
    },
    "Rv0063a": {
        "no_confident_hit": True,
        "top_hits": ["PF25184 YxzE putative bacteriocin (70.2%, E=2.1)", "PF14030 DUF4245 (69.4%, E=3.9)", "PF04971 Phage_holin_2_1 P21 holin (68.3%, E=1.5)", "PF30068 Daisho antimicrobial peptide (65%, E=1.8)", "converging small membrane/secreted-peptide + T3SS/flagellar export-gate hits (FliQ, SpaQ, MafB2 immunity, conotoxin, YhcB) all 43-58%, E>1"],
        "note": "No confident hit (best 70.2%, YxzE putative bacteriocin, E=2.1; all 25 hits are 43-70% with E-values >1). The hits converge thematically — NOT on a single named fold — toward small membrane/secreted peptides (bacteriocins, phage holins, antimicrobial peptides such as Daisho/MafB2-immunity/conotoxin) and T3SS/flagellar export-gate membrane proteins (FliQ, SpaQ), which is consistent with the DeepTMHMM-predicted signal peptide (SP) and the small size (53 aa). Points to a small secreted or membrane-embedded peptide, but no fold is assignable at the confident level. Remains dark; the small-secreted-peptide lean is recorded, not asserted.",
    },
    "Rv3430a": {
        "no_confident_hit": True,
        "top_hits": ["6ZPJ_B KKT4 kinetochore Leishmania mexicana (80.5%, E=4.3)", "6ZPM_B KKT4 kinetochore Trypanosoma cruzi (78%, E=4.3)", "5U59_A designed dimeric coiled-coil peptide (77.7%, E=3.4)", "centrosomin / GCN4 / kinesin coiled-coils (44-57%, E>4)"],
        "note": "No confident hit. All top hits are coiled-coil templates at 78-80.5% (below the 95% confident threshold) with E-values in the noise range (3.4-4.3), and are EXCLUSIVELY eukaryotic/synthetic (KKT4 kinetochore, designed peptide, centrosomin, GCN4, kinesin). This is the classic coiled-coil composition artefact: coiled-coils cross-match by heptad periodicity, not by homology. No Foldseek / AlphaFold-DB corroboration (struct_af empty), pLDDT only 70. The C-terminal sequence is plausibly α-helical, but a coiled-coil is not a function and cannot be assigned from these hits. Remains dark (guarded against over-calling; cf. KB coiled-coil-artefact rule 2026-06-10).",
    },
    "Rv1434": {
        "no_confident_hit": True,
        "top_hits": ["3O6Q_B Stage II sporulation protein SB / SpoIISB (25.6%, E=1.8)", "PF14185 SpoIISB antitoxin type II toxin-antitoxin (25.4%, E=1.9)", "7CGP_J Tim10 mitochondrial translocase (22.8%, E=2.4)", "SARS-CoV-2 M protein (2 hits, 22-22.4%, E>2)"],
        "note": "No confident hit (best 25.6%, SpoIISB antitoxin from Bacillus sporulation, E=1.8; all hits are 22-26%, E>1.5). The top hits do not converge on a single fold; the antitoxin SpoIISB is biologically plausible if Rv1434 were a TA component, but at <26% probability it is below the threshold for assignment. The remaining hits (Tim10 translocase, SARS-CoV-2 M protein) are implausible and scattered. Previously noted as a weak neighbour of hsr1 in the literature survey (P7.2), so the lack of HHpred signal is consistent. Remains dark; the tentative antitoxin lean is recorded but not asserted.",
    },
}


def _validate_hits(hits):
    """Garde-fou P7.1c : chaque entrée est SOIT un no-hit, SOIT une requalif bien formée.
    Empêche une entrée malformée (les deux à la fois, ou requalif sans function_revised/verdict/
    confidence, ou top_hits vide) de passer silencieusement lors d'un hand-off itératif."""
    errors = []
    for rv, h in hits.items():
        no_hit = h.get("no_confident_hit", False)
        has_verdict = "verdict" in h
        if no_hit and has_verdict:
            errors.append(f"{rv}: no_confident_hit ET verdict (ambigu)")
        if not no_hit and not has_verdict:
            errors.append(f"{rv}: ni no_confident_hit ni verdict")
        if not no_hit:
            if not h.get("function_revised"):
                errors.append(f"{rv}: requalif sans function_revised")
            if h.get("verdict") not in ("requalified", "family_assigned"):
                errors.append(f"{rv}: verdict inattendu {h.get('verdict')!r}")
            if not h.get("confidence"):
                errors.append(f"{rv}: requalif sans confidence")
        if not h.get("top_hits"):
            errors.append(f"{rv}: top_hits manquant/vide")
    if errors:
        raise SystemExit("phase34 : dict HITS invalide —\n  " + "\n  ".join(errors))


def main():
    print("== phase34 : requalification HHpred (P7.1) ==")
    _validate_hits(HITS)
    out = {rv: {**h, "references": [HHPRED_REF], "source": "HHpred profile-profile (MPI Toolkit web)"} for rv, h in HITS.items()}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "hhpred_curation.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))

    n = n_req = n_nohit = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        h = HITS.get(d["rv"])
        if not h:
            continue
        d["hhpred"] = {"top_hits": h["top_hits"], "source": "HHpred (MPI Toolkit)",
                       "no_confident_hit": h.get("no_confident_hit", False),
                       "note": h.get("note", "")}
        if h.get("no_confident_hit"):
            n_nohit += 1
            # un no-hit HHpred ne requalifie jamais : rétablir l'état dark canonique si une
            # passe antérieure avait sur-classé ce gène (rétrogradation propre et idempotente).
            # No-op pour les gènes déjà dark (verdict inchangé, function_revised préservé, y c. queue Foldseek).
            if d.get("verdict") != "dark":
                d["verdict"] = "dark"
                d["function_revised"] = "Conserved hypothetical protein; no recognised domain. Function unknown."
                d["confidence"] = "low"
                d["auto"] = True
                d["needs_review"] = False
                print(f"  {d['rv']} [no confident hit] RÉTROGRADÉ -> dark")
            else:
                print(f"  {d['rv']} [no confident hit] reste dark")
        else:
            d["function_revised"] = h["function_revised"]
            d["verdict"] = h["verdict"]
            d["confidence"] = h["confidence"]
            d["auto"] = False
            d["needs_review"] = False
            refs = d.get("references") or []
            if HHPRED_REF["doi"] not in {r.get("doi") for r in refs}:
                refs.append(HHPRED_REF)
            d["references"] = refs
            n_req += 1
            print(f"  {d['rv']} [{h['verdict']}/{h['confidence']}]: {h['function_revised'][:70]}")
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
    print(f"HHpred fusionné dans {n} fiche(s) : {n_req} requalifié(s), {n_nohit} sans hit (restent dark)")


if __name__ == "__main__":
    main()
