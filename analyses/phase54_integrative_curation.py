#!/usr/bin/env python3
"""phase54_integrative_curation.py -- P8 étape 4 : matérialisation des verdicts intégratifs par batch.

Curation ciblée au squelette des scripts existants (phase39/41/49) : dict CURATIONS keyé par Rv, écrit EN
PLACE dans site/content/genes/Rv*.json (NE PAS re-runner phase4, destructif). Politique conservatrice à 2
niveaux (plan P8, décidée avec CG) :
  - flip dark→family_assigned/low : SEULEMENT si ≥2 handles INDÉPENDANTS convergents ont survécu à la vérif
    adversariale → écrit verdict/confidence/function_revised + couche integrative_lead.
  - lead documenté sans flip : verdict=dark maintenu, function_revised INCHANGÉE, mais couche integrative_lead
    écrite (hypothèse tracée). Le garde-fou P6.2 d'ingest ne recale PAS (function_revised reste conclusive-dark).
Chaque entrée provient d'un batch APPROUVÉ par CG (synthèse + réfutation adversariale). auto=False.
Run: python analyses/phase54_integrative_curation.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
DATE = "2026-07-05"

# Rempli batch par batch APRÈS approbation CG. Schéma d'une entrée :
#   "RvXXXX": {
#       "flip": True/False,                     # True => verdict devient family_assigned ; False => reste dark
#       "confidence": "low"|"medium",           # utilisé seulement si flip
#       "function_revised": "candidate <voie>-associated component, ... (à valider)",  # seulement si flip
#       "integrative_lead": {
#           "hypothesis": "...",
#           "handles": ["operon: co-transcribed with X (pathway)", "STRING fusion+cooccurrence with Y", ...],
#           "n_independent_handles": 2,
#           "real_gene": "conserved (snp_sites=.. / 145k), MS-detected, ...",
#           "adversarial": "skeptic: not refuted — handles orthogonaux confirmés",
#           "source": f"integrative multi-layer synthesis (P8), {DATE}",
#       },
#       "refs": [ {authors, year, title, journal, doi}, ... ],   # optionnel
#   }
CONKLE = {"authors": "Conkle-Gutierrez D, et al.", "year": 2022,
          "title": "Mutations Associated with Capreomycin Resistance in Mycobacterium tuberculosis",
          "journal": "Antimicrob Agents Chemother", "doi": "10.1128/aac.02075-21"}

# Batch 1 (P8.1) — 5 HIGH, tous VRAIS gènes mais AUCUN flip (< 2 handles fonctionnels indépendants convergents).
# Leads documentés (versions corrigées par la vérification adversariale). Verdict dark maintenu.
CURATIONS: dict = {
    "Rv2680": {
        "flip": False,
        "integrative_lead": {
            "hypothesis": "candidate accessory/associated component of the RNA-metabolism operon Rv2680-Rv2681 (rnd / RNase D), to validate — no nuclease activity of its own",
            "handles": ["operon: co-transcribed with rnd/Rv2681 (RNase D, gap 1 bp) — single syntenic handle"],
            "n_independent_handles": 1,
            "real_gene": "conserved (snp_sites=1/145209, 0 nonsense/frameshift), MS-detected (15 datasets), globular cytoplasmic (AF pLDDT 86); DUF3000, no significant Foldseek hit, no M-CSA",
            "adversarial": "skeptic corrected an over-sell: removed the 'exoribonuclease' label (co-transcribed with rnd != is rnd; the RNase D is Rv2681, M-CSA 5/5). Citation fixed to Conkle-Gutierrez et al. AAC 2022 (PMID 35532237). SigH regulon, ethambutol/in-vivo phenotype and promoter->capreomycin-resistance are condition/expression, not function.",
        },
        "refs": [CONKLE],
    },
    "Rv0039c": {
        "flip": False,
        "integrative_lead": {
            "hypothesis": "candidate SigH-regulon-associated small integral membrane protein, plausibly in envelope/stress-response membrane physiology during infection, to validate",
            "handles": [],
            "n_independent_handles": 0,
            "real_gene": "conserved (snp_sites=4/145209, 4 missense, 0 nonsense/frameshift), 115 aa (>micro-ORF), 3 TM helices; MS-negative expected for a small hydrophobic membrane protein",
            "adversarial": "skeptic not refuted (dark confirmed, even reinforced). leuS STRING link = neighborhood/synteny AND aminoacyl-tRNA-synthetase hub artefact; mazG/hadC text-mining only. SigH regulon + in-vivo fitness = when, not what. No fold/domain handle.",
        },
    },
    "Rv0184": {
        "flip": False,
        "integrative_lead": {
            "hypothesis": "candidate redox-stress (WhiB4-regulon) associated component, possibly linked to the lprB/lprC lipoproteins, to validate",
            "handles": [],
            "n_independent_handles": 0,
            "real_gene": "conserved (snp_sites=16/145209, pN/pS=0.086), MS-detected (10 datasets), confident globular fold (pLDDT 87), 249 aa",
            "adversarial": "skeptic not refuted. Operon neighbours functionally heterogeneous (ytpA lysophospholipase / Rv0185 metallohydrolase / bglS beta-glucosidase) -> no coherent pathway; radA/recA coexpression == WhiB4 regulon context + stress hubs; HypF Foldseek non-significant (fold != enzyme); tandem DUF2786+DUF7168 = unknown function.",
        },
    },
    "Rv0381c": {
        "flip": False,
        "integrative_lead": {
            "hypothesis": "candidate cell-envelope/surface lipoprotein required for in-vivo fitness, genomically linked to but functionally distinct from pyrE, to validate",
            "handles": [],
            "n_independent_handles": 0,
            "real_gene": "conserved (snp_sites=9/145209, pN/pS=0.188 strong purifying), MS-detected (8 datasets), 302 aa; predicted secreted lipoprotein (signal peptide + lipobox)",
            "adversarial": "skeptic hardened the verdict: the secreted-lipoprotein localisation is incompatible with the cytoplasmic pyrimidine-biosynthesis role of operon partner pyrE, so the single syntenic pyrE handle is NOT a valid functional handle (0 functional handles). PTP Foldseek hits all non-significant, no conserved CX5R active-site motif (fold != enzyme). Strong in-vivo fitness defect = when, not what.",
        },
    },
    "Rv0448c": {
        "flip": False,
        "integrative_lead": {
            "hypothesis": "candidate SigK-regulon / ufaA1-lipid-neighbourhood associated component, to validate",
            "handles": [],
            "n_independent_handles": 1,
            "real_gene": "conserved (snp_sites=2/145209, 0 nonsense/frameshift, non-pseudogene), MS-detected, globular (pLDDT ESM 92 / AF 82.8), 221 aa; DUF1365",
            "adversarial": "skeptic not refuted. ufaA1 STRING fusion/cooccurrence downgraded as operon co-inheritance artefact (same fact as neighbourhood); operon heterogeneous (ufaA1 cyclopropane-lipid vs sigK/rskA sigma cassette) -> no single pathway; Foldseek #1 = engineered Anticalin (eukaryotic design artefact), all hits non-significant; no M-CSA.",
        },
    },
    # --- Batch 2 (P8.2) — 5 ortho-paires, tous VRAIS gènes mais 0 flip (opérons hétérogènes / partenaire externe = membre de cluster / voies divergentes) ---
    "Rv0514": {
        "flip": False,
        "integrative_lead": {
            "hypothesis": "candidate cell-envelope lipoprotein co-regulated with the heme cluster / kstR regulon, to validate — NOT a heme-pathway enzyme",
            "handles": [],
            "n_independent_handles": 0,
            "real_gene": "conserved (snp_sites=5/145209, no nonsense/frameshift), MS-detected (7 datasets), 99 aa, predicted envelope lipoprotein (signal peptide + lipobox)",
            "adversarial": "skeptic not refuted. hemA-D STRING = neighbourhood/synteny = same fact as operon; the only external partner aroE points to shikimate (divergent pathway, B2 fails); lipoprotein localisation incompatible with cytoplasmic tetrapyrrole enzymes; Foldseek Mediator (eukaryote) non-significant (coiled-coil noise).",
        },
    },
    "Rv3103c": {
        "flip": False,
        "integrative_lead": {
            "hypothesis": "candidate cell-division-associated membrane component (adjacent to the ftsEX module), to validate",
            "handles": [],
            "n_independent_handles": 0,
            "real_gene": "conserved (snp_sites=6/145209, no disruption), 145 aa, membrane (1 TM), MS-negative (expected); DUF-less 'proline-rich' legacy",
            "adversarial": "skeptic not refuted. 7-gene syntenic block functionally heterogeneous (smpB/prfB translation, ftsEX division, mscS, fprA); all high STRING = neighbourhood echo (same fact); external partners are PE_PGRS paralogy-hub artefacts or dark hypotheticals (nothing to transfer); no catalytic/structural evidence; Osman 2017 phage-secretome screen only, TnSeq contradicts its 'essential' prediction.",
        },
    },
    "Rv1025": {
        "flip": False,
        "integrative_lead": {
            "hypothesis": "candidate accessory component of the conserved eno-divIC-Rv1025-ppx2 gene block (no pathway assignable); ALSO a candidate drug target (essential + CRISPRi VI=-5.7) warranting a dedicated deep-dive",
            "handles": [],
            "n_independent_handles": 0,
            "real_gene": "REAL, essential (TnSeq), MS-detected (10 datasets), highly conserved (snp_sites=5/145209), well-folded (AF pLDDT 95), 155 aa, DUF501",
            "adversarial": "skeptic not refuted. operon eno/divIC/ppx2 functionally heterogeneous; whiB1/whiB2 STRING cooccurrence = trivial phylogenetic co-inheritance between essential genes (not a functional partnership); Foldseek non-significant, no M-CSA. Orthogonal applied value: druggable deep-dive target (see piste P8.2-deepdive).",
        },
    },
    "Rv0219": {
        "flip": False,
        "integrative_lead": {
            "hypothesis": "candidate accessory membrane component (anchor/permease/transporter) of a lipid-remodelling/degradation module, to validate — NOT a lipid enzyme itself",
            "handles": [],
            "n_independent_handles": 0,
            "real_gene": "conserved (snp_sites=6/145209, intact ORF, pLDDT_AF 84), 4-TM polytopic membrane protein, MS-negative (expected)",
            "adversarial": "skeptic not refuted. operon Rv0218-Rv0222 = lipid metabolism (lipC/echA1) + kstR/whiB3 regulon = coherent context, BUT the only external STRING partner (aao/galE3 cooccurrence) does not converge on the lipid pathway (B2 fails); 4-TM fold is non-enzymatic.",
        },
    },
    "Rv0477": {
        "flip": False,
        "integrative_lead": {
            "hypothesis": "candidate conserved secreted small protein whose cooccurrence profile (cutinases/ESX-3) suggests an envelope/surface association rather than the deoxyribose metabolism of its operon neighbour deoC, to validate",
            "handles": [],
            "n_independent_handles": 0,
            "real_gene": "REAL, strongly constrained (snp_sites=2/145209), 148 aa, MS-detected (9 datasets), predicted secreted (signal peptide), high pLDDT but no significant Foldseek; DUF2599",
            "adversarial": "skeptic REFUTED the proposed 2-handle claim: secreted localisation is incompatible with the co-transcribed cytoplasmic deoC (aldolase) -> the operon link is a syntenic passenger, not co-function; the cutinase/ESX cooccurrence is a weak phylogenetic-profile channel, insufficient alone. Rv0081 (hypoxia regulator) = condition, not function.",
        },
    },
}


def main():
    n_flip = n_lead = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        c = CURATIONS.get(d["rv"])
        if not c:
            continue
        lead = dict(c["integrative_lead"])
        lead.setdefault("source", f"integrative multi-layer synthesis (P8), {DATE}")
        d["integrative_lead"] = lead
        d["auto"] = False
        if c.get("flip"):
            d["verdict"] = "family_assigned"
            d["confidence"] = c.get("confidence", "low")
            d["function_revised"] = c["function_revised"]
            d["needs_review"] = False
            n_flip += 1
        else:
            n_lead += 1  # verdict/function inchangés, lead documenté seulement
        if c.get("refs"):
            refs = d.get("references") or []
            have = {r.get("doi") for r in refs}
            for r in c["refs"]:
                if r.get("doi") not in have:
                    refs.append(r)
            d["references"] = refs
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
    print(f"phase54 : {n_flip} requalifiés (flip dark→family_assigned), {n_lead} leads documentés (dark maintenu).")
    if not CURATIONS:
        print("(CURATIONS vide — remplir après approbation CG du batch.)")


if __name__ == "__main__":
    main()
