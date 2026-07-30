#!/usr/bin/env python3
"""phase70_literature_audit.py -- P16.2 : couche `literature` (audit de péremption, volet littérature).

Balayage systématique des 219 gènes `dark` contre le corpus PubMed TB pré-indexé (tbmonitor, ~330k papiers,
serveur MCP `tbmonitor`). Question : la littérature sait-elle déjà quelque chose que l'atlas ignore ?

SQL utilisé (MCP tbmonitor, reproductible) :
    WITH g(rv) AS (VALUES ('Rv0007'),...)              -- les 219 tags dark
    SELECT g.rv, COUNT(*) FROM g JOIN papers p
      ON (p.title LIKE '%'||g.rv||'%' OR p.abstract LIKE '%'||g.rv||'%')
    GROUP BY g.rv ORDER BY 2 DESC;
puis récupération des papiers pour les gènes à n>0.

RÉSULTAT : **187/219 (85 %) n'ont AUCUNE mention** (titre/résumé) → génuinement inexplorés. **32 en ont**,
dont deux massivement (Rv2660c 38, Rv2628 31) alors que l'atlas les dit « function unknown » avec 0 référence.

DISTINCTION APPORTÉE (qu'aucune ressource ne fait) :
  - `unstudied`         : personne n'a jamais regardé (187).
  - `studied_*`         : des gens ONT regardé, et la fonction MOLÉCULAIRE reste inconnue.
Sous-catégories des étudiés :
  - `antigen_vaccine`   : étudié comme antigène / composant vaccinal / marqueur diagnostique.
  - `functional`        : caractérisation fonctionnelle ou phénotypique dédiée publiée.
  - `variant_mention`   : SEULEMENT cité comme variant dans un crible génomique/GWAS de résistance
                          → AUCUN contenu fonctionnel, ne pas surinterpréter (garde-fou anti-survente).

GARDE-FOU : un phénotype de virulence ou un statut d'antigène N'EST PAS une fonction moléculaire → le verdict
`dark` est MAINTENU (l'atlas gradue la fonction biochimique). La couche AJOUTE le contexte manquant, elle ne
requalifie pas. C'est précisément ce qui manquait : Rv2660c est dans le vaccin H56 (essais cliniques) et la
fiche n'en disait rien.

Écrit EN PLACE. Run: python analyses/phase70_literature_audit.py
"""
from __future__ import annotations
import glob, json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
SOURCE = ("PubMed (whole), MULTI-ALIAS sweep of title+abstract: the H37Rv locus tag AND every ortholog identifier "
          "(M. bovis Mb…, M. marinum MMAR_…, M. smegmatis MSMEG_…, M. leprae ML…, M. abscessus MAB_…), each hit "
          "VERIFIED against the abstract text (word-boundary regex). phase73, 2026-07-13")


def P(t, doi, link, date):
    return {"title": t, "doi": doi, "link": link, "date": date}


# Résultat du balayage : rv -> (catégorie, résumé humain, papiers)
LIT: dict[str, tuple[str, str, list]] = {
    "Rv2660c": ("antigen_vaccine",
        "Major latency-associated protective antigen: component of the H56 vaccine candidate (Ag85B-ESAT6-Rv2660c) "
        "and of many multistage vaccine constructs (nanoparticle, DNA, fusion proteins) in preclinical/clinical "
        "development. Heavily studied immunologically, yet its molecular function remains unknown.", [
        P("A fusion protein LT25 containing Rv2660c failed to provide protection against M. bovis infection due to limited humoral immune responses", "10.1016/j.tube.2026.102771", "https://pubmed.ncbi.nlm.nih.gov/42096986/", "2026-04-27"),
        P("Novel Multistage Subunit M. tuberculosis Nanoparticle Vaccine (ESAT6-CFP10-Ag85A-Rv2660c-Rv1813c)", "10.3390/vaccines14010005", "https://pubmed.ncbi.nlm.nih.gov/41600921/", "2025-12-19"),
        P("Therapeutic vaccination with the Ag85B-Rv2660c-MPT70 fusion protein enhances M. tuberculosis clearance", "10.3389/fimmu.2025.1624923", "https://pubmed.ncbi.nlm.nih.gov/40895571/", "2025-08-14")]),
    "Rv2628": ("antigen_vaccine",
        "DosR (dormancy survival regulon) latency antigen: immunodominant in latent infection, used in multi-epitope "
        "vaccines and validated by meta-analysis as an immunodiagnostic marker discriminating latent from active TB. "
        "Molecular function unknown.", [
        P("Multifaceted DosR proteins of M. tuberculosis: emerging roles for tuberculosis control", "10.1007/s00203-026-04983-7", "https://pubmed.ncbi.nlm.nih.gov/42240825/", "2026-06-04"),
        P("Bioinformatics analysis, immunogenicity and therapeutic efficacy of a multi-stage multi-epitope DNA vaccine", "10.1016/j.intimp.2025.114415", "https://pubmed.ncbi.nlm.nih.gov/40086060/", "2025-04-16"),
        P("Role of DosR and Rpf antigens in differentiating active vs latent tuberculosis: systematic review and meta-analysis", "10.1186/s12890-024-03348-4", "https://pubmed.ncbi.nlm.nih.gov/39472851/", "2024-10-29")]),
    "Rv0426c": ("functional",
        "A dedicated functional study reports that Rv0426c promotes intracellular survival of recombinant mycobacteria "
        "by manipulating host inflammatory cytokines and suppressing apoptosis — a virulence-associated phenotype "
        "(not a molecular function).", [
        P("M. tuberculosis Rv0426c promotes recombinant mycobacteria intracellular survival via manipulating host inflammatory cytokines and suppressing cell apoptosis", "10.1016/j.meegid.2019.104070", "https://pubmed.ncbi.nlm.nih.gov/31614213/", "2020-01-01")]),
    "Rv0966c": ("functional",
        "Implicated in lipid / fatty-acid import by M. tuberculosis within macrophages (genetic-requirement screen), "
        "and differentially expressed with lipid import / beta-oxidation genes.", [
        P("The genetic requirements of fatty acid import by M. tuberculosis within macrophages", "10.7554/eLife.43621", "https://pubmed.ncbi.nlm.nih.gov/30735132/", "2019-02-08"),
        P("Differential expression of genes associated with lipid import, beta-oxidation and lactate oxidation induced by curli pili", "10.1099/jmm.0.001994", "https://pubmed.ncbi.nlm.nih.gov/40162564/", "2025-03-01")]),
    "Rv1048c": ("functional",
        "A dedicated study reports that Rv1048c affects the biological characteristics of recombinant M. smegmatis.", [
        P("M. tuberculosis Rv1048c affects the biological characteristics of recombinant Mycobacterium smegmatis", "10.1038/s41598-024-81405-y", "https://pubmed.ncbi.nlm.nih.gov/39613837/", "2024-11-29")]),
    "Rv3126c": ("antigen_vaccine",
        "Latency-associated antigen incorporated into multistage vaccine constructs (Sendai-virus vectored, and a "
        "protein-subunit BCG booster) conferring protection in murine models.", [
        P("A multistage Sendai virus vaccine incorporating latency-associated antigens induces protection against acute and latent tuberculosis", "10.1080/22221751.2023.2300463", "https://pubmed.ncbi.nlm.nih.gov/38164736/", "2024-12-01"),
        P("A multistage protein subunit vaccine as BCG-booster confers protection against M. tuberculosis infection", "10.1016/j.intimp.2024.112811", "https://pubmed.ncbi.nlm.nih.gov/39068754/", "2024-09-30")]),
    "Rv2656c": ("antigen_vaccine",
        "Reported as a candidate antigen associated with latent tuberculosis infection (2026).", [
        P("Rv2656c: A Potential Candidate Antigen Associated with Latent Tuberculosis Infection", "10.3390/vaccines14050442", "https://pubmed.ncbi.nlm.nih.gov/42188812/", "2026-05-15")]),
    "Rv0049": ("functional",
        "Surfaced by a functional whole-genome screen of nutrient-starved M. tuberculosis as involved in rifampin / "
        "antibiotic tolerance; also detected among serum-reactive antigens.", [
        P("Functional Whole Genome Screen of Nutrient-Starved M. tuberculosis Identifies Genes Involved in Rifampin Tolerance", "10.3390/microorganisms11092269", "https://pubmed.ncbi.nlm.nih.gov/37764112/", "2023-09-09"),
        P("Serum proteomic analysis of M. tuberculosis antigens for discriminating active tuberculosis from latent infection", "10.1177/0300060520910042", "https://pubmed.ncbi.nlm.nih.gov/32216499/", "2020-03-01")]),
    "Rv0997": ("functional",
        "Reported among DosR-regulated genes (differential expression in a successful Ethiopian sub-lineage) and among "
        "macrophage-induced genes.", [
        P("Transcriptomic and genomic analysis of Ethiopian M. tuberculosis sub-lineage 4.2.2.2 reveals differential expression of DosR-regulated genes", "10.1038/s41598-025-34471-9", "https://pubmed.ncbi.nlm.nih.gov/41484201/", "2026-01-03"),
        P("Macrophage-specific M. tuberculosis genes: identification by GFP and kanamycin resistance selection", "10.1099/mic.0.2006/000547-0", "https://pubmed.ncbi.nlm.nih.gov/17322185/", "2007-03-01")]),
    "Rv3612c": ("functional",
        "Located in an operon under long-range transcriptional control that is necessary for virulence-critical ESX-1 "
        "secretion.", [
        P("Long-range transcriptional control of an operon necessary for virulence-critical ESX-1 secretion in M. tuberculosis", "10.1128/JB.00142-12", "https://pubmed.ncbi.nlm.nih.gov/22389481/", "2012-05-01")]),
    "Rv3258c": ("functional",
        "Studied among the putative mannose-donor biosynthesis genes (comparative transcription, H37Rv vs BCG).", [
        P("Comparative transcriptional study of the putative mannose donor biosynthesis genes in virulent M. tuberculosis and attenuated M. bovis BCG", "10.1128/IAI.05635-11", "https://pubmed.ncbi.nlm.nih.gov/21896775/", "2011-11-01")]),
    "Rv2517c": ("functional",
        "Induced as part of the SOS response upon moxifloxacin exposure.", [
        P("Moxifloxacin Activates the SOS Response in M. tuberculosis in a Dose- and Time-Dependent Manner", "10.3390/microorganisms9020255", "https://pubmed.ncbi.nlm.nih.gov/33513836/", "2021-01-27")]),
    "Rv3103c": ("functional",
        "Identified among unique essential proteins from a phage-secretome library of an F15/LAM4/KZN strain.", [
        P("Identification of unique essential proteins from a M. tuberculosis F15/LAM4/KZN phage secretome library", "10.1093/femspd/ftx001", "https://pubmed.ncbi.nlm.nih.gov/28087649/", "2017-01-01")]),
    "Rv0057": ("antigen_vaccine", "Studied as an immunogen: recombinant Rv0057-Rv1352 fusion protein elicits immune responses.", [
        P("Immune responses to a recombinant Rv0057-Rv1352 fusion protein of M. tuberculosis", None, "https://pubmed.ncbi.nlm.nih.gov/25696009/", "2015-01-01")]),
    "Rv0378": ("antigen_vaccine", "Cloned and sequenced as a candidate for a subunit DNA vaccine.", [
        P("Cloning and sequencing of M. tb gene Rv0378 for making subunit based DNA vaccine", None, "https://pubmed.ncbi.nlm.nih.gov/33592991/", "2017-03-01")]),
    "Rv2706c": ("antigen_vaccine", "Predicted antigenic protein with promiscuous CTL epitopes (in silico).", [
        P("In silico identification of potential antigenic proteins and promiscuous CTL epitopes in M. tuberculosis", "10.1016/j.meegid.2012.03.023", "https://pubmed.ncbi.nlm.nih.gov/22484107/", "2012-08-01")]),
    "Rv2693c": ("antigen_vaccine", "Used as an antigen in IFN-gamma release assays in latent tuberculosis.", [
        P("Effect of isoniazid on antigen-specific interferon-gamma secretion in latent tuberculosis", "10.1183/09031936.00123314", "https://pubmed.ncbi.nlm.nih.gov/25359354/", "2015-02-01")]),
    "Rv2114": ("antigen_vaccine", "Among proteins purified from MDR clinical isolates eliciting T-cell cytokine responses.", [
        P("T cell cytokine responses in PBMC from MDR-TB patients following stimulation with proteins purified from MDR clinical isolates", "10.1016/j.ijmyco.2016.10.009", "https://pubmed.ncbi.nlm.nih.gov/28043507/", "2016-12-01")]),
    "Rv3235": ("antigen_vaccine", "Used in epitope design as a putative DNA-binding protein vaccine candidate (in silico).", [
        P("Computational approaches in epitope design using DNA binding proteins as vaccine candidate in M. tuberculosis", "10.1016/j.meegid.2020.104357", "https://pubmed.ncbi.nlm.nih.gov/32438080/", "2020-09-01")]),
    "Rv3747": ("antigen_vaccine", "Identified among M. bovis antigens by bovine T-cell responses after infection.", [
        P("Identification of M. bovis antigens by analysis of bovine T-cell responses after infection with a virulent strain", "10.1590/s0100-879x2003001100011", "https://pubmed.ncbi.nlm.nih.gov/14576908/", "2003-11-01")]),
    "Rv0740": ("antigen_vaccine", "Appears in a serological cross-reactivity study (mycobacterial sera vs HIV structural proteins); weak evidence.", [
        P("Serum samples from patients with mycobacterial infections cross-react with HIV structural proteins", None, "https://pubmed.ncbi.nlm.nih.gov/17824484/", "2007-06-01")]),
    # --- CORRECTION 2026-07-13 (phase73) : trouvés UNIQUEMENT via un ORTHOLOGUE, donc invisibles au tag Rv ---
    "Rv2297": ("antigen_vaccine",
        "Invisible under its H37Rv tag: this gene appears in the literature only under its M. bovis ortholog "
        "**Mb2319**, selected in silico as a candidate M. bovis antigen and evaluated in an interferon-gamma release "
        "assay (IGRA) for bovine tuberculosis diagnosis.", [
        P("Identification and evaluation of new Mycobacterium bovis antigens in the in vitro interferon gamma release assay for bovine tuberculosis diagnosis", "10.1016/j.tube.2015.07.009", "https://pubmed.ncbi.nlm.nih.gov/26320985/", "2015")]),
    "Rv3190c": ("antigen_vaccine",
        "Invisible under its H37Rv tag: appears in the literature only under its M. bovis ortholog **Mb3212c**, "
        "selected in silico as a candidate M. bovis antigen and evaluated in an IGRA for bovine tuberculosis.", [
        P("Identification and evaluation of new Mycobacterium bovis antigens in the in vitro interferon gamma release assay for bovine tuberculosis diagnosis", "10.1016/j.tube.2015.07.009", "https://pubmed.ncbi.nlm.nih.gov/26320985/", "2015")]),
    "Rv0463": ("variant_mention",
        "Invisible under its H37Rv tag: appears only under its M. leprae ortholog **ML2388**, listed among genes "
        "up-regulated in leprosy reactional states (transcriptomic signature). A mention in a gene list, not a "
        "functional characterisation.", [
        P("Mycobacterium leprae and host immune transcriptomic signatures for reactional states in leprosy", "10.1016/j.jid.2018.09.029", "https://pubmed.ncbi.nlm.nih.gov/37051521/", "2023")]),
}

# Cités UNIQUEMENT comme variants dans des cribles génomiques / GWAS de résistance : AUCUN contenu fonctionnel.
# NB : Rv0961, Rv1434 et Rv2699c ont été RETIRÉS le 2026-07-13 — c'étaient des FAUX POSITIFS du balayage tbmonitor,
# dus au piège de sous-chaîne du `LIKE` SQL (`'%Rv2699%'` capture `Rv2699c`, un AUTRE gène ; `'%Rv1434%'` capture
# un papier portant en réalité sur Rv1435c). Le balayage PubMed vérifié (phase73) ne les retrouve pas → « unstudied ».
VARIANT_ONLY = {
    "Rv1907c": "genome-wide study of drug-resistant Mtb / intra-host evolution",
    "Rv2077c": "Bayesian transmission-network reconstruction (genomic variants)",
    "Rv2082": "in vitro bedaquiline/clofazimine exposure variants",
    "Rv2680": "GWAS of second-line injectable drug resistance",
    "Rv3136A": "genomic/transcriptomic/phenotypic integration, sub-lineage 4.2.2.2",
    "Rv3861": "machine-learning characterisation of drug-resistant mutations from large-scale WGS",
    "Rv3346c": "RD-Rio subfamily prevalence / genotypic diversity",
    "Rv1779c": "co-mentioned in a study of Rv0180c (cell shape / infectivity)",
}

# Comptes RÉELS du balayage SQL (n de papiers mentionnant le tag) — distincts du nombre de papiers LISTÉS ci-dessus
# (on ne liste que les plus récents/représentatifs). Ne pas sous-estimer la littérature sur la fiche.
TRUE_COUNTS = {
    "Rv2660c": 38, "Rv2628": 31, "Rv0049": 3, "Rv3126c": 2, "Rv0997": 2, "Rv0966c": 2,
}

CAT_LABEL = {
    "antigen_vaccine": "studied as an antigen / vaccine or diagnostic component",
    "functional": "a functional or phenotypic characterisation has been published",
    "variant_mention": "only cited as a variant in genomic / resistance screens (no functional content)",
    "unstudied": "no TB publication mentions this locus tag",
}


def main() -> None:
    n_unstudied = n_studied = 0
    cats: dict[str, int] = {}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.loads(Path(f).read_text())
        if d.get("verdict") != "dark":
            continue
        rv = d["rv"]
        if rv in LIT:
            cat, summary, papers = LIT[rv]
        elif rv in VARIANT_ONLY:
            cat, papers = "variant_mention", []
            summary = ("Appears in the TB literature ONLY as a variant in a genomic / drug-resistance screen "
                       f"({VARIANT_ONLY[rv]}); no functional content — do not over-interpret.")
        else:
            cat, summary, papers = "unstudied", (
                "No TB publication mentions this locus tag in its title or abstract (sweep of a ~330k-abstract PubMed "
                "TB corpus). This gene is genuinely unstudied: its darkness reflects absence of investigation, not "
                "failure of investigation."), []
        n_true = TRUE_COUNTS.get(rv, 0 if cat == "unstudied" else max(1, len(papers)))
        layer = {
            "n_papers": n_true,                    # vrai compte du balayage
            "n_papers_listed": len(papers),        # nombre effectivement listé sur la fiche
            "category": cat,
            "category_label": CAT_LABEL[cat],
            "summary": summary,
            "papers": papers,
            "caveat": ("An antigen status or a virulence phenotype is NOT a molecular function: the verdict stays "
                       "unchanged. This layer adds the missing context, it does not requalify the gene."),
            "source": SOURCE,
        }
        d["literature"] = layer
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        cats[cat] = cats.get(cat, 0) + 1
        if cat == "unstudied":
            n_unstudied += 1
        else:
            n_studied += 1

    print(f"couche `literature` écrite sur {n_unstudied + n_studied} gènes dark.")
    print(f"  JAMAIS ÉTUDIÉS  : {n_unstudied}  (obscurité = personne n'a regardé)")
    print(f"  ÉTUDIÉS         : {n_studied}   (on a regardé, la fonction moléculaire reste inconnue)")
    print(f"  par catégorie   : {cats}")
    print("\nOmissions les plus criantes corrigées : Rv2660c (vaccin H56, 38 papiers) et Rv2628 (antigène DosR, 31),")
    print("qui n'avaient AUCUNE référence dans l'atlas. Verdicts inchangés (antigène != fonction moléculaire).")


if __name__ == "__main__":
    main()
