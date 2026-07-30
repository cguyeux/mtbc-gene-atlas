from __future__ import annotations

from datetime import datetime

from sqlalchemy import JSON, Boolean, DateTime, Integer, String, Text
from sqlalchemy.orm import Mapped, mapped_column

from .db import Base


class Gene(Base):
    """One MTBC gene, anchored on the ancestral MTBC0 protein and joined to H37Rv.

    The catalogue (all joined genes) is loaded from gene_xref.tsv; enriched
    genes additionally carry a curated function, a verdict, Pfam domains,
    optional structural hits (Foldseek), ESM signal, evidence and references.
    Rich nested data is kept as JSON columns: the site filters on the scalar
    columns (rv, gene_name, products, verdict, is_hypothetical, enriched).
    """
    __tablename__ = "genes"

    rv: Mapped[str] = mapped_column(String(16), primary_key=True)
    mtbc0: Mapped[str | None] = mapped_column(String(24), index=True)
    np_id: Mapped[str | None] = mapped_column(String(24))
    gene_name: Mapped[str | None] = mapped_column(String(32), index=True)
    len_aa: Mapped[int | None] = mapped_column(Integer)
    strand: Mapped[str | None] = mapped_column(String(1))
    start_mtbc0: Mapped[int | None] = mapped_column(Integer)
    end_mtbc0: Mapped[int | None] = mapped_column(Integer)

    product_h37rv: Mapped[str | None] = mapped_column(Text)        # legacy
    product_pgap: Mapped[str | None] = mapped_column(Text)         # PGAP re-annotation
    function_revised: Mapped[str | None] = mapped_column(Text)     # curated
    verdict: Mapped[str | None] = mapped_column(String(24), index=True)  # requalified|family_assigned|dark
    confidence: Mapped[str | None] = mapped_column(String(16))
    curation_note: Mapped[str | None] = mapped_column(Text)  # provenance of a manual curation/adjudication (e.g. P16.14)

    is_hypothetical: Mapped[bool] = mapped_column(Boolean, default=False, index=True)
    enriched: Mapped[bool] = mapped_column(Boolean, default=False, index=True)
    auto: Mapped[bool] = mapped_column(Boolean, default=False, index=True)         # rule-based curation
    needs_review: Mapped[bool] = mapped_column(Boolean, default=False, index=True)  # flagged for a human pass

    protein_mtbc0: Mapped[str | None] = mapped_column(Text)
    esm: Mapped[dict | None] = mapped_column(JSON)
    domains: Mapped[list | None] = mapped_column(JSON)
    struct_hits: Mapped[list | None] = mapped_column(JSON)
    eggnog: Mapped[dict | None] = mapped_column(JSON)  # orthology: COG/EC/KEGG/GO/CAZy
    uniprot: Mapped[dict | None] = mapped_column(JSON)  # curated: SwissProt function/EC/PE level
    conservation: Mapped[dict | None] = mapped_column(JSON)  # intra-MTBC pN/pS + disruption
    plddt: Mapped[dict | None] = mapped_column(JSON)    # ESMFold model confidence (gates Foldseek)
    struct_af: Mapped[dict | None] = mapped_column(JSON)  # genome-wide AlphaFold-model Foldseek hits (P3.2)
    plddt_af: Mapped[dict | None] = mapped_column(JSON)   # AlphaFold DB model confidence (mean pLDDT)
    string: Mapped[dict | None] = mapped_column(JSON)   # STRING v12 functional interaction network
    mcsa: Mapped[dict | None] = mapped_column(JSON)     # M-CSA catalytic-site verification (P3.1)
    pfam_tentative: Mapped[dict | None] = mapped_column(JSON)  # sub-threshold Pfam lead for dark genes (P3.3)
    localization: Mapped[dict | None] = mapped_column(JSON)  # predicted localisation (Kyte-Doolittle TM + lipobox, P3.5)
    essentiality: Mapped[dict | None] = mapped_column(JSON)  # Tn-seq essentiality (DeJesus 2017 / Griffin 2011)
    proteomics: Mapped[dict | None] = mapped_column(JSON)  # MS existence proof + abundance (PaxDb 5.0, P5.1)
    protparam: Mapped[dict | None] = mapped_column(JSON)  # physico-chemical properties (Biopython ProtParam, P5.3)
    resistance: Mapped[dict | None] = mapped_column(JSON)  # WHO drug-resistance association per gene (P5.9)
    funccat: Mapped[dict | None] = mapped_column(JSON)  # TubercuList functional category (P5.4)
    orthologs: Mapped[dict | None] = mapped_column(JSON)  # species-named orthologs, RBH DIAMOND (P5.5)
    mutant_phenotypes: Mapped[dict | None] = mapped_column(JSON)  # conditional Tn-seq phenotypes, MtbTnDB (P5.2)
    pdb: Mapped[dict | None] = mapped_column(JSON)  # experimental PDB structures via PDBe/SIFTS (P5.8)
    genomic_context: Mapped[dict | None] = mapped_column(JSON)  # neighbors + predicted operon (P5.7)
    regulation: Mapped[dict | None] = mapped_column(JSON)  # TF regulators + regulon, signed TRN (P5.6)
    mycobrowser: Mapped[dict | None] = mapped_column(JSON)  # legacy record + concordance vs Mycobrowser (P5.11)
    rd: Mapped[dict | None] = mapped_column(JSON)  # Regions of Difference: per-gene deletion across lineages (P5.10)
    ec_override: Mapped[dict | None] = mapped_column(JSON)  # curated EC override (eggNOG->UniProt) (P5.11c)
    phenotype_lead: Mapped[dict | None] = mapped_column(JSON)  # phenotype-driven functional hypothesis (P5.2b)
    hhpred: Mapped[dict | None] = mapped_column(JSON)  # HHpred profile-profile top hits / no-hit record (P7.1)
    struct_cluster: Mapped[dict | None] = mapped_column(JSON)  # all-vs-all structural cluster / guilt-by-structure (P7.3)
    context_lead: Mapped[dict | None] = mapped_column(JSON)  # operon-context functional lead (co-transcription) (P7.10)
    vulnerability: Mapped[dict | None] = mapped_column(JSON)  # CRISPRi vulnerability index, Bosch 2021 (P7.6)
    ptm: Mapped[dict | None] = mapped_column(JSON)  # post-translational modifications (phospho/PTM), UniProt (P7.8)
    expression: Mapped[dict | None] = mapped_column(JSON)  # conditional co-expression context: iModulons, Yoo 2022 (P7.9)
    integrative_lead: Mapped[dict | None] = mapped_column(JSON)  # multi-layer integrative synthesis for dark genes (P8)
    outgroup: Mapped[dict | None] = mapped_column(JSON)  # M. canettii dN/dS + NTM presence + cross-genus phylostratum (P10)
    rbtnseq_fitness: Mapped[dict | None] = mapped_column(JSON)  # RB-TnSeq 95-condition fitness phenotypes, PLoS Biol 2026 (P16.9)
    literature: Mapped[dict | None] = mapped_column(JSON)  # TB-literature sweep: studied-vs-unstudied + antigen/vaccine context (P16.2)
    disorder: Mapped[dict | None] = mapped_column(JSON)  # intrinsic disorder (metapredict + pLDDT); flags IDR-orphan dark genes (P16.13)
    prospect: Mapped[dict | None] = mapped_column(JSON)  # PROSPECT chemical-genetic: hypomorph tool strain + knockdown fitness + druggability (Nat Commun 2025) (P16.10)
    abpp: Mapped[dict | None] = mapped_column(JSON)  # activity-based protein profiling: FP-reactive active serine hydrolase (Cell Chem Biol 2021) (P16.12)
    evidence: Mapped[list | None] = mapped_column(JSON)
    refs: Mapped[list | None] = mapped_column(JSON)


class Microprotein(Base):
    """A non-canonical, MS-proven microprotein (smORF) excluded from the H37Rv annotation.

    Separate track (P16.5): proteogenomic microproteome of de Souza et al. 2024. These
    are NOT canonical genes and are counted apart from the 3906 CDS and the dark stock.
    overlap_class is classified by GEOMETRY vs the H37Rv CDS set (intergenic / antisense /
    sense_overlap), never trusting the source's own tORF/gORF label. Existence is proven
    (high-confidence peptide-spectrum match); function is unknown.
    """
    __tablename__ = "microproteins"

    id: Mapped[str] = mapped_column(String(32), primary_key=True)   # e.g. tORF_35175
    kind: Mapped[str | None] = mapped_column(String(24))            # source label (informational only)
    overlap_class: Mapped[str | None] = mapped_column(String(16), index=True)  # intergenic|antisense|sense_overlap
    overlap_locus: Mapped[str | None] = mapped_column(String(16))   # H37Rv Rv it overlaps (geometry)
    tier: Mapped[str | None] = mapped_column(String(4))             # T1|T2 (MS high-confidence only)
    ms_evidence: Mapped[str | None] = mapped_column(Text)           # existence-evidence string (any method)
    detection_method: Mapped[str | None] = mapped_column(String(64), index=True)  # MS compendium | VapC4 ribosome-stalling | ...
    extra: Mapped[dict | None] = mapped_column(JSON)                # source-specific fields (e.g. cysteines, cross-validation)
    start: Mapped[int | None] = mapped_column(Integer)
    end: Mapped[int | None] = mapped_column(Integer)
    strand: Mapped[str | None] = mapped_column(String(1))
    length_aa: Mapped[int | None] = mapped_column(Integer, index=True)
    start_codon: Mapped[str | None] = mapped_column(String(4))
    shine_dalgarno: Mapped[str | None] = mapped_column(Text)
    folding_free_energy: Mapped[float | None] = mapped_column()
    sequence: Mapped[str | None] = mapped_column(Text)
    localization: Mapped[str | None] = mapped_column(String(32))
    localization_reliability: Mapped[str | None] = mapped_column(String(32))
    essentiality: Mapped[str | None] = mapped_column(String(48))
    essentiality_code: Mapped[str | None] = mapped_column(String(4), index=True)  # ES|GD|GA
    conserved_in: Mapped[list | None] = mapped_column(JSON)
    host_gene: Mapped[str | None] = mapped_column(String(32))       # source-reported host gene name
    host_rv: Mapped[str | None] = mapped_column(String(16), index=True)
    surrounding_genes: Mapped[str | None] = mapped_column(Text)
    source: Mapped[str | None] = mapped_column(Text)


class Feedback(Base):
    __tablename__ = "feedback"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    gene_rv: Mapped[str | None] = mapped_column(String(16), index=True)
    email: Mapped[str | None] = mapped_column(String(255))
    message: Mapped[str] = mapped_column(Text)
    status: Mapped[str] = mapped_column(String(16), default="new")  # new|reviewed|integrated|rejected
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow, index=True)
