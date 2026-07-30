"""Public REST API for the MTBC gene-annotation atlas.

Mounted as a self-contained FastAPI sub-application at /api/v1, so it ships its own
OpenAPI schema (/api/v1/openapi.json) and interactive docs (/api/v1/docs) — a stable,
citable machine interface to the 3906 enriched MTBC0/H37Rv gene records and their
~35 evidence layers. Read-only, CORS-open for GET (usable from any origin / notebook).
"""
from __future__ import annotations

from fastapi import Depends, FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from sqlalchemy import and_, func, or_, select
from sqlalchemy.orm import Session

from .db import get_db
from .models import Gene, Microprotein

API_VERSION = "1.0"
CITATION = ("Guyeux C. MTBC Gene Annotation Atlas: a re-annotation of the "
            "Mycobacterium tuberculosis complex proteome. https://mtbc.gclab.fr "
            "(REST API https://mtbc.gclab.fr/api/v1).")

# Scalar columns returned in full records / summaries.
SCALAR_FIELDS = (
    "rv", "mtbc0", "np_id", "gene_name", "len_aa", "strand", "start_mtbc0", "end_mtbc0",
    "product_h37rv", "product_pgap", "function_revised", "verdict", "confidence",
    "is_hypothetical", "enriched", "auto", "needs_review",
)

# JSON evidence layers (name -> one-line description). Order = pipeline order.
LAYERS: dict[str, str] = {
    "esm": "ESM Atlas signal (exploratory SAE features / ESMFold)",
    "domains": "Pfam domains (hmmscan --cut_ga)",
    "struct_hits": "Foldseek hits on the ESMFold model (dark genes)",
    "eggnog": "eggNOG orthology: COG / EC / KEGG KO / pathway / GO / CAZy",
    "uniprot": "UniProt/SwissProt curated function, EC, protein-existence level",
    "conservation": "intra-MTBC selection: pN/pS + disruption load over 145k strains",
    "plddt": "ESMFold model confidence (gates Foldseek)",
    "struct_af": "genome-wide AlphaFold-model Foldseek hits vs AFDB",
    "plddt_af": "AlphaFold DB model confidence (mean pLDDT)",
    "string": "STRING v12 functional interaction network (guilt-by-association)",
    "mcsa": "M-CSA catalytic-site verification (fold vs active enzyme)",
    "pfam_tentative": "sub-threshold Pfam lead for dark genes",
    "localization": "predicted subcellular localisation (DeepTMHMM: TM/SP/lipoprotein)",
    "essentiality": "Tn-seq essentiality (DeJesus 2017 / Griffin 2011)",
    "proteomics": "MS existence proof + abundance (PaxDb 5.0)",
    "protparam": "physico-chemical properties (Biopython ProtParam)",
    "resistance": "WHO drug-resistance association per gene",
    "funccat": "TubercuList functional category",
    "orthologs": "species-named orthologs (RBH DIAMOND)",
    "mutant_phenotypes": "conditional Tn-seq phenotypes (MtbTnDB)",
    "pdb": "experimental PDB structures (PDBe/SIFTS)",
    "genomic_context": "genomic neighbours + predicted operon",
    "regulation": "transcription-factor regulators + regulon (signed TRN)",
    "mycobrowser": "legacy Mycobrowser record + field-by-field concordance",
    "rd": "Regions of Difference: per-gene deletion across lineages",
    "ec_override": "curated EC override (eggNOG -> UniProt)",
    "phenotype_lead": "phenotype-driven functional hypothesis",
    "hhpred": "HHpred profile-profile top hits / no-hit record",
    "struct_cluster": "all-vs-all structural cluster (guilt-by-structure)",
    "context_lead": "operon-context functional lead (co-transcription)",
    "vulnerability": "CRISPRi vulnerability index (Bosch 2021)",
    "ptm": "post-translational modifications / phosphosites (UniProt)",
    "expression": "conditional co-expression context: iModulons (Yoo 2022)",
    "integrative_lead": "multi-layer integrative synthesis for dark genes",
    "outgroup": "deep-divergence conservation: M. canettii dN/dS vs MTBC; genus-wide NTM presence; cross-genus phylostratum / gene age (P10)",
    "rbtnseq_fitness": "RB-TnSeq 95-condition conditional fitness phenotypes (PLoS Biol 2026); experimental, condition-specific (P16.9)",
    "literature": "TB-literature sweep (~330k PubMed abstracts): studied vs genuinely unstudied, antigen/vaccine and functional context (P16.2)",
    "disorder": "intrinsic disorder (metapredict + AlphaFold pLDDT); flags dark genes that are dark-by-construction (IDR-orphan) (P16.13)",
    "prospect": "PROSPECT chemical-genetic target-ID resource: TetON hypomorph tool strain, baseline knockdown fitness, and drug-target/MOA cross-reference (Bond 2025) (P16.10)",
    "abpp": "activity-based protein profiling: biochemically confirmed FP-reactive active serine hydrolase (census + prioritized targets), probe enrichment / covalent-inhibitor competition, and concordance with the atlas family call (Li 2021) (P16.12)",
    "evidence": "evidence trail for the revised function",
    "refs": "cited references",
}

VERDICTS = ("requalified", "family_assigned", "dark")

api = FastAPI(
    title="MTBC Gene Annotation Atlas API",
    version=API_VERSION,
    description=(
        "Read-only REST access to a re-annotation of the *Mycobacterium tuberculosis* "
        "complex proteome (3906 genes, anchored on the ancestral MTBC0 genome and joined "
        "to H37Rv NC_000962.3). Each gene carries a curated function, a verdict "
        "(`requalified` / `family_assigned` / `dark`) and up to ~35 evidence layers "
        "(orthology, structure, conservation, interaction network, essentiality, "
        "CRISPRi vulnerability, expression, phenotypes, ...). "
        f"Cite as: {CITATION}"
    ),
    docs_url="/docs", redoc_url="/redoc",
)
api.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["GET"], allow_headers=["*"])


def _summary(g: Gene) -> dict:
    return {
        "rv": g.rv, "gene_name": g.gene_name, "mtbc0": g.mtbc0,
        "product": g.product_pgap or g.product_h37rv,
        "function_revised": g.function_revised,
        "verdict": g.verdict, "confidence": g.confidence,
        "is_hypothetical": g.is_hypothetical, "enriched": g.enriched,
        "len_aa": g.len_aa,
        "url": f"/gene/{g.rv}",
    }


def _full(g: Gene) -> dict:
    rec = {k: getattr(g, k) for k in SCALAR_FIELDS}
    rec["protein_mtbc0"] = g.protein_mtbc0
    rec["layers"] = {name: getattr(g, name) for name in LAYERS if getattr(g, name, None)}
    return rec


@api.get("/", summary="API root: version, dataset size, citation, endpoint index")
def root(db: Session = Depends(get_db)) -> dict:
    total = db.scalar(select(func.count()).select_from(Gene)) or 0
    return {
        "name": "MTBC Gene Annotation Atlas API",
        "version": API_VERSION,
        "genes": total,
        "reference": "MTBC0 (ancestral) joined to H37Rv NC_000962.3",
        "citation": CITATION,
        "website": "https://mtbc.gclab.fr",
        "docs": "/api/v1/docs",
        "openapi": "/api/v1/openapi.json",
        "license": "CC-BY 4.0",
        "endpoints": {
            "GET /api/v1/stats": "dataset counts (verdicts, layer coverage)",
            "GET /api/v1/layers": "list of the ~35 evidence layers + descriptions",
            "GET /api/v1/genes": "search / filter genes (paginated summaries)",
            "GET /api/v1/genes/{rv}": "full record for one gene (all layers)",
            "GET /api/v1/genes/{rv}/{layer}": "one evidence layer for one gene",
        },
    }


@api.get("/stats", summary="Dataset counts: verdicts and per-layer coverage")
def stats(db: Session = Depends(get_db)) -> dict:
    total = db.scalar(select(func.count()).select_from(Gene)) or 0
    by_verdict = dict(
        db.execute(select(Gene.verdict, func.count())
                   .where(Gene.verdict.is_not(None)).group_by(Gene.verdict)).all()
    )
    hypo = db.scalar(select(func.count()).where(Gene.is_hypothetical.is_(True))) or 0
    coverage = {
        name: (db.scalar(select(func.count()).where(getattr(Gene, name).is_not(None))) or 0)
        for name in LAYERS
    }
    # Separate microproteome track (P16.5): counted APART from the canonical CDS / dark stock.
    mp_total = db.scalar(select(func.count()).select_from(Microprotein)) or 0
    mp_by_class = dict(
        db.execute(select(Microprotein.overlap_class, func.count())
                   .group_by(Microprotein.overlap_class)).all()
    ) if mp_total else {}
    return {"total_genes": total, "hypothetical": hypo,
            "by_verdict": by_verdict, "layer_coverage": coverage,
            "microproteins": {"total": mp_total, "by_overlap_class": mp_by_class,
                              "note": "MS-proven non-canonical microproteome (de Souza 2024); "
                                      "separate track, NOT counted in total_genes or dark"}}


@api.get("/layers", summary="List the evidence layers and what each one is")
def layers() -> dict:
    return {"count": len(LAYERS), "layers": LAYERS}


@api.get("/genes", summary="Search / filter genes; returns paginated summaries")
def list_genes(
    q: str = Query("", description="free-text on rv / gene name / MTBC0 / product / revised function; a pure number or 1.5Mb/1500kb is treated as an MTBC0 genomic position"),
    verdict: str = Query("", description="requalified | family_assigned | dark"),
    hypothetical: bool | None = Query(None, description="restrict to H37Rv-hypothetical genes"),
    has_layer: str = Query("", description="only genes that carry this evidence layer (e.g. vulnerability)"),
    limit: int = Query(50, ge=1, le=500),
    offset: int = Query(0, ge=0),
    db: Session = Depends(get_db),
) -> dict:
    if verdict and verdict not in VERDICTS:
        raise HTTPException(400, f"verdict must be one of {VERDICTS}")
    if has_layer and has_layer not in LAYERS:
        raise HTTPException(400, f"unknown layer '{has_layer}'; see GET /api/v1/layers")
    stmt = select(Gene)
    if q:
        like = f"%{q.strip()}%"
        conds = [Gene.rv.ilike(like), Gene.gene_name.ilike(like),
                 Gene.mtbc0.ilike(like), Gene.product_h37rv.ilike(like),
                 Gene.product_pgap.ilike(like), Gene.function_revised.ilike(like)]
        s = q.strip().lower().replace(",", ".").replace(" ", "")
        pos = None
        try:
            if s.endswith("mb"):
                pos = int(float(s[:-2]) * 1_000_000)
            elif s.endswith("kb"):
                pos = int(float(s[:-2]) * 1_000)
            elif s.rstrip("bp").isdigit():
                pos = int(s.rstrip("bp"))
        except ValueError:
            pos = None
        if pos is not None:
            conds.append(and_(Gene.start_mtbc0 <= pos, Gene.end_mtbc0 >= pos))
        stmt = stmt.where(or_(*conds))
    if verdict:
        stmt = stmt.where(Gene.verdict == verdict)
    if hypothetical is not None:
        stmt = stmt.where(Gene.is_hypothetical.is_(hypothetical))
    if has_layer:
        stmt = stmt.where(getattr(Gene, has_layer).is_not(None))
    total = db.scalar(select(func.count()).select_from(stmt.subquery())) or 0
    rows = db.scalars(stmt.order_by(Gene.rv).limit(limit).offset(offset)).all()
    return {"total": total, "limit": limit, "offset": offset,
            "results": [_summary(g) for g in rows]}


@api.get("/genes/{rv}", summary="Full record for one gene (all scalar fields + all layers)")
def get_gene(rv: str, db: Session = Depends(get_db)) -> dict:
    g = db.get(Gene, rv)
    if not g:
        raise HTTPException(404, f"gene '{rv}' not found")
    return _full(g)


@api.get("/genes/{rv}/{layer}", summary="One evidence layer (JSON) for one gene")
def get_layer(rv: str, layer: str, db: Session = Depends(get_db)) -> dict:
    if layer not in LAYERS:
        raise HTTPException(404, f"unknown layer '{layer}'; see GET /api/v1/layers")
    g = db.get(Gene, rv)
    if not g:
        raise HTTPException(404, f"gene '{rv}' not found")
    return {"rv": rv, "layer": layer, "description": LAYERS[layer], "data": getattr(g, layer)}


# --- Microproteome track (P16.5): MS-proven non-canonical smORFs, separate from the CDS set ---

MP_FIELDS = ("id", "overlap_class", "overlap_locus", "tier", "length_aa", "strand",
             "start", "end", "start_codon", "shine_dalgarno", "folding_free_energy",
             "sequence", "localization", "localization_reliability", "essentiality",
             "essentiality_code", "conserved_in", "host_gene", "host_rv",
             "surrounding_genes", "ms_evidence", "source")
OVERLAP_CLASSES = ("intergenic", "antisense", "sense_overlap")


def _mp(m: Microprotein, full: bool = True) -> dict:
    if not full:
        return {"id": m.id, "overlap_class": m.overlap_class, "overlap_locus": m.overlap_locus,
                "length_aa": m.length_aa, "localization": m.localization,
                "essentiality": m.essentiality, "url": f"/microprotein/{m.id}"}
    return {k: getattr(m, k) for k in MP_FIELDS}


@api.get("/microproteins", summary="List MS-proven non-canonical microproteins (separate track)")
def list_microproteins(
    overlap_class: str = Query("", description="intergenic | antisense | sense_overlap"),
    essential: bool | None = Query(None, description="restrict to TnSeq-essential (ES) microproteins"),
    host_rv: str = Query("", description="microproteins overlapping this H37Rv locus"),
    limit: int = Query(100, ge=1, le=500),
    offset: int = Query(0, ge=0),
    db: Session = Depends(get_db),
) -> dict:
    if overlap_class and overlap_class not in OVERLAP_CLASSES:
        raise HTTPException(400, f"overlap_class must be one of {OVERLAP_CLASSES}")
    stmt = select(Microprotein)
    if overlap_class:
        stmt = stmt.where(Microprotein.overlap_class == overlap_class)
    if essential is True:
        stmt = stmt.where(Microprotein.essentiality_code == "ES")
    if host_rv:
        stmt = stmt.where(Microprotein.host_rv == host_rv)
    total = db.scalar(select(func.count()).select_from(stmt.subquery())) or 0
    rows = db.scalars(stmt.order_by(Microprotein.id).limit(limit).offset(offset)).all()
    return {"total": total, "limit": limit, "offset": offset,
            "note": "MS-proven non-canonical microproteome (de Souza 2024); separate from the CDS set",
            "results": [_mp(m, full=False) for m in rows]}


@api.get("/microproteins/{mp_id}", summary="Full record for one microprotein")
def get_microprotein(mp_id: str, db: Session = Depends(get_db)) -> dict:
    m = db.get(Microprotein, mp_id)
    if not m:
        raise HTTPException(404, f"microprotein '{mp_id}' not found")
    return _mp(m)
