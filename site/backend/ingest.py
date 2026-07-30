"""Populate the database from the pipeline artefacts.

Two layers:
  1. Catalogue   -- every joined gene from gene_xref.tsv (3022 rows).
  2. Enriched    -- per-gene JSON records from content/genes/<rv>.json,
                    which add curated function, verdict, Pfam domains,
                    Foldseek hits, ESM signal, evidence and references.

Run standalone:  python -m backend.ingest
"""
from __future__ import annotations

import csv
import json

from .config import settings
from .db import SessionLocal, init_db
from .models import Gene, Microprotein

_MP_FIELDS = (
    "kind", "overlap_class", "overlap_locus", "tier", "ms_evidence", "detection_method",
    "extra", "start", "end", "strand", "length_aa", "start_codon", "shine_dalgarno",
    "folding_free_energy", "sequence", "localization", "localization_reliability",
    "essentiality", "essentiality_code", "conserved_in", "host_gene", "host_rv",
    "surrounding_genes", "source",
)


def ingest_microproteins(db) -> int:
    """Load the separate microproteome track (P16.5), all sources. Idempotent; no-op if absent.

    Loads every content/microproteins*.json (de Souza MS compendium, Troian VapC4
    ribosome-stalling, ...); entries keyed by id, get-or-create."""
    import glob
    n = 0
    for path in sorted(glob.glob(str(settings.microproteins_content.parent / "microproteins*.json"))):
        for rec in json.load(open(path)):
            mp = db.get(Microprotein, rec["id"]) or Microprotein(id=rec["id"])
            for k in _MP_FIELDS:
                setattr(mp, k, rec.get(k))
            db.add(mp)
            n += 1
    return n


def ingest() -> dict:
    init_db()
    xref = settings.resolved_xref()
    db = SessionLocal()
    try:
        n = 0
        seen: set[str] = set()
        with open(xref) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rv = r["rv"]
                # A few H37Rv loci map to >1 MTBC0 gene (paralogues annotated by
                # homology to the same RefSeq). Keep the first occurrence; the
                # catalogue is H37Rv-centric (one row per Rv).
                if rv in seen:
                    continue
                seen.add(rv)
                g = db.get(Gene, rv) or Gene(rv=rv)
                g.mtbc0 = r.get("mtbc0")
                g.np_id = r.get("np_id")
                g.gene_name = r.get("gene") or None
                g.len_aa = int(r["len_aa"]) if r.get("len_aa") else None
                g.strand = r.get("strand")
                g.start_mtbc0 = int(r["start_mtbc0"]) if r.get("start_mtbc0") else None
                g.end_mtbc0 = int(r["end_mtbc0"]) if r.get("end_mtbc0") else None
                g.product_h37rv = r.get("product_h37rv")
                g.product_pgap = r.get("product_mtbc0_pgap") or None
                g.is_hypothetical = r.get("is_hypothetical_h37rv") == "True"
                g.protein_mtbc0 = r.get("protein_mtbc0")
                db.add(g)
                n += 1
        db.commit()

        m = 0
        if settings.genes_content.exists():
            for p in sorted(settings.genes_content.glob("*.json")):
                rec = json.loads(p.read_text())
                g = db.get(Gene, rec["rv"]) or Gene(rv=rec["rv"])
                # Hand-curated fiches may name a previously-unnamed locus; let the
                # curation layer's gene symbol override the catalogue (else keep it).
                if rec.get("gene"):
                    g.gene_name = rec["gene"]
                # The enriched fiche carries the AUTHORITATIVE MTBC0 coordinates/sequence:
                # gene_xref.tsv can hold a stale/duplicate row for a locus (e.g. Rv2082 had a
                # 54-aa fragment row + the real 721-aa row), which would misplace the gene in the
                # browser and break position search. The per-gene fiche is the source of truth.
                if rec.get("start_mtbc0"):
                    g.start_mtbc0 = rec["start_mtbc0"]
                if rec.get("end_mtbc0"):
                    g.end_mtbc0 = rec["end_mtbc0"]
                if rec.get("strand"):
                    g.strand = rec["strand"]
                if rec.get("protein_mtbc0"):
                    g.protein_mtbc0 = rec["protein_mtbc0"]
                    g.len_aa = len(rec["protein_mtbc0"])
                g.function_revised = rec.get("function_revised") or None
                g.verdict = rec.get("verdict") or None
                g.confidence = rec.get("confidence") or None
                g.curation_note = rec.get("curation_note") or None
                g.esm = rec.get("esm")
                g.domains = rec.get("domains")
                g.struct_hits = rec.get("struct_hits")
                g.eggnog = rec.get("eggnog") or None
                g.uniprot = rec.get("uniprot") or None
                g.conservation = rec.get("conservation") or None
                g.plddt = rec.get("plddt") or None
                g.struct_af = rec.get("struct_af") or None
                g.plddt_af = rec.get("plddt_af") or None
                g.string = rec.get("string") or None
                g.mcsa = rec.get("mcsa") or None
                g.pfam_tentative = rec.get("pfam_tentative") or None
                g.localization = rec.get("localization") or None
                g.essentiality = rec.get("essentiality") or None
                g.proteomics = rec.get("proteomics") or None
                g.protparam = rec.get("protparam") or None
                g.resistance = rec.get("resistance") or None
                g.funccat = rec.get("funccat") or None
                g.orthologs = rec.get("orthologs") or None
                g.mutant_phenotypes = rec.get("mutant_phenotypes") or None
                g.pdb = rec.get("pdb") or None
                g.genomic_context = rec.get("genomic_context") or None
                g.regulation = rec.get("regulation") or None
                g.mycobrowser = rec.get("mycobrowser") or None
                g.rd = rec.get("rd") or None
                g.ec_override = rec.get("ec_override") or None
                g.phenotype_lead = rec.get("phenotype_lead") or None
                g.hhpred = rec.get("hhpred") or None
                g.struct_cluster = rec.get("struct_cluster") or None
                g.context_lead = rec.get("context_lead") or None
                g.vulnerability = rec.get("vulnerability") or None
                g.ptm = rec.get("ptm") or None
                g.expression = rec.get("expression") or None
                g.integrative_lead = rec.get("integrative_lead") or None
                g.outgroup = rec.get("outgroup") or None
                g.rbtnseq_fitness = rec.get("rbtnseq_fitness") or None
                g.literature = rec.get("literature") or None
                g.disorder = rec.get("disorder") or None
                g.prospect = rec.get("prospect") or None
                g.abpp = rec.get("abpp") or None
                g.evidence = rec.get("evidence")
                g.refs = rec.get("references")
                g.auto = bool(rec.get("auto", False))
                g.needs_review = bool(rec.get("needs_review", False))
                # P6.2 -- coherence guard verdict<->function. A hand-curated fiche
                # (auto=false) that carries an assigned function but a verdict left at
                # "dark" is a back-propagation inconsistency (a companion layer wrote
                # function_revised without updating verdict). Recale it at ingest so the
                # atlas never serves "dark" over an already-resolved function again.
                if g.verdict == "dark" and not g.auto:
                    fr = (g.function_revised or "").strip().lower()
                    dark_conclusive = (
                        not fr
                        or fr.startswith("conserved hypothetical")
                        or fr.startswith("conserved protein of unknown")
                        or "function unknown" in fr
                        or "function remains unknown" in fr
                    )
                    if not dark_conclusive:
                        g.verdict = "family_assigned"
                        g.confidence = g.confidence or "low"
                g.enriched = True
                db.add(g)
                m += 1
        mp = ingest_microproteins(db)
        db.commit()
        return {"catalogue": n, "enriched": m, "microproteins": mp}
    finally:
        db.close()


if __name__ == "__main__":
    print(ingest())
