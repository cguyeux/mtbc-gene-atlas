import re
from pathlib import Path

from fastapi import Depends, FastAPI, Form, HTTPException, Request
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from sqlalchemy import and_, func, or_, select
from sqlalchemy.orm import Session

from . import resistance, tbannotator
from .config import settings
from .db import get_db, init_db
from .models import Gene, Microprotein

SITE_ROOT = Path(__file__).resolve().parent.parent
FRONTEND = SITE_ROOT / "frontend"

app = FastAPI(title=settings.site_title, version="0.1.0")
templates = Jinja2Templates(directory=str(FRONTEND / "templates"))
app.mount("/static", StaticFiles(directory=str(FRONTEND / "static")), name="static")

# Public, citable REST API (self-contained sub-app: docs at /api/v1/docs).
from .api import api as api_v1  # noqa: E402
app.mount("/api/v1", api_v1)

VERDICT_LABEL = {
    "requalified": "Resolved",
    "family_assigned": "Family assigned",
    "dark": "Still unknown",
}

# COG functional category single-letter codes (eggNOG / NCBI COG).
COG_LABEL = {
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall / membrane / envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion and vesicular transport",
    "O": "Post-translational modification, protein turnover, chaperones",
    "X": "Mobilome: prophages, transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
}


def render(name: str, request: Request, **ctx):
    return templates.TemplateResponse(
        request, name,
        {"site_title": settings.site_title, "verdict_label": VERDICT_LABEL,
         "cog_label": COG_LABEL, **ctx},
    )


def _stats(db: Session) -> dict:
    total = db.scalar(select(func.count()).select_from(Gene)) or 0
    hypo = db.scalar(select(func.count()).where(Gene.is_hypothetical.is_(True))) or 0
    enriched = db.scalar(select(func.count()).where(Gene.enriched.is_(True))) or 0
    by_verdict = dict(
        db.execute(
            select(Gene.verdict, func.count()).where(Gene.verdict.is_not(None)).group_by(Gene.verdict)
        ).all()
    )
    # How many H37Rv-hypothetical already get a non-hypothetical PGAP product.
    pgap_rescued = db.scalar(
        select(func.count()).where(
            Gene.is_hypothetical.is_(True),
            Gene.product_pgap.is_not(None),
            func.lower(Gene.product_pgap) != "hypothetical protein",
        )
    ) or 0
    auto = db.scalar(select(func.count()).where(Gene.auto.is_(True))) or 0
    review = db.scalar(select(func.count()).where(Gene.needs_review.is_(True))) or 0
    return {"total": total, "hypothetical": hypo, "enriched": enriched,
            "by_verdict": by_verdict, "pgap_rescued": pgap_rescued,
            "auto": auto, "hand": enriched - auto, "needs_review": review}


@app.on_event("startup")
def _startup() -> None:
    init_db()
    if settings.auto_ingest:
        from .db import SessionLocal
        with SessionLocal() as db:
            if (db.scalar(select(func.count()).select_from(Gene)) or 0) == 0:
                from .ingest import ingest
                ingest()


@app.get("/")
def home(request: Request, db: Session = Depends(get_db)):
    stats = _stats(db)
    stats["microproteins"] = db.scalar(select(func.count()).select_from(Microprotein)) or 0
    showcase = db.scalars(
        select(Gene).where(Gene.enriched.is_(True)).order_by(Gene.rv).limit(12)
    ).all()
    return render("home.html", request, stats=stats, showcase=showcase)


def _parse_coord(q: str) -> int | None:
    """Interprète une requête comme une position génomique MTBC0 (bp). Accepte 2367322, 1.5Mb, 1500kb."""
    s = q.strip().lower().replace(",", ".").replace(" ", "")
    try:
        if s.endswith("mb"):
            return int(float(s[:-2]) * 1_000_000)
        if s.endswith("kb"):
            return int(float(s[:-2]) * 1_000)
        if s.endswith("bp"):
            s = s[:-2]
        if s.isdigit():
            return int(s)
    except ValueError:
        return None
    return None


@app.get("/genes")
def genes(request: Request, q: str = "", verdict: str = "", hypo: str = "",
          enriched: str = "", review: str = "", page: int = 1, db: Session = Depends(get_db)):
    per_page = 50
    stmt = select(Gene)
    if q:
        like = f"%{q.strip()}%"
        conds = [Gene.rv.ilike(like), Gene.gene_name.ilike(like),
                 Gene.mtbc0.ilike(like), Gene.product_h37rv.ilike(like),
                 Gene.product_pgap.ilike(like), Gene.function_revised.ilike(like)]
        pos = _parse_coord(q)
        if pos is not None:   # recherche par POSITION génomique (MTBC0) : gène chevauchant
            conds.append(and_(Gene.start_mtbc0 <= pos, Gene.end_mtbc0 >= pos))
        stmt = stmt.where(or_(*conds))
    if verdict:
        stmt = stmt.where(Gene.verdict == verdict)
    if hypo == "1":
        stmt = stmt.where(Gene.is_hypothetical.is_(True))
    if enriched == "1":
        stmt = stmt.where(Gene.enriched.is_(True))
    if review == "1":
        stmt = stmt.where(Gene.needs_review.is_(True))
    total = db.scalar(select(func.count()).select_from(stmt.subquery())) or 0
    page = max(1, page)
    rows = db.scalars(stmt.order_by(Gene.rv).limit(per_page).offset((page - 1) * per_page)).all()
    return render("genes_list.html", request, rows=rows, total=total, page=page,
                  per_page=per_page, q=q, verdict=verdict, hypo=hypo, enriched=enriched,
                  review=review, n_pages=(total + per_page - 1) // per_page)


@app.get("/gene/{rv}")
def gene(rv: str, request: Request, db: Session = Depends(get_db)):
    g = db.get(Gene, rv)
    if not g:
        return render("not_found.html", request, rv=rv)
    micro = db.scalars(
        select(Microprotein).where(Microprotein.host_rv == rv).order_by(Microprotein.id)
    ).all()
    return render("gene.html", request, g=g, micro=micro,
                  track=_browse_track(db, center_rv=rv, span=24000, width=760),
                  net=_network_graph(db, rv, top=10, link_mode="gene"))


# --- Microproteome track (P16.5): MS-proven non-canonical smORFs, separate from the 3906 CDS ---

OVERLAP_LABEL = {
    "intergenic": "Intergenic (novel locus)",
    "antisense": "Antisense to a gene",
    "sense_overlap": "Same-strand overlap (alternative frame)",
}


@app.get("/microproteins")
def microproteins(request: Request, overlap: str = "", essential: str = "", method: str = "",
                  q: str = "", page: int = 1, db: Session = Depends(get_db)):
    per_page = 60
    stmt = select(Microprotein)
    if overlap in OVERLAP_LABEL:
        stmt = stmt.where(Microprotein.overlap_class == overlap)
    if essential == "1":
        stmt = stmt.where(Microprotein.essentiality_code == "ES")
    if method:
        stmt = stmt.where(Microprotein.detection_method == method)
    if q:
        like = f"%{q.strip()}%"
        stmt = stmt.where(or_(Microprotein.id.ilike(like), Microprotein.host_gene.ilike(like),
                              Microprotein.host_rv.ilike(like), Microprotein.localization.ilike(like)))
    total = db.scalar(select(func.count()).select_from(stmt.subquery())) or 0
    by_class = dict(db.execute(select(Microprotein.overlap_class, func.count())
                               .group_by(Microprotein.overlap_class)).all())
    by_method = dict(db.execute(select(Microprotein.detection_method, func.count())
                                .group_by(Microprotein.detection_method)).all())
    n_ess = db.scalar(select(func.count()).where(Microprotein.essentiality_code == "ES")) or 0
    page = max(1, page)
    rows = db.scalars(stmt.order_by(Microprotein.overlap_class, Microprotein.id)
                      .limit(per_page).offset((page - 1) * per_page)).all()
    return render("microproteins_list.html", request, rows=rows, total=total, page=page,
                  per_page=per_page, overlap=overlap, essential=essential, method=method, q=q,
                  by_class=by_class, by_method=by_method, n_ess=n_ess, overlap_label=OVERLAP_LABEL,
                  n_pages=(total + per_page - 1) // per_page)


@app.get("/microprotein/{mp_id}")
def microprotein(mp_id: str, request: Request, db: Session = Depends(get_db)):
    m = db.get(Microprotein, mp_id)
    if not m:
        return render("not_found.html", request, rv=mp_id)
    host = db.get(Gene, m.host_rv) if m.host_rv else None
    return render("microprotein.html", request, m=m, host=host, overlap_label=OVERLAP_LABEL)


# --- Genome browser (linear track, MTBC0 coordinates, H37Rv locus tags) ---

VERDICT_COLOR = {"requalified": "#2e7d32", "family_assigned": "#1565c0", "dark": "#8a8a8a"}
BROWSE_W = 1000          # SVG viewBox width
BROWSE_SPAN = 40000      # default window (bp)


def _browse_track(db: Session, center_rv: str = "", center_pos: int = 0,
                  span: int = BROWSE_SPAN, width: int = BROWSE_W) -> dict:
    """Build a linear genome-track view around a gene or position (server-rendered SVG geometry)."""
    span = max(3000, min(int(span), 400000))
    center = None
    cg = db.get(Gene, center_rv) if center_rv else None
    if cg and cg.start_mtbc0 and cg.end_mtbc0:
        center = (cg.start_mtbc0 + cg.end_mtbc0) // 2
    if center is None:
        center = int(center_pos) or 2_200_000
    half = span // 2
    w0 = max(0, center - half)
    w1 = w0 + span
    genes = db.scalars(
        select(Gene).where(Gene.start_mtbc0.is_not(None), Gene.end_mtbc0.is_not(None),
                           Gene.start_mtbc0 <= w1, Gene.end_mtbc0 >= w0)
        .order_by(Gene.start_mtbc0)
    ).all()
    scale = width / span
    AH, TIP, Y_PLUS, Y_MINUS = 15, 7, 26, 63   # arrow height, tip length, top-y of +/- rows
    items = []
    for g in genes:
        x = max(0.0, (g.start_mtbc0 - w0) * scale)
        w = max(5.0, (min(g.end_mtbc0, w1) - max(g.start_mtbc0, w0)) * scale)
        plus = (g.strand or "+") != "-"
        y = Y_PLUS if plus else Y_MINUS
        tip = min(TIP, w)
        if plus:   # pointe à droite
            pts = f"{x:.1f},{y} {x+w-tip:.1f},{y} {x+w:.1f},{y+AH/2:.1f} {x+w-tip:.1f},{y+AH} {x:.1f},{y+AH}"
        else:      # pointe à gauche
            pts = f"{x+w:.1f},{y} {x+tip:.1f},{y} {x:.1f},{y+AH/2:.1f} {x+tip:.1f},{y+AH} {x+w:.1f},{y+AH}"
        items.append({
            "rv": g.rv, "name": g.gene_name or g.rv,
            "verdict": g.verdict, "color": VERDICT_COLOR.get(g.verdict, "#bdbdbd"),
            "x": round(x, 1), "w": round(w, 1), "pts": pts,
            "cx": round(x + w / 2, 1), "ty": y - 4 if plus else y + AH + 11,
            "is_center": bool(center_rv) and g.rv == center_rv,
            "product": (g.product_pgap or g.product_h37rv or "")[:60],
            "label": w >= 26,
        })
    # ruler ticks (~6)
    import math
    step = max(500, int(round(span / 6 / 500)) * 500)
    ticks = []
    t = (w0 // step + 1) * step
    while t < w1:
        ticks.append({"pos": t, "x": round((t - w0) * scale, 1), "kb": round(t / 1000)})
        t += step
    return {"items": items, "w0": w0, "w1": w1, "span": span, "center": center,
            "center_rv": center_rv, "width": width,
            "left": max(0, center - span), "right": center + span,
            "zin": max(3000, span // 2), "zout": min(400000, span * 2), "ticks": ticks}


@app.get("/browse")
def browse(request: Request, center: str = "", pos: int = 0, span: int = BROWSE_SPAN,
           db: Session = Depends(get_db)):
    track = _browse_track(db, center_rv=center.strip(), center_pos=pos, span=span, width=BROWSE_W)
    return render("browse.html", request, track=track)


# --- Functional interaction network (STRING ego-network, server-rendered SVG) ---

import math as _math

NET_W, NET_H = 900, 560


def _edge_kind(p: dict) -> tuple[str, str]:
    """(kind, colour) d'une arête d'après ses canaux STRING."""
    ch = p.get("channels") or {}
    if (ch.get("experimental") or 0) >= 400 or (ch.get("database") or 0) >= 400:
        return "experimental", "#2e7d32"          # mesuré
    if p.get("context_driven") or (ch.get("neighborhood") or 0) >= 400 or (ch.get("cooccurence") or 0) >= 400 \
            or (ch.get("fusion") or 0) >= 400:
        return "context", "#ef6c00"               # contexte génomique
    return "other", "#9e9e9e"                     # coexpression / text-mining


def _top_partners(g: Gene, k: int) -> list[dict]:
    ps = (g.string or {}).get("partners") or []
    return sorted(ps, key=lambda p: -(p.get("combined_no_tm") or 0))[:k]


def _network_graph(db: Session, center_rv: str, top: int = 12, hops: int = 1,
                   link_mode: str = "recenter") -> dict | None:
    cg = db.get(Gene, center_rv)
    if not cg or not (cg.string and cg.string.get("partners")):
        return None
    top = max(3, min(top, 24))
    hops = 2 if hops >= 2 else 1
    MAX_NODES = 46
    # BFS des nœuds par saut
    hop = {cg.rv: 0}
    order = [cg.rv]
    h1 = _top_partners(cg, top)
    for p in h1:
        if p["rv"] not in hop:
            hop[p["rv"]] = 1
            order.append(p["rv"])
    if hops == 2:
        h1genes = {g.rv: g for g in db.scalars(select(Gene).where(Gene.rv.in_([p["rv"] for p in h1]))).all()}
        for p in h1:
            g1 = h1genes.get(p["rv"])
            if not g1:
                continue
            for q in _top_partners(g1, 6):
                if q["rv"] not in hop and len(hop) < MAX_NODES:
                    hop[q["rv"]] = 2
                    order.append(q["rv"])
    genes = {g.rv: g for g in db.scalars(select(Gene).where(Gene.rv.in_(order))).all()}
    order = [rv for rv in order if rv in genes]
    nodeset = set(order)

    cx, cy = NET_W / 2, NET_H / 2
    R1, R2 = min(NET_W, NET_H) / 2 - 130, min(NET_W, NET_H) / 2 - 40
    href = (lambda rv: f"/gene/{rv}") if link_mode == "gene" else (lambda rv: f"/network?center={rv}")
    h1list = [rv for rv in order if hop[rv] == 1]
    h2list = [rv for rv in order if hop[rv] == 2]
    pos = {cg.rv: (cx, cy)}
    for i, rv in enumerate(h1list):
        a = 2 * _math.pi * i / max(1, len(h1list)) - _math.pi / 2
        pos[rv] = (cx + R1 * _math.cos(a), cy + R1 * _math.sin(a))
    for i, rv in enumerate(h2list):
        a = 2 * _math.pi * i / max(1, len(h2list)) - _math.pi / 2
        pos[rv] = (cx + R2 * _math.cos(a), cy + R2 * _math.sin(a))

    nodes = []
    for rv in order:
        g = genes[rv]
        x, y = pos[rv]
        nodes.append({
            "rv": rv, "name": g.gene_name or rv, "x": round(x, 1), "y": round(y, 1),
            "r": 22 if rv == cg.rv else (13 if not g.is_hypothetical else 12),
            "color": VERDICT_COLOR.get(g.verdict, "#bdbdbd"), "verdict": g.verdict,
            "is_center": rv == cg.rv, "hop": hop[rv], "hypothetical": bool(g.is_hypothetical),
            "href": f"/gene/{rv}" if rv == cg.rv else href(rv),
            "tx": round(x + (16 if x >= cx else -16), 1), "anchor": "start" if x >= cx else "end",
            "title": f"{g.gene_name or rv} ({g.verdict or 'n/a'}) — {(g.product_pgap or g.product_h37rv or '')[:50]}",
        })
    # arêtes entre nœuds inclus, dérivées des listes de partenaires (score/canal conservés)
    edges = []
    seen = set()
    for a in order:
        for q in ((genes[a].string or {}).get("partners") or []):
            b = q["rv"]
            if b not in nodeset or a == b:
                continue
            key = tuple(sorted((a, b)))
            if key in seen:
                continue
            seen.add(key)
            kind, col = _edge_kind(q)
            sc = q.get("combined_no_tm") or 0
            x1, y1 = pos[a]; x2, y2 = pos[b]
            edges.append({"source": key[0], "target": key[1], "kind": kind, "score": sc, "color": col,
                          "x1": round(x1, 1), "y1": round(y1, 1), "x2": round(x2, 1), "y2": round(y2, 1),
                          "w": round(0.8 + 3.0 * sc / 1000, 2), "op": round(0.22 + 0.5 * sc / 1000, 2)})
    return {"nodes": nodes, "edges": edges, "center_rv": cg.rv, "n_partners": len(h1list),
            "n_nodes": len(nodes), "width": NET_W, "height": NET_H, "top": top, "hops": hops,
            "anchor": (cg.string.get("anchor") or {}).get("gene")}


@app.get("/network")
def network(request: Request, center: str = "", top: int = 12, hops: int = 1,
            db: Session = Depends(get_db)):
    graph = _network_graph(db, center.strip(), top=top, hops=hops, link_mode="recenter") if center.strip() else None
    return render("network.html", request, graph=graph, center=center.strip(), hops=2 if hops >= 2 else 1)


@app.get("/api/genes")
def api_genes(q: str = "", limit: int = 50, db: Session = Depends(get_db)):
    like = f"%{q.strip()}%"
    rows = db.scalars(
        select(Gene).where(or_(Gene.rv.ilike(like), Gene.gene_name.ilike(like),
                               Gene.product_pgap.ilike(like))).order_by(Gene.rv).limit(limit)
    ).all()
    return JSONResponse([
        {"rv": g.rv, "gene": g.gene_name, "product": g.product_pgap or g.product_h37rv,
         "verdict": g.verdict, "enriched": g.enriched} for g in rows
    ])


# Feedback is now a client-side mailto: link on each gene page (honest delivery on a scale-to-zero
# container: no ephemeral DB black hole, no unread inbox). The old POST /feedback endpoint was removed.


@app.get("/changelog.csv")
def changelog_csv():
    """Per-gene changelog: Mycobrowser baseline -> atlas result + what is new (one row per gene).
    Baked into the image by analyses/phase65_atlas_changelog.py."""
    path = settings.genes_content.parent / "atlas_changelog.csv"
    if not path.exists():
        raise HTTPException(status_code=404, detail="changelog not built")
    return FileResponse(path, media_type="text/csv", filename="mtbc_atlas_changelog.csv")


@app.get("/about")
def about(request: Request):
    return render("about.html", request)


# --- Resistance tester (second tool, deterministic WHO-catalogue baseline) ---

def _looks_like_vcf(text: str) -> bool:
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        return len(f) >= 5 and f[1].isdigit()
    return False


_SRA_RE = re.compile(r"^(SRR|ERR|DRR|SRS|ERS|DRS|SAMN|SAMEA|SAMD)\w+$", re.I)


def _looks_like_sra(text: str) -> bool:
    return bool(_SRA_RE.match(text.strip()))


def _looks_like_mutation(text: str) -> bool:
    t = text.strip()
    return ("\n" not in t) and ("NC_000962.3:" not in t) and len(t.split()) <= 2 and any(c.isdigit() for c in t)


def _mdr_summary(prof: dict) -> dict:
    r = lambda d: prof.get(d, {}).get("verdict") == "R"
    mdr = r("isoniazid") and r("rifampicin")
    fq = r("moxifloxacin") or r("levofloxacin")
    return {"mdr": mdr, "pre_xdr": mdr and fq}


@app.get("/resistance")
def resistance_form(request: Request):
    return render("resistance.html", request, stats=resistance.catalogue_stats())


@app.post("/resistance")
def resistance_run(request: Request, payload: str = Form(""), mode: str = Form("auto")):
    text = (payload or "").strip()
    spdis: set[str] = set()
    mutation_info = sra_info = error = None
    parsed_as = mode
    if text:
        if mode == "sra" or (mode == "auto" and _looks_like_sra(text)):
            parsed_as = "sra"
            try:
                spdis = set(tbannotator.fetch_spdi(text))
                sra_info = {"accession": text.strip(), "n": len(spdis)}
                if not spdis:
                    error = (f"No variants found for {text.strip()} in TBannotator "
                             "(unknown or not yet annotated accession).")
            except tbannotator.TBannotatorError as e:
                sra_info = {"accession": text.strip(), "n": 0}
                error = (f"Could not resolve {text.strip()}: {e}. "
                         "You can still paste the strain's variants directly.")
        elif mode == "mutation" or (mode == "auto" and _looks_like_mutation(text)):
            parsed_as = "mutation"
            spdi, ok = resistance.resolve_mutation(text)
            mutation_info = {"input": text, "spdi": spdi, "recognised": ok}
            if spdi:
                spdis = {spdi}
        elif mode == "vcf" or (mode == "auto" and _looks_like_vcf(text)):
            parsed_as = "vcf"
            spdis = resistance.parse_vcf_text(text)
        else:
            parsed_as = "spdi"
            spdis = resistance.parse_spdi_text(text)
    prof = resistance.profile(spdis) if spdis else {}
    resistant = {d: v for d, v in prof.items() if v["verdict"] == "R"}
    susceptible = [d for d, v in prof.items() if v["verdict"] != "R"]
    return render("resistance.html", request, stats=resistance.catalogue_stats(),
                  payload=text, mode=mode, parsed_as=parsed_as, n_input=len(spdis),
                  prof=prof, resistant=resistant, susceptible=susceptible,
                  mutation_info=mutation_info, sra_info=sra_info, error=error,
                  mdr=_mdr_summary(prof), ran=True)
