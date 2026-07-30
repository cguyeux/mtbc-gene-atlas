#!/usr/bin/env python3
"""phase13_pipeline_schema.py -- schéma du pipeline d'annotation (Figure 1 du manuscrit).

Flux : protéine ancestrale ancrée MTBC0 -> 7 couches de preuve -> fiche consolidée -> ressource servie.
Sortie : article/figures/fig_pipeline.{pdf,png} (vectoriel + PNG 600 DPI). Matplotlib pur.
"""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = Path(__file__).resolve().parent.parent / "article" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.size": 8.5, "figure.dpi": 150, "savefig.dpi": 600, "savefig.bbox": "tight"})

C_IN = "#4C72B0"; C_OUT = "#C44E52"; C_REC = "#55A868"
BANDS = [
    ("Sequence & structure", "#EAF0F7", "#4C72B0", [
        "Pfam\n(HMMER3)",
        "ESMFold\n+ Foldseek",
        "ESM Atlas\n(PLM)"]),
    ("Orthology & curation", "#FBEEE6", "#DD8452", [
        "eggNOG-mapper\n(COG/EC/KO/GO)",
        "UniProt\n(SwissProt)"]),
    ("Population & network", "#EAF3EC", "#55A868", [
        "STRING v12\n(interactions)",
        "Intra-MTBC pN/pS\n(145,209 genomes)"]),
]


def box(ax, x, y, w, h, text, fc, ec, fs=8.5, bold=False, lw=1.1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                                fc=fc, ec=ec, lw=lw))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", color="#222222")


def arrow(ax, x0, y0, x1, y1, color="#888888"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=11,
                                 lw=1.2, color=color, shrinkA=0, shrinkB=0))


fig, ax = plt.subplots(figsize=(9.2, 4.5))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

# ── Entrée ────────────────────────────────────────────────────────────────────
box(ax, 0.5, 40, 15, 20, "Ancestral protein\nanchored on MTBC0\n\n3,906 H37Rv CDS\n→ MTBC0 orthologue",
    "#DCE6F2", C_IN, fs=8.0, bold=False)

# ── Bandes de couches ──────────────────────────────────────────────────────────
bx, bw = 21, 38
band_y = [66, 36, 6]; band_h = 26
for (title, fc, ec, items), by in zip(BANDS, band_y):
    box(ax, bx, by, bw, band_h, "", fc, ec, lw=1.2)
    ax.text(bx + bw / 2, by + band_h - 3.0, title, ha="center", va="center",
            fontsize=8.6, fontweight="bold", color=ec)
    n = len(items); gap = 2.2; iw = (bw - gap * (n + 1)) / n
    for k, it in enumerate(items):
        ix = bx + gap * (k + 1) + iw * k
        box(ax, ix, by + 2.4, iw, band_h - 9.2, it, "white", ec, fs=7.2)
    arrow(ax, 15.5, 50, bx, by + band_h / 2)

# ── Fiche consolidée ─────────────────────────────────────────────────────────
rx = 63
box(ax, rx, 33, 18, 34,
    "Consolidated\nper-gene record\n\nverdict · confidence\ndated sources\ncross-layer note",
    "#E3F0E7", C_REC, fs=8.0, bold=False)
for by in band_y:
    arrow(ax, bx + bw, by + band_h / 2, rx, 50)

# ── Ressource servie ────────────────────────────────────────────────────────
ox = 84.5
box(ax, ox, 36, 15, 28,
    "Containerised\nliving resource\n\nFastAPI\n/gene/<Rv>\nJSON API",
    "#F7E3E5", C_OUT, fs=8.0, bold=False)
arrow(ax, rx + 18, 50, ox, 50, color=C_REC)

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"fig_pipeline.{ext}")
print("wrote fig_pipeline.pdf / .png")
