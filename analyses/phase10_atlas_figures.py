#!/usr/bin/env python3
"""phase10_atlas_figures.py -- figures publication du manuscrit atlas (depuis les fiches du site).

Produit dans article/figures/ (PDF vectoriel + PNG 600 DPI) :
  fig1_coverage.{pdf,png}     (A) couverture par couche de preuve, protéome entier vs hypothétiques ;
                              (B) devenir des gènes historiquement « hypothetical ».
  fig2_conservation.{pdf,png} (A) distribution pN/pS intra-MTBC des hypothétiques ;
                              (B) confiance pLDDT des 350 modèles ESMFold (gating à 70).

Source de vérité : site/content/genes/*.json (mêmes chiffres que claim_check.md). Stdlib + matplotlib.
"""
from __future__ import annotations
import glob, json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "article" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.dpi": 600, "savefig.bbox": "tight",
})
C_ALL, C_HYP = "#4C72B0", "#C44E52"  # protéome entier / hypothétiques


def ne(x):
    return bool(x) and str(x).strip() not in ("", "[]", "{}")


def load():
    g = [json.load(open(f)) for f in glob.glob(str(GENES / "*.json"))]
    hyp = [d for d in g if "hypothetical" in (d.get("product_h37rv") or "").lower()]
    return g, hyp


def layer_flags(d):
    egg = d.get("eggnog") or {}
    up = d.get("uniprot") or {}
    return {
        "Pfam domain": ne(d.get("domains")),
        "Structural hit": any(h.get("significant") for h in (d.get("struct_hits") or [])),
        "eggNOG (COG/EC/KO)": any(ne(egg.get(k)) for k in ("cog_cat", "ec", "kegg_ko")),
        "UniProt curated": ne(up.get("function")),
        "STRING anchor": ne((d.get("string") or {}).get("anchor")),
    }


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}")
    plt.close(fig)
    print(f"  wrote {name}.pdf / .png")


def fig_coverage(g, hyp):
    layers = list(layer_flags(g[0]).keys())
    pct_all = [100 * sum(layer_flags(d)[L] for d in g) / len(g) for L in layers]
    pct_hyp = [100 * sum(layer_flags(d)[L] for d in hyp) / len(hyp) for L in layers]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.0, 3.6))
    y = range(len(layers))
    h = 0.38
    axA.barh([i + h/2 for i in y], pct_all, height=h, color=C_ALL, label=f"All genes (n={len(g)})")
    axA.barh([i - h/2 for i in y], pct_hyp, height=h, color=C_HYP, label=f"Hypothetical (n={len(hyp)})")
    axA.set_yticks(list(y)); axA.set_yticklabels(layers); axA.invert_yaxis()
    axA.set_xlabel("Genes covered (%)"); axA.set_xlim(0, 100)
    axA.legend(loc="upper right", bbox_to_anchor=(0.99, 0.82), frameon=False, fontsize=8)
    axA.set_title("(A) Functional-annotation coverage by layer")
    for i, v in enumerate(pct_all):
        axA.text(v + 1, i + h/2, f"{v:.0f}", va="center", fontsize=7, color=C_ALL)
    for i, v in enumerate(pct_hyp):
        axA.text(v + 1, i - h/2, f"{v:.0f}", va="center", fontsize=7, color=C_HYP)

    # (B) fate of the hypotheticals
    H = len(hyp)
    handle = sum(1 for d in hyp if any(ne((d.get("eggnog") or {}).get(k)) for k in ("cog_cat", "ec", "kegg_ko"))
                 or ne((d.get("uniprot") or {}).get("function")))
    purif = sum(1 for d in hyp if "purifying" in ((d.get("conservation") or {}).get("selection") or "").lower())
    anchored = sum(1 for d in hyp if ne((d.get("string") or {}).get("anchor")))
    pseudo = sum(1 for d in hyp if (d.get("conservation") or {}).get("pseudogene_flag") is True)
    cats = ["Functional handle\n(COG/EC/KO or curated)", "STRING-anchored", "Under purifying\nselection", "Pseudogene candidate"]
    vals = [handle, anchored, purif, pseudo]
    axB.barh(range(len(cats)), vals, color=C_HYP, alpha=0.85)
    axB.set_yticks(range(len(cats))); axB.set_yticklabels(cats, fontsize=8); axB.invert_yaxis()
    axB.set_xlabel(f"Hypothetical genes (of {H})"); axB.set_xlim(0, H)
    axB.set_title("(B) Fate of the historically hypothetical genes")
    for i, v in enumerate(vals):
        axB.text(v + H*0.01, i, f"{v} ({round(100*v/H)}%)", va="center", fontsize=8)
    fig.tight_layout()
    save(fig, "fig1_coverage")


def fig_conservation(g, hyp):
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.0, 3.4))
    # (A) pN/pS of hypotheticals
    pnps = [(d.get("conservation") or {}).get("pN_pS") for d in hyp]
    pnps = [x for x in pnps if isinstance(x, (int, float)) and x >= 0]
    clipped = [min(x, 3.0) for x in pnps]
    axA.hist(clipped, bins=40, color=C_HYP, alpha=0.85, edgecolor="white", linewidth=0.3)
    axA.axvline(1.0, color="black", ls="--", lw=1)
    purif = sum(1 for d in hyp if "purifying" in ((d.get("conservation") or {}).get("selection") or "").lower())
    axA.set_xlabel("pN/pS (intra-MTBC, clipped at 3)"); axA.set_ylabel("Hypothetical genes")
    axA.set_title("(A) Cross-strain selection on hypothetical genes")
    axA.annotate(f"{purif}/{len(hyp)} under\npurifying selection\n(pN/pS < 1)",
                 xy=(0.3, 0.92), xycoords="axes fraction", va="top", fontsize=8,
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=C_HYP, lw=0.8))
    axA.text(1.05, axA.get_ylim()[1]*0.96, "neutral", fontsize=7, color="black", va="top")

    # (B) pLDDT of the 350 ESMFold models
    pl = [(d.get("plddt") or {}).get("mean_plddt") for d in g]
    pl = [x for x in pl if isinstance(x, (int, float))]
    below = sum(1 for x in pl if x < 70)
    axB.hist(pl, bins=30, color=C_ALL, alpha=0.85, edgecolor="white", linewidth=0.3)
    axB.axvline(70, color="black", ls="--", lw=1)
    axB.set_xlabel("Mean pLDDT (ESMFold model)"); axB.set_ylabel("Dark genes modelled")
    axB.set_title("(B) Structure-model confidence and gating")
    axB.annotate(f"{below}/{len(pl)} below 70\n(Foldseek hits discounted)",
                 xy=(0.04, 0.92), xycoords="axes fraction", va="top", fontsize=8,
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=C_ALL, lw=0.8))
    fig.tight_layout()
    save(fig, "fig2_conservation")


def main():
    g, hyp = load()
    print(f"Loaded {len(g)} fiches ({len(hyp)} hypothetical). Writing figures to {OUT} ...")
    fig_coverage(g, hyp)
    fig_conservation(g, hyp)
    print("Done.")


if __name__ == "__main__":
    main()
