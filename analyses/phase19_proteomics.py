#!/usr/bin/env python3
"""phase19_proteomics.py -- couche protéomique / spectrométrie de masse (P5.1).

Comble le plus grand manque du gap analysis Mycobrowser : une PREUVE D'EXISTENCE
expérimentale du produit protéique, orthogonale à la conservation et à l'essentialité,
et particulièrement précieuse pour les "hypothetical". Mycobrowser indiquait, même pour
un hypothétique, s'il avait été détecté par MS ; on fait mieux en agrégeant plusieurs
jeux MS indépendants et en donnant l'abondance quantitative.

Source : PaxDb 5.0 (Huang et al. 2023, MCP), organisme M. tuberculosis H37Rv (83332).
  - 1 fichier intégré  : 83332-WHOLE_ORGANISM-integrated.txt (abondance ppm, moyenne pondérée)
  - 16 jeux individuels : chacun un dataset ProteomeXchange/GPM/PeptideAtlas ré-analysé
    uniformément par PaxDb (dont Schubert 2013 Cell Host Microbe complete proteome,
    Albrethsen 2013 MCP log/starvation, Cortes 2013, whole-cell lysate PXD011081, ...).
  Les IDs PaxDb sont "83332.Rvxxxx" -> locus tag H37Rv, mapping DIRECT sur nos fiches.

Signaux dérivés par gène :
  - detected           : présent dans l'intégré OU dans >=1 des 16 jeux
  - abundance_ppm      : abondance intégrée (parties par million de masse protéique)
  - abundance_rank     : rang (1 = plus abondant) parmi les protéines quantifiées
  - abundance_percentile : 100 = plus abondant
  - n_datasets / 16    : dans combien de jeux MS indépendants la protéine est détectée
                         (détection reproductible = preuve d'existence robuste)

Sortie : résultats/phase19_proteomics/proteomics.json (keyed Rv) + fusion directe dans
les fiches (patron phase8/9/14, idempotent, évite un re-run phase4 qui écraserait les
curations TA/dark_enzymes).

Stdlib pur. Cache réutilisable sous data/proteomics/paxdb/ (téléchargé une fois).
Run: python analyses/phase19_proteomics.py
"""
from __future__ import annotations
import json, re, glob, urllib.request
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "proteomics" / "paxdb"
OUT = ROOT / "résultats" / "phase19_proteomics"
GENES = ROOT / "site" / "content" / "genes"
BASE_URL = "https://pax-db.org/downloads/latest/datasets/83332"

INTEGRATED = "83332-WHOLE_ORGANISM-integrated.txt"
# les 16 jeux individuels intégrés par PaxDb (cf. header #weights de l'intégré)
DATASETS = [
    "83332-GPM_201408.txt",
    "83332-M.tuberculosis_atlas_build_330.txt",
    "83332-Mtuberculosis_Albrethsen_2013.txt",
    "83332-My_PXD000483_Cortes_2013.txt",
    "83332-Mycobacterium-tuberculosis_PA_2013-7.txt",
    "83332-MycobacteriumTuberculosis.txt",
    "83332-PA_201307.txt",
    "83332-PXD000111_Albrethsen_Mol_Cell_Proteomics_2013_M-tuberculosis_log_phase.txt",
    "83332-PXD000111_Albrethsen_Mol_Cell_Proteomics_2013_M-tuberculosis_starvation_phase.txt",
    "83332-PXD000259_Schubert_Cell_Host_Microbe_2014_Mtuberculosis_complete_proteome.txt",
    "83332-PXD009676_late_exponential_Mtb_culture.txt",
    "83332-PXD011081_whole_cell_lysate.txt",
    "83332-PXD013677_Mtb_H37Rv_control.txt",
    "83332-PXD013677_Mtb_SAMMtb_control.txt",
    "83332-PXD016006_ctrl.txt",
    "83332-PXD020383_mtbc_lineage_pool.txt",
]
N_DATASETS = len(DATASETS)

SOURCE = "PaxDb 5.0 (M. tuberculosis H37Rv, org 83332)"
REFS = [
    {"authors": "Huang Q, Szklarczyk D, Wang M, Simonovic M, von Mering C", "year": 2023,
     "title": "PaxDb 5.0: Curated Protein Quantification Data Suggests Adaptive Proteome Changes in Yeasts",
     "journal": "Molecular & Cellular Proteomics", "doi": "10.1016/j.mcpro.2023.100640"},
    {"authors": "Schubert OT, Mouritsen J, Ludwig C, et al.", "year": 2013,
     "title": "The Mtb proteome library: a resource of assays to quantify the complete proteome of Mycobacterium tuberculosis",
     "journal": "Cell Host & Microbe", "doi": "10.1016/j.chom.2013.04.008"},
    {"authors": "Albrethsen J, Agner J, Piersma SR, et al.", "year": 2013,
     "title": "Proteomic profiling of Mycobacterium tuberculosis identifies nutrient-starvation-responsive toxin-antitoxin systems",
     "journal": "Molecular & Cellular Proteomics", "doi": "10.1074/mcp.M112.018846"},
]

RV = re.compile(r"^Rv\d{4}[A-Bc]?$")


def ensure_cache():
    CACHE.mkdir(parents=True, exist_ok=True)
    for f in [INTEGRATED] + DATASETS:
        p = CACHE / f
        if p.exists() and p.stat().st_size > 0:
            continue
        url = f"{BASE_URL}/{f}"
        print(f"  download {f} ...")
        try:
            urllib.request.urlretrieve(url, p)
        except Exception as e:
            print(f"  !! échec {f}: {e}")


def parse_paxdb(path: Path):
    """Retourne {rv: abundance_float} pour un fichier PaxDb (3 col, header #)."""
    out = {}
    if not path.exists():
        return out
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        sid = parts[1]
        rv = sid.split(".", 1)[1] if "." in sid else sid
        try:
            ab = float(parts[2])
        except ValueError:
            continue
        if ab > 0:
            out[rv] = ab
    return out


def main():
    print("== phase19 : couche protéomique PaxDb (P5.1) ==")
    ensure_cache()

    integ = parse_paxdb(CACHE / INTEGRATED)
    print(f"Intégré : {len(integ)} protéines quantifiées")

    # comptage de détection reproductible sur les 16 jeux
    det_count = Counter()
    for f in DATASETS:
        d = parse_paxdb(CACHE / f)
        for rv in d:
            det_count[rv] += 1
    print(f"Détection cumulée sur {N_DATASETS} jeux : {len(det_count)} protéines vues au moins une fois")

    # rang / percentile sur l'abondance intégrée (décroissant)
    ranked = sorted(integ.items(), key=lambda kv: kv[1], reverse=True)
    n_ranked = len(ranked)
    rank_of = {rv: i + 1 for i, (rv, _) in enumerate(ranked)}

    # univers des protéines détectées = intégré U (au moins un jeu)
    detected = set(integ) | set(det_count)

    rec = {}
    for rv in detected:
        ab = integ.get(rv)
        r = rank_of.get(rv)
        rec[rv] = {
            "detected": True,
            "abundance_ppm": round(ab, 2) if ab is not None else None,
            "abundance_rank": r,
            "n_ranked": n_ranked,
            "abundance_percentile": round(100.0 * (1 - (r - 1) / n_ranked), 1) if r else None,
            "n_datasets": det_count.get(rv, 0),
            "n_datasets_total": N_DATASETS,
            "source": SOURCE,
            "refs": REFS,
        }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "proteomics.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))
    print(f"Écrit {OUT/'proteomics.json'} ({len(rec)} protéines détectées)")

    # ── fusion directe dans les fiches (idempotent, patron phase8/9/14) ──
    n_written = n_det = n_undet = 0
    hyp_total = hyp_det = 0
    by_loc_det = Counter(); by_loc_tot = Counter()
    by_loc_ab = {}  # classe localisation -> liste d'abondances (protéines détectées)
    for path in glob.glob(str(GENES / "*.json")):
        d = json.load(open(path))
        rv = d["rv"]
        if rv in rec:
            d["proteomics"] = rec[rv]
            n_det += 1
        else:
            d["proteomics"] = {
                "detected": False, "abundance_ppm": None, "abundance_rank": None,
                "n_datasets": 0, "n_datasets_total": N_DATASETS,
                "source": SOURCE, "refs": REFS,
            }
            n_undet += 1
        Path(path).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n_written += 1

        is_hyp = "hypothetical" in (d.get("product_h37rv") or "").lower()
        if is_hyp:
            hyp_total += 1
            if rv in rec:
                hyp_det += 1
        # cross-check vs localisation DeepTMHMM
        loc = (d.get("localization") or {}).get("type")
        if loc:
            by_loc_tot[loc] += 1
            if rv in rec:
                by_loc_det[loc] += 1
                ab = rec[rv].get("abundance_ppm")
                if ab is not None:
                    by_loc_ab.setdefault(loc, []).append(ab)

    print(f"Fusionné 'proteomics' dans {n_written} fiches ({n_det} détectés, {n_undet} non détectés)")

    # ── bilan mission ──
    print("\n── Preuve d'existence MS ──")
    print(f"Protéome : {n_det}/{n_written} détectés par MS ({100*n_det/n_written:.1f} %)")
    print(f"HYPOTHÉTIQUES : {hyp_det}/{hyp_total} détectés par MS ({100*hyp_det/max(hyp_total,1):.1f} %)  <<< chiffre-phare")
    # robustesse : distribution du nombre de jeux
    dist = Counter(v["n_datasets"] for v in rec.values())
    multi = sum(c for k, c in dist.items() if k >= 5)
    print(f"Détectés dans >=5 jeux MS indépendants : {multi} protéines (détection très robuste)")

    print("\n── Cross-check vs localisation DeepTMHMM (taux de détection & abondance médiane) ──")
    import statistics
    for loc in ("GLOB", "TM", "SP", "SP+TM"):
        if by_loc_tot.get(loc):
            rate = 100 * by_loc_det[loc] / by_loc_tot[loc]
            abs_ = by_loc_ab.get(loc, [])
            med = statistics.median(abs_) if abs_ else 0
            print(f"  {loc:6s} : {by_loc_det[loc]:4d}/{by_loc_tot[loc]:4d} détectés ({rate:5.1f} %) | abondance médiane {med:8.1f} ppm")


if __name__ == "__main__":
    main()
