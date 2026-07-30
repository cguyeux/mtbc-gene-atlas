#!/usr/bin/env python3
"""phase81 — ABPP active serine-hydrolase layer (P16.12), census-based.

Source: Li M, Patel HV, Cognetta AB 3rd, ..., Canaan S, Aldridge BB et al.
"Identification of cell wall synthesis inhibitors active against Mycobacterium
tuberculosis by competitive activity-based protein profiling." Cell Chem Biol
2021 (PMC8964833; doi:10.1016/j.chembiol.2021.09.002). OPEN ACCESS.

Two data sources, merged per gene:
  - Table S3 "Active serine hydrolase proteome" (supplementary mmc4.xlsx,
    downloaded via the els-cdn OA path because the PMC bin/ path serves HTML):
    the CENSUS of proteins enriched by a fluorophosphonate (FP) serine-hydrolase
    activity probe, with the Probe/No-probe enrichment and the competition ratio
    against covalent serine-hydrolase inhibitors (Cmpd/Vehicle). This is the broad
    list CG asked for (~44 genes) and includes many still-vaguely-annotated ones.
  - Table 2 (inline main text, from the OA full-text XML): 18 high-confidence
    PRIORITIZED serine hydrolase targets of AA692, with condition dependence
    (active at pH 6.6 vs pH 5.0, under hypoxia) and covalent-inhibitor sensitivity.

Serine-hydrolase criterion (per gene, conservative, anti-oversell): FP-reactive AND
[ competed by a covalent inhibitor (max Cmpd/Vehicle >= 2) OR annotated as a
hydrolase-family enzyme OR strongly probe-enriched (Probe/No-probe = 20 cap) ],
EXCLUDING proteins annotated as a canonical NON-hydrolase (oxidoreductase,
dehydrogenase, methyltransferase, regulator, lyase) that are NOT competed — those
are background. CAVEAT retained in the layer: covalent probes/inhibitors can also
label active-site cysteines, so a competed protein annotated as a non-hydrolase
(e.g. an oxidoreductase) is flagged DISCORDANT (likely covalent off-target), and
the atlas assignment is retained rather than overwritten.

Atlas concordance: for each census gene we compare with the atlas' own
`function_revised`. Measured result: 0 dark, 0 hypothetical in the atlas; ~24
carry a still-vague H37Rv product ("hypothetical protein"/"hydrolase"/"membrane
protein") yet the atlas already assigned a serine-hydrolase family independently
-> ABPP is EXPERIMENTAL CORROBORATION of those computational calls, not a
dark-matter mover. This is logged as `concordance: confirms`.

Writes layer `abpp` in place into site/content/genes/Rv*.json. Read-only on
verdicts. Idempotent (clears any previous `abpp` before writing).
"""
from __future__ import annotations
import xml.etree.ElementTree as ET
import json, re, os, glob, warnings
warnings.filterwarnings("ignore")
from openpyxl import load_workbook

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ABPP = os.path.join(ROOT, "résultats", "phase81_abpp")
XML = os.path.join(ABPP, "pmc8964833.xml")
S3 = os.path.join(ABPP, "census", "li_mmc4.xlsx")
CG = os.path.join(ROOT, "site", "content", "genes")
CITATION = ("Li M, Patel HV, ..., Canaan S, Aldridge BB et al., Cell Chem Biol 2021 "
            "(PMC8964833; doi:10.1016/j.chembiol.2021.09.002); competitive "
            "activity-based protein profiling of FP-reactive serine hydrolases")

RV = re.compile(r"Rv\d{4}[A-Za-z]?")
HYD = re.compile(r"hydrolase|esterase|lipase|protease|cutinase|thioesterase|"
                 r"peroxidase|amidase|phospholip|lactamase|NTE family|peptidase|"
                 r"patatin|serine hydrolase", re.I)
NONH = re.compile(r"oxidoreductase|dehydrogenase|methyltransfer|transcription|"
                  r"beta-lyase|phosphoribosyl|\breductase\b|synthase|kinase", re.I)
CATNAME = {"CH": "conserved hypotheticals", "IMR": "intermediary metabolism & respiration",
           "LM": "lipid metabolism", "CWP": "cell wall & cell processes",
           "RP": "regulatory proteins", "IP": "information pathways",
           "VAPBS": "virulence/detox/adaptation", "IA": "insertion seqs & phages"}
COMP_COLS = [4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19]


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


# ---- Table 2 (inline, condition-rich, 18 genes) ----
def parse_table2():
    root = ET.parse(XML).getroot()
    COND = ["pH 6.6", "pH 5.0"]
    FLAGS = ["in_vitro_essential", "in_vivo_essential", "active_hypoxia"]
    INHIB = ["EZ120", "Lalistat", "CyC17", "THL"]
    out = {}
    for tw in root.iter():
        if not tw.tag.endswith("table-wrap"):
            continue
        lab = next((("".join(c.itertext())).strip() for c in tw if c.tag.endswith("label")), "")
        if lab != "Table 2":
            continue
        for tr in tw.iter():
            if not tr.tag.endswith("tr"):
                continue
            cells = ["".join(td.itertext()).strip() for td in tr if td.tag.endswith(("td", "th"))]
            if len(cells) >= 3 and RV.fullmatch(cells[1]):
                c = cells + [""] * (12 - len(cells))
                out[cells[1]] = {
                    "prioritized_at_pH": [COND[i] for i in (0, 1) if "●" in c[3 + i]],
                    **{FLAGS[i]: ("●" in c[5 + i]) for i in range(3)},
                    "covalent_inhibitor_targets": [INHIB[i] for i in range(4) if "●" in c[8 + i]],
                }
    return out


# ---- Table S3 (census: enrichment + competition) ----
def parse_s3():
    ws = load_workbook(S3, data_only=True)["Table S3"]
    rows = list(ws.iter_rows(values_only=True))
    out = {}
    for r in rows[7:]:
        if not r[1]:
            continue
        m = RV.search(str(r[21] or r[1]))
        if not m:
            continue
        rv = m.group(0)
        pnp = fnum(r[2])
        comps = [v for c in COMP_COLS if (v := fnum(r[c])) is not None]
        out[rv] = {
            "probe_enrichment": pnp,                    # Probe/No-probe (20 = capped max)
            "max_inhibitor_competition": max(comps) if comps else None,
            "source_annotation": str(r[1]).strip(),
            "functional_category": r[22],
        }
    return out


def sh_criterion(rec):
    """Conservative serine-hydrolase call from S3 fields. Returns (is_sh, is_nonhyd_competed)."""
    ann = rec["source_annotation"]
    pnp = rec["probe_enrichment"]
    mc = rec["max_inhibitor_competition"]
    ish = bool(HYD.search(ann))
    inonh = bool(NONH.search(ann))
    competed = (mc is not None and mc >= 2)
    if inonh and not ish:
        # non-hydrolase annotation: only a covalent-off-target flag, not a SH call
        return (False, competed)
    is_sh = ish or competed or (pnp == 20)
    return (is_sh, False)


def main():
    t2 = parse_table2()
    s3 = parse_s3()
    genes = sorted(set(t2) | set(s3))

    # clear any previous abpp layer
    for f in glob.glob(os.path.join(CG, "Rv*.json")):
        d = json.load(open(f))
        if d.get("abpp") is not None:
            d.pop("abpp", None)
            json.dump(d, open(f, "w"), ensure_ascii=False, indent=2)

    written, confirms, discord, offtarget = 0, 0, 0, []
    for rv in genes:
        fn = os.path.join(CG, f"{rv}.json")
        if not os.path.exists(fn):
            continue
        d = json.load(open(fn))
        rec3 = s3.get(rv)
        rec2 = t2.get(rv)

        is_sh, nonhyd_competed = (True, False)
        if rec3:
            is_sh, nonhyd_competed = sh_criterion(rec3)
        # a Table-2 prioritized target is by construction a serine hydrolase
        if rec2:
            is_sh = True

        # atlas concordance vs its own family call
        fr = (d.get("function_revised") or "")
        atlas_is_hyd = bool(HYD.search(fr))
        atlas_is_nonhyd = bool(NONH.search(fr)) and not atlas_is_hyd
        # DISCORDANT when the FP hit contradicts a non-hydrolase family call — from the
        # paper's annotation (nonhyd_competed) OR the atlas' own call (atlas_is_nonhyd).
        # Covalent probes/inhibitors also label active-site cysteines, so this is treated
        # as a likely off-target and the atlas assignment (structure/homology-based) wins.
        if nonhyd_competed or atlas_is_nonhyd:
            concordance = "discordant"
            is_sh = False
            offtarget.append(rv)
            discord += 1
        elif atlas_is_hyd:
            concordance = "confirms"          # atlas already called a serine-hydrolase family
            confirms += 1
        else:
            concordance = "consistent"

        layer = {
            "fp_active_serine_hydrolase": bool(is_sh and concordance != "discordant"),
            "source_annotation": (rec3 or {}).get("source_annotation"),
            "functional_category": CATNAME.get(str((rec3 or {}).get("functional_category") or ""),
                                               (rec3 or {}).get("functional_category")),
            "probe_enrichment": (rec3 or {}).get("probe_enrichment"),
            "max_inhibitor_competition": (rec3 or {}).get("max_inhibitor_competition"),
            "in_census": rv in s3,
            "prioritized_target": rv in t2,
            "atlas_family_call": fr[:80] if fr else None,
            "concordance": concordance,
            "source": CITATION,
        }
        if rec2:
            layer["prioritized_at_pH"] = rec2["prioritized_at_pH"]
            layer["active_hypoxia"] = rec2["active_hypoxia"]
            layer["in_vitro_essential"] = rec2["in_vitro_essential"]
            layer["in_vivo_essential"] = rec2["in_vivo_essential"]
            layer["covalent_inhibitor_targets"] = rec2["covalent_inhibitor_targets"]
        if concordance == "discordant":
            layer["note"] = ("FP-reactive / inhibitor-competed but annotated as a non-hydrolase; "
                             "covalent probes can also label active-site cysteines, so this is flagged "
                             "as a likely off-target rather than a serine-hydrolase call. Atlas assignment retained.")
        elif concordance == "confirms":
            layer["note"] = ("experimental (activity-based) confirmation of the atlas serine-hydrolase "
                             "family assignment; proves the enzyme is catalytically active in live M. tuberculosis.")
        else:
            layer["note"] = ("FP-reactive serine hydrolase (activity-based protein profiling); "
                             "proves active-enzyme status, does not assign a physiological substrate.")

        d["abpp"] = layer
        json.dump(d, open(fn, "w"), ensure_ascii=False, indent=2)
        written += 1

    print(f"phase81 (census): couche abpp écrite sur {written} gènes "
          f"(union Table S3 census [{len(s3)}] + Table 2 prioritaires [{len(t2)}])")
    print(f"  concordance atlas : confirms={confirms} (corroboration expérimentale), "
          f"discordant/off-target={discord} {offtarget}")
    # dark/hypo overlap
    dark = hypo = 0
    for rv in genes:
        fn = os.path.join(CG, f"{rv}.json")
        if os.path.exists(fn):
            d = json.load(open(fn))
            if d.get("verdict") == "dark":
                dark += 1
            if d.get("is_hypothetical"):
                hypo += 1
    print(f"  apport matière noire : dark={dark}, hypothetical={hypo} "
          f"(0/0 attendu = validation externe, pas mouvement dark)")


if __name__ == "__main__":
    main()
