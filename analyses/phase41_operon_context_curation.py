#!/usr/bin/env python3
"""phase41_operon_context_curation.py -- P7.10 : requalification par CONTEXTE d'opéron.

Généralise P7.10a (co-transcription > STRING) à TOUT le stock dark : requalifie les dark
co-transcrits (opéron réel, couche `genomic_context` = distance intergénique co-directionnelle)
avec un voisin NOMMÉ de VOIE métabolique connue, quand la localisation et la conservation
concordent. Complémentaire de P5.2b (phenotype_lead), qui couvrait les dark à PHÉNOTYPE ;
ici la valeur ajoutée = les dark SANS phénotype (88/110 du gisement).

Garde-fous (KB P7.10a) : « co-transcrit avec X » ≠ « est X » → formulation « candidate
<voie>-associated », fonction non figée, confidence low ; exiger cohérence localisation
(membranaire pour un transporteur) + sélection purifiante (gène réel). Voisins non caractérisés
(RvXXXX sans fonction) ou fonction trop centrale/générique (lpd) ÉCARTÉS.

Émet aussi un rapport `operon_context_leads.tsv` du gisement complet (110 dark co-transcrits avec
un voisin nommé) pour instruction ultérieure. Couche `context_lead`, hand-review `auto:false`.
Run: python analyses/phase41_operon_context_curation.py
"""
from __future__ import annotations
import json, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "résultats" / "phase41_operon_context"

CURATIONS = {
    "Rv2219": {"operon": "lipB-lipA-Rv2219", "pathway": "lipoyl-cofactor biosynthesis",
        "function_revised": (
            "Membrane protein (2 predicted TM helices) co-transcribed with lipB and lipA — the "
            "octanoyltransferase and lipoyl synthase of lipoic-acid cofactor biosynthesis — in the "
            "lipB-lipA-Rv2219 operon, under purifying selection. Contextual candidate: a component "
            "associated with lipoyl-cofactor biosynthesis/utilisation. Its molecular role is not established.")},
    "Rv1231c": {"operon": "Rv1231c-mgtE", "pathway": "magnesium transport",
        "function_revised": (
            "Membrane protein (2 predicted TM helices) co-transcribed with mgtE (a CBS-domain Mg2+ "
            "transporter) as a two-gene operon, under purifying selection. Contextual candidate: associated "
            "with magnesium transport / homeostasis. Its molecular role is not established.")},
    "Rv0970": {"operon": "ctpV-Rv0970", "pathway": "metal (copper) transport",
        "function_revised": (
            "Polytopic membrane protein (6 predicted TM helices) co-transcribed with ctpV, a P-type "
            "metal-transporting ATPase (copper), under purifying selection. Contextual candidate: associated "
            "with metal (copper) transport / homeostasis. Its molecular role is not established.")},
    "Rv0633c": {"operon": "echA3-Rv0633c", "pathway": "fatty-acid beta-oxidation",
        "function_revised": (
            "Membrane protein co-transcribed with echA3 (enoyl-CoA hydratase, fatty-acid beta-oxidation), "
            "under purifying selection. Contextual candidate: associated with fatty-acid / lipid "
            "beta-oxidation metabolism. Its molecular role is not established.")},
    "Rv2767c": {"operon": "Rv2766c-Rv2767c", "pathway": "fatty-acid synthesis",
        "function_revised": (
            "Conserved protein co-transcribed with Rv2766c, an enoyl-(acyl-carrier-protein) reductase of "
            "fatty-acid synthesis, under purifying selection. Contextual candidate: associated with "
            "fatty-acid metabolism. Its molecular role is not established.")},
    "Rv3489": {"operon": "Rv3489-otsA", "pathway": "trehalose biosynthesis",
        "function_revised": (
            "Conserved protein co-transcribed with otsA (trehalose-6-phosphate synthase) under purifying "
            "selection. Contextual candidate: associated with trehalose biosynthesis / metabolism. Its "
            "molecular role is not established.")},
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    G = {}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f)); G[d["rv"]] = (f, d)

    print("== phase41 : requalification par contexte d'opéron (P7.10) ==")
    n = 0
    for rv, c in CURATIONS.items():
        f, d = G[rv]
        d["verdict"] = "family_assigned"
        d["confidence"] = "low"
        d["function_revised"] = c["function_revised"]
        d["auto"] = False
        d["needs_review"] = False
        d["context_lead"] = {"basis": f"co-transcription in operon {c['operon']} + purifying selection + coherent localisation",
                             "operon": c["operon"], "pathway": c["pathway"],
                             "source": "operon-context curation (P7.10): co-transcription with a named pathway gene"}
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n += 1
        print(f"  {rv} -> family_assigned/low [{c['operon']}] :: {c['pathway']}")

    # rapport du gisement restant (dark en opéron avec voisin nommé, non requalifié ici)
    def named_fn(rv):
        e = (G.get(rv, (None, {}))[1]).get("eggnog") or {}
        u = (G.get(rv, (None, {}))[1]).get("uniprot") or {}
        nm = e.get("preferred_name") or ""
        desc = (e.get("description") or u.get("protein_name") or "")
        if (G.get(rv, (None, {}))[1]).get("verdict") == "dark": return None
        if nm and nm not in ("", "-"): return f"{nm}|{desc[:40]}"
        if desc and "hypothetical" not in desc.lower() and "unknown" not in desc.lower(): return desc[:45]
        return None

    rows = []
    for rv, (f, d) in G.items():
        if d.get("verdict") != "dark":
            continue
        op = (d.get("genomic_context") or {}).get("operon") or []
        if len(op) < 2:
            continue
        named = [(m.get("locus"), named_fn(m.get("locus"))) for m in op if m.get("locus") != rv]
        named = [(l, fn) for l, fn in named if fn]
        if not named:
            continue
        cons = (d.get("conservation") or {}).get("selection", "")
        loc = (d.get("localization") or {}).get("type", "")
        pheno = "Y" if d.get("mutant_phenotypes") else "n"
        rows.append((rv, loc, cons, pheno, "; ".join(f"{l}:{fn}" for l, fn in named[:2])))
    rows.sort(key=lambda r: (("purifying" not in r[2]), r[0]))
    with open(OUT / "operon_context_leads.tsv", "w") as fh:
        fh.write("rv\tlocalization\tconservation\tphenotype\tnamed_cotranscribed_neighbours\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"{n} requalifié(s). Gisement restant documenté : {len(rows)} dark co-transcrits avec voisin nommé "
          f"-> {OUT/'operon_context_leads.tsv'}")


if __name__ == "__main__":
    main()
