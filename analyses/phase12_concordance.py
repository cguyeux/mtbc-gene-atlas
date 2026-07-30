#!/usr/bin/env python3
"""phase12_concordance.py -- validation de la ré-annotation par concordance EC.

Question posée par le relecteur : « la couverture est haute, mais est-elle juste ? »
On mesure l'accord entre l'EC assigné par notre pipeline (eggNOG-mapper, orthologie)
et l'EC curé de UniProt/SwissProt, sur les gènes où les DEUX existent. UniProt/SwissProt
sert d'étalon-or indépendant : la couche orthologie n'y puise pas.

Trois niveaux d'accord :
  - exact          : intersection non vide des EC complets (a.b.c.d)
  - sous-sous-classe: intersection non vide des préfixes à 3 champs (a.b.c)
  - classe         : intersection non vide du 1er champ (a)

Nuance quantifiée : deux révisions récentes de la nomenclature EC produisent des
« désaccords » qui n'en sont pas. (i) En août 2018, la classe 7 (« translocases ») a été
créée en déplaçant de nombreuses entrées 3.6.3.- (ATPases de transport) et 1.6/1.7.-.
(ii) Plus récemment, la sous-classe 5.6.2 (« motor proteins acting on nucleic acids »)
a absorbé les hélicases/exonucléases ADN de la classe 3 (3.6.4.12, 3.1.11.-). Un désaccord
où exactement un côté relève de ces classes nouvelles est un artefact de version de
nomenclature, PAS un désaccord fonctionnel. On les compte à part.

Source de vérité : site/content/genes/*.json. Stdlib pure.
Sorties : résultats/phase12_concordance.json (résumé) + article/supplementary_materials/table_s5_ec_discordances.tsv
Run: python analyses/phase12_concordance.py
"""
from __future__ import annotations
import glob, json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
OUTJSON = ROOT / "résultats" / "phase12_concordance.json"
OUTTSV = ROOT / "article" / "supplementary_materials" / "table_s5_ec_discordances.tsv"


def eclist(x):
    return sorted({e for e in (x or []) if e and e != "-" and e.count(".") == 3})


def pref(ec, n):
    return ".".join(ec.split(".")[:n])


def classof(ec):
    return ec.split(".")[0]


def load():
    out = {}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        out[d["rv"]] = d
    return out


def main():
    g = load()
    rows = []
    for rv, d in g.items():
        ee = eclist((d.get("eggnog") or {}).get("ec"))
        ue = eclist((d.get("uniprot") or {}).get("ec"))
        if ee and ue:
            rows.append((rv, d.get("gene") or "", ee, ue))

    n = len(rows)
    exact = tri = cls = 0
    discord = []            # ni exact ni classe -> vrai candidat désaccord
    tri_only = []           # accord 3-niveaux, pas exact
    reclass = []            # désaccord de classe attribuable à la classe 7 (translocase)
    for rv, gene, ee, ue in rows:
        se, su = set(ee), set(ue)
        if se & su:
            exact += 1; tri += 1; cls += 1; continue
        if {pref(e, 3) for e in ee} & {pref(e, 3) for e in ue}:
            tri += 1; cls += 1; tri_only.append((rv, gene, ee, ue)); continue
        if {classof(e) for e in ee} & {classof(e) for e in ue}:
            cls += 1
            discord.append((rv, gene, ee, ue, "same_class_diff_subclass"))
            continue
        # désaccord de classe -- écarter les reclassements de nomenclature
        cl_e = {classof(e) for e in ee}; cl_u = {classof(e) for e in ue}
        sub_e = {pref(e, 2) for e in ee}; sub_u = {pref(e, 2) for e in ue}
        # (i) classe 7 (translocases, 2018) présente d'un seul côté
        is_translocase = ("7" in (cl_e | cl_u)) and ("7" not in (cl_e & cl_u))
        # (ii) sous-classe 5.6 (motor proteins, absorbe les hélicases classe 3) d'un seul côté
        motor = ("5.6" in (sub_e | sub_u)) and ("5.6" not in (sub_e & sub_u))
        is_motor = motor and ("3" in (cl_e | cl_u))
        if is_translocase:
            tag = "translocase_reclass"; reclass.append((rv, gene, ee, ue, tag))
        elif is_motor:
            tag = "motor_protein_reclass"; reclass.append((rv, gene, ee, ue, tag))
        else:
            tag = "class_mismatch"
        discord.append((rv, gene, ee, ue, tag))

    genuine = [d for d in discord if d[4] == "class_mismatch"]

    summary = {
        "genes_with_both_EC": n,
        "exact": exact, "exact_pct": round(100 * exact / n, 1),
        "three_level": tri, "three_level_pct": round(100 * tri / n, 1),
        "class_level": cls, "class_level_pct": round(100 * cls / n, 1),
        "class_mismatches_total": n - cls,
        "nomenclature_reclassification": len(reclass),
        "genuine_class_mismatch": len(genuine),
        "genuine_class_mismatch_pct": round(100 * len(genuine) / n, 1),
    }
    print(json.dumps(summary, indent=2))
    print("\nExemples reclassement translocase (classe 7):")
    for rv, gene, ee, ue, _ in reclass[:8]:
        print(f"  {rv} {gene}: eggNOG {ee} vs UniProt {ue}")
    print("\nExemples désaccord de classe résiduel (genuine):")
    for rv, gene, ee, ue, _ in genuine[:12]:
        print(f"  {rv} {gene}: eggNOG {ee} vs UniProt {ue}")

    OUTJSON.parent.mkdir(parents=True, exist_ok=True)
    OUTJSON.write_text(json.dumps(summary, indent=2))
    OUTTSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTTSV, "w") as fh:
        fh.write("locus\tgene\teggNOG_EC\tUniProt_EC\tcategory\n")
        for rv, gene, ee, ue, tag in sorted(discord):
            fh.write(f"{rv}\t{gene}\t{';'.join(ee)}\t{';'.join(ue)}\t{tag}\n")
    print(f"\nWrote {OUTJSON} and {OUTTSV} ({len(discord)} discordant rows)")


if __name__ == "__main__":
    main()
