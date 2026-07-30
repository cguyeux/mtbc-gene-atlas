#!/usr/bin/env python3
"""phase72_pdb_delta.py -- P16.2c : delta PDB (volet « péremption » de l'audit de fraîcheur).

De nouvelles structures expérimentales ont-elles été déposées depuis notre passage (phase25) ? Question posée
en priorité sur les gènes `dark` : une structure expérimentale d'un dark serait un événement.

Interroge PDBe/SIFTS `best_structures` en POST par lot sur les accessions UniProt des gènes SANS structure.

PIÈGE MAJEUR (et sa parade) — « le silence n'est pas le succès » :
l'API renvoie **HTTP 404** quand AUCUNE accession du lot n'a de structure. Un 404 est donc ambigu : « rien
trouvé » ou « appel cassé » ? Conclure « 0 » depuis un appel cassé serait une faute.
PARADE : glisser un **TÉMOIN POSITIF** (Ag85B, P9WQP1, 8 structures connues) dans CHAQUE lot. Si le témoin
revient, l'appel a abouti et un lot sans autre hit signifie vraiment « aucune structure ». Si le témoin manque,
le lot est suspect et on ne conclut pas.

RÉSULTAT (2026-07-13) : 4/4 lots validés par le témoin → **0 des 207 gènes dark sans structure n'en a acquis une**.
Un seul des 219 dark possède une structure expérimentale, tous PDB confondus.
Converge avec les deux autres volets : littérature (187/219 jamais cités) et UniProt (0/219 recurés).

Lecture SEULE (rien à écrire : aucun hit). Run: python analyses/phase72_pdb_delta.py
"""
from __future__ import annotations
import glob, json, time, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENES = ROOT / "site" / "content" / "genes"
URL = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/"
CTRL = "P9WQP1"      # Ag85B : témoin positif, 8 structures — DOIT revenir dans chaque lot
BATCH = 60


def post(accs: list[str]) -> dict | None:
    req = urllib.request.Request(URL, data=",".join(accs).encode(), headers={"Content-Type": "text/plain"})
    try:
        return json.load(urllib.request.urlopen(req, timeout=60))
    except Exception:
        return None      # inclut le 404 « aucune donnée pour ce lot »


def main() -> None:
    todo = []
    n_with = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.loads(Path(f).read_text())
        if d.get("verdict") != "dark":
            continue
        acc = (d.get("uniprot") or {}).get("acc")
        if d.get("pdb"):
            n_with += 1
        elif acc:
            todo.append((d["rv"], acc))
    rv_of = {a: r for r, a in todo}
    print(f"dark avec structure (snapshot) : {n_with} | dark à re-interroger : {len(todo)}")

    found, ok, suspect = {}, 0, 0
    accs = [a for _, a in todo]
    for i in range(0, len(accs), BATCH):
        lot = accs[i:i + BATCH]
        d = post([CTRL] + lot)
        if not d or not d.get(CTRL):
            suspect += 1
            print(f"  lot {i//BATCH}: TÉMOIN ABSENT -> appel suspect, on NE conclut PAS sur ce lot")
            continue
        ok += 1
        for acc, hits in d.items():
            if acc != CTRL and hits:
                found[acc] = hits
        print(f"  lot {i//BATCH}: OK (témoin {len(d[CTRL])} structures) — {len(lot)} testés, cumul {len(found)}")
        time.sleep(1)

    print(f"\nLots validés par le témoin : {ok} | suspects : {suspect}")
    print(f">>> gènes DARK ayant acquis une structure expérimentale : {len(found)}")
    for acc, h in found.items():
        print(f"   {rv_of[acc]} -> PDB {h[0].get('pdb_id')} (cov {h[0].get('coverage')})")
    if not found and not suspect:
        print("    AUCUN — négatif propre et VALIDÉ (le témoin prouve que les appels ont abouti).")


if __name__ == "__main__":
    main()
