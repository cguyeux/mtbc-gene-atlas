#!/usr/bin/env python3
"""phase25_pdb.py -- structures PDB EXPÉRIMENTALES résolues (P5.8).

Champ de parité Mycobrowser : les structures cristallographiques/cryo-EM réellement résolues
pour un gène, en complément de nos modèles PRÉDITS (ESMFold, AlphaFold). Une structure
expérimentale est une preuve d'existence + le meilleur support structural possible.

Source : PDBe SIFTS via l'API `mappings/best_structures/{acc}` (mapping UniProt->PDB curé),
interrogée en POST par LOT (l'endpoint accepte une liste d'accessions séparées par virgule).
Accession UniProt = champ `uniprot.acc` des fiches (3883/3906 fiches en ont une).

Par gène : structures uniques (dédupliquées par pdb_id), méthode, résolution, couverture de
la séquence, triées par couverture puis résolution.

Sortie : résultats/phase25_pdb/pdb.json (keyed Rv). Fusion directe dans les fiches (idempotent)
+ enregistrement dans phase4. Cache brut réutilisable raw.json. Stdlib.
Run: python analyses/phase25_pdb.py
"""
from __future__ import annotations
import json, glob, time, urllib.request, urllib.error
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "phase25_pdb"
GENES = ROOT / "site" / "content" / "genes"
API = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/"
BATCH = 150
CAP = 8  # structures affichées max par gène

SOURCE = "PDBe / SIFTS (experimental structures mapped from UniProt)"
REFS = [
    {"authors": "Dana JM, Gutmanas A, Tyagi N, et al.", "year": 2019,
     "title": "SIFTS: updated Structure Integration with Function, Taxonomy and Sequences resource allows 40-fold increase in coverage of structure-based annotations for proteins",
     "journal": "Nucleic Acids Research", "doi": "10.1093/nar/gky1114"},
]


def post_batch(accs):
    data = ",".join(accs).encode()
    req = urllib.request.Request(API, data=data,
                                 headers={"Content-Type": "application/x-www-form-urlencoded"})
    for _ in range(3):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return {}   # aucune des accessions n'a de structure
            time.sleep(2)
        except Exception:
            time.sleep(2)
    return {}


def main():
    print("== phase25 : structures PDB expérimentales via PDBe/SIFTS (P5.8) ==")
    OUT.mkdir(parents=True, exist_ok=True)

    # acc -> [rv...] (une accession peut être partagée par des paralogues)
    acc_to_rv = {}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        a = (d.get("uniprot") or {}).get("acc")
        if a:
            acc_to_rv.setdefault(a, []).append(d["rv"])
    accs = sorted(acc_to_rv)
    print(f"{len(accs)} accessions UniProt distinctes à interroger")

    raw_file = OUT / "raw.json"
    raw = json.loads(raw_file.read_text()) if raw_file.exists() else {}
    todo = [a for a in accs if a not in raw]
    print(f"{len(todo)} à télécharger ({len(raw)} en cache)")
    for i in range(0, len(todo), BATCH):
        chunk = todo[i:i + BATCH]
        res = post_batch(chunk)
        for a in chunk:
            raw[a] = res.get(a, [])   # [] = pas de structure (mémorisé pour ne pas re-tenter)
        print(f"  lot {i//BATCH+1}/{(len(todo)+BATCH-1)//BATCH} ({len(chunk)} acc)")
        raw_file.write_text(json.dumps(raw))

    # condenser par accession -> structures uniques par pdb_id
    def condense(structs):
        by_pdb = {}
        for s in structs:
            pid = s.get("pdb_id")
            if not pid:
                continue
            cov = s.get("coverage") or 0
            cur = by_pdb.get(pid)
            if cur is None or cov > cur["coverage"]:
                by_pdb[pid] = {
                    "pdb_id": pid,
                    "method": s.get("experimental_method") or "",
                    "resolution": s.get("resolution"),
                    "coverage": round(cov, 2),
                }
        out = sorted(by_pdb.values(),
                     key=lambda x: (-x["coverage"], x["resolution"] if x["resolution"] else 99))
        return out

    rec = {}
    for a, structs in raw.items():
        cond = condense(structs)
        if not cond:
            continue
        for rv in acc_to_rv.get(a, []):
            rec[rv] = {
                "n_structures": len(cond),
                "structures": cond[:CAP],
                "acc": a,
                "source": SOURCE, "refs": REFS,
            }

    (OUT / "pdb.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))
    print(f"Écrit {OUT/'pdb.json'} ({len(rec)} gènes avec structure expérimentale)")

    # fusion dans les fiches
    n_written = n_hit = hyp_hit = 0
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        d["pdb"] = rec.get(d["rv"], {})
        Path(f).write_text(json.dumps(d, indent=2, ensure_ascii=False))
        n_written += 1
        if d["rv"] in rec:
            n_hit += 1
            if "hypothetical" in (d.get("product_h37rv") or "").lower():
                hyp_hit += 1
    print(f"Fusionné 'pdb' dans {n_written} fiches ({n_hit} avec structure PDB, dont {hyp_hit} hypothétiques)")
    tot = sum(r["n_structures"] for r in rec.values())
    print(f"Total structures uniques cataloguées : {tot}")


if __name__ == "__main__":
    main()
