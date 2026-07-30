#!/usr/bin/env python3
"""phase86c — UniProt curated layer for the 68 Mycobrowser genes (P16.5d-cont).

The 68 genes added in phase86b each carry a Mycobrowser UniProt accession. This
queries the UniProt REST API by accession and writes the `uniprot` layer (same
shape as phase2e): reviewed/SwissProt status, protein name, EC, curated function,
GO, KEGG, protein-existence level. Homology-independent, high value, feasible now
(no heavy toolchain). Idempotent; polite rate limiting; positive-control check.
"""
from __future__ import annotations
import json, glob, os, time, urllib.request, urllib.error

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CG = os.path.join(ROOT, "site", "content", "genes")
API = "https://rest.uniprot.org/uniprotkb/{}.json"
REF = {"authors": "The UniProt Consortium", "year": 2023,
       "title": "UniProt: the Universal Protein Knowledgebase in 2023",
       "journal": "Nucleic Acids Res 51:D523-D531"}


def fetch(acc):
    req = urllib.request.Request(API.format(acc),
                                 headers={"User-Agent": "mtbc-atlas/1.0 (guyeux@gmail.com)"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def parse(d):
    desc = d.get("proteinDescription", {})
    rec = desc.get("recommendedName") or (desc.get("submissionNames") or [{}])[0]
    name = (rec.get("fullName", {}) or {}).get("value")
    ec = [e["value"] for e in (rec.get("ecNumbers") or [])]
    func = None
    for c in d.get("comments", []):
        if c.get("commentType") == "FUNCTION":
            texts = c.get("texts") or []
            if texts:
                func = texts[0].get("value")
            break
    go, kegg = [], []
    for x in d.get("uniProtKBCrossReferences", []):
        if x.get("database") == "GO":
            go.append(x["id"])
        elif x.get("database") == "KEGG":
            kegg.append(x["id"])
    return {
        "acc": d.get("primaryAccession"),
        "reviewed": d.get("entryType", "").startswith("UniProtKB reviewed"),
        "protein_name": name,
        "ec": ec,
        "function": func,
        "go": go,
        "kegg": kegg,
        "protein_existence": d.get("proteinExistence"),
        "source": "UniProt REST API (curated)",
        "refs": [REF],
    }


def main():
    targets = []
    for f in glob.glob(os.path.join(CG, "Rv*.json")):
        d = json.load(open(f))
        if str(d.get("curation_note", "")).startswith("Added P16.5d"):
            acc = (d.get("mycobrowser") or {}).get("uniprot_ac")
            if acc:
                targets.append((f, d, acc.split(";")[0].strip()))
    # positive control: a well-known accession must return a function
    try:
        ctrl = parse(fetch("P9WGR1"))       # KatG
        assert ctrl.get("function"), "positive control returned no function"
        print(f"[contrôle+] P9WGR1 KatG: {ctrl['protein_name']} — OK")
    except Exception as e:
        print(f"[ATTENTION] contrôle positif échoué ({e}); l'API UniProt est peut-être cassée — arrêt.")
        return

    n_func = n_ec = written = 0
    for f, d, acc in targets:
        try:
            layer = parse(fetch(acc))
        except urllib.error.HTTPError as e:
            print(f"  {d['rv']} ({acc}): HTTP {e.code}")
            continue
        except Exception as e:
            print(f"  {d['rv']} ({acc}): {e}")
            continue
        d["uniprot"] = layer
        # if UniProt gives a real function and the fiche is still dark, upgrade it (curated evidence)
        if d.get("verdict") == "dark" and (layer.get("function") or layer.get("ec")):
            d["verdict"] = "family_assigned"
            d["confidence"] = d.get("confidence") or "low"
            if not d.get("function_revised") and layer.get("function"):
                d["function_revised"] = layer["function"]
        json.dump(d, open(f, "w"), ensure_ascii=False, indent=2)
        written += 1
        n_func += bool(layer.get("function"))
        n_ec += bool(layer.get("ec"))
        time.sleep(0.3)
    print(f"phase86c-uniprot: couche uniprot écrite sur {written}/{len(targets)} gènes")
    print(f"  avec fonction curée: {n_func} | avec EC: {n_ec}")


if __name__ == "__main__":
    main()
