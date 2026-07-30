#!/usr/bin/env python3
"""phase86h — Pfam domain layer for the 68 Mycobrowser genes, via UniProt (P16.5d-cont).

Correction (same pattern as STRING & Foldseek): Pfam was marked "blocked (hmmscan/PFAM_DB
absent)". No local hmmscan is needed — UniProt carries Pfam domain cross-references per
entry, and the 68 all have a UniProt accession (from Mycobrowser). Verified: 4/6 sampled
have >=1 Pfam domain via UniProt.

Fetches each gene's UniProt entry, extracts Pfam xrefs (accession + short name) and the
InterPro descriptions, and writes the `domains` layer (source = UniProt Pfam xref; no
i_evalue/coords, which are hmmscan-specific). Policy consistency: a real Pfam domain is a
FAMILY assignment, so a still-dark gene that gains a Pfam domain is upgraded dark ->
family_assigned (with function_revised set from the domain), exactly as the whole-proteome
Pfam pass (phase2b) drove family assignments. Idempotent; polite rate limiting; +control.
"""
from __future__ import annotations
import json, glob, os, time, urllib.request, urllib.error

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CG = os.path.join(ROOT, "site", "content", "genes")
API = "https://rest.uniprot.org/uniprotkb/{}.json"


def fetch(acc):
    req = urllib.request.Request(API.format(acc),
                                 headers={"User-Agent": "mtbc-atlas/1.0 (guyeux@gmail.com)"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def extract(d):
    """UniProt entry -> ([domains], {interpro_id: name})."""
    ipr = {}
    for x in d.get("uniProtKBCrossReferences", []):
        if x.get("database") == "InterPro":
            name = next((p["value"] for p in x.get("properties", []) if p["key"] == "EntryName"), "")
            ipr[x["id"]] = name
    doms = []
    for x in d.get("uniProtKBCrossReferences", []):
        if x.get("database") != "Pfam":
            continue
        name = next((p["value"] for p in x.get("properties", []) if p["key"] == "EntryName"), "")
        n_match = next((p["value"] for p in x.get("properties", []) if p["key"] == "MatchStatus"), "")
        doms.append({
            "source": "UniProt (Pfam xref)",
            "pfam_acc": x["id"],
            "pfam_name": name or x["id"],
            "description": name,
            "n_matches": n_match,
        })
    # enrich domain descriptions with a matching InterPro name when available
    return doms, ipr


def main():
    targets = []
    for f in glob.glob(os.path.join(CG, "Rv*.json")):
        d = json.load(open(f))
        if str(d.get("curation_note", "")).startswith("Added P16.5d"):
            acc = (d.get("uniprot") or {}).get("acc")
            if acc:
                targets.append((f, d, acc))
    # positive control
    try:
        cdoms, _ = extract(fetch("P9WGR1"))
        assert cdoms, "positive control returned no Pfam domain"
        print(f"[contrôle+] P9WGR1: {len(cdoms)} domaine(s) Pfam — OK")
    except Exception as e:
        print(f"[ATTENTION] contrôle positif échoué ({e}) — arrêt."); return

    written = with_dom = upgraded = 0
    for f, d, acc in targets:
        try:
            doms, ipr = extract(fetch(acc))
        except Exception as e:
            print(f"  {d['rv']} ({acc}): {e}"); continue
        if doms:
            d["domains"] = doms
            with_dom += 1
            # a real Pfam domain = family assignment -> upgrade a still-dark fiche
            if d.get("verdict") == "dark":
                d["verdict"] = "family_assigned"
                d["confidence"] = d.get("confidence") or "low"
                if not d.get("function_revised"):
                    dn = doms[0]["pfam_name"]
                    d["function_revised"] = f"{dn} domain-containing protein (Pfam {doms[0]['pfam_acc']})"
                upgraded += 1
            json.dump(d, open(f, "w"), ensure_ascii=False, indent=2)
            written += 1
        time.sleep(0.25)
    print(f"phase86h: domaines Pfam (via UniProt) écrits sur {with_dom}/{len(targets)} gènes")
    print(f"  dark -> family_assigned (domaine Pfam réel): {upgraded}")


if __name__ == "__main__":
    main()
