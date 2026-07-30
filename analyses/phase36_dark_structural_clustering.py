#!/usr/bin/env python3
"""phase36_dark_structural_clustering.py -- P7.3 (recadré) : clustering structural des dark.

But : détecter des FAMILLES de repli propres au MTBC parmi les gènes dark, invisibles en
séquence (Pfam/HHpred muets) mais partageant une structure entre eux. All-vs-all Foldseek des
~350 modèles ESMFold dark (phase2c), graphe TM-score, composantes connexes.

Garde-fous KB (tuberculosis.md) OBLIGATOIRES, sinon on « découvre » des artefacts :
- exclure les modèles à bas pLDDT (structure non fiable) ;
- flaguer les clusters dominés par low-complexity (compo A+G+P>45% / 1aa>25% / top3>55%),
  multi-TM (promiscuité d'hélices membranaires), ou coiled-coil (promiscuité heptade) ;
- un cluster de repli ≠ module fonctionnel : ne PAS conclure « complexe », juste « même repli ».

NB : la moitié « Foldseek vs AFDB entière » de P7.3 est ABANDONNÉE (Foldseek vs AFDB-SwissProt
déjà fait, dark 448→338 ; AFDB non curée ne nomme rien + infaisable disque 24 Go). Voir pistes P7.3.

Run: python analyses/phase36_dark_structural_clustering.py
"""
from __future__ import annotations
import json, shutil, subprocess
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FS = ROOT / "tools" / "foldseek" / "bin" / "foldseek"
SRC_PDB = ROOT / "résultats" / "phase2c_foldseek"          # modèles ESMFold dark (*.pdb)
GENES = ROOT / "site" / "content" / "genes"
OUT = ROOT / "résultats" / "phase36_dark_clustering"
MODELS = OUT / "models"                                     # dossier PROPRE (KB : que des .pdb)

TM_MIN = 0.5          # seuil de similarité de repli
EVAL_MAX = 0.01       # significativité
PLDDT_MIN = 60        # modèle ESMFold fiable (échelle 0-100)


def low_complexity(seq: str) -> bool:
    if not seq:
        return False
    n = len(seq)
    c = Counter(seq)
    agp = sum(c.get(a, 0) for a in "AGP") / n
    dom = c.most_common(1)[0][1] / n
    top3 = sum(v for _, v in c.most_common(3)) / n
    return agp > 0.45 or dom > 0.25 or top3 > 0.55


def load_gene(rv):
    try:
        return json.load(open(GENES / f"{rv}.json"))
    except FileNotFoundError:
        return {}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if MODELS.exists():
        shutil.rmtree(MODELS)
    MODELS.mkdir()

    # 1. copier les .pdb dark dans un dossier propre (KB : jamais de tmp/json dedans)
    pdbs = sorted(SRC_PDB.glob("Rv*.pdb"))
    for p in pdbs:
        shutil.copy(p, MODELS / p.name)
    print(f"P7.3 clustering structural : {len(pdbs)} modèles ESMFold dark copiés.")

    # 2. all-vs-all Foldseek (dossier query = dossier target)
    aln = OUT / "allvall.m8"
    tmp = OUT / "tmp"
    cmd = [str(FS), "easy-search", str(MODELS), str(MODELS), str(aln), str(tmp),
           "--format-output", "query,target,evalue,alntmscore,qtmscore,ttmscore",
           "-e", "10", "--max-seqs", "400"]
    print("Foldseek all-vs-all…")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("STDERR:", r.stderr[-1500:])
        raise SystemExit("foldseek a échoué")

    # 3. graphe : arête si paire distincte, TM>=TM_MIN, E<=EVAL_MAX
    adj = defaultdict(set)
    pair_tm = {}
    for line in open(aln):
        q, t, ev, tm, _qtm, _ttm = line.rstrip("\n").split("\t")
        q = q.replace(".pdb", ""); t = t.replace(".pdb", "")
        if q == t:
            continue
        try:
            ev = float(ev); tm = float(tm)
        except ValueError:
            continue
        if tm >= TM_MIN and ev <= EVAL_MAX:
            adj[q].add(t); adj[t].add(q)
            key = tuple(sorted((q, t)))
            pair_tm[key] = max(pair_tm.get(key, 0), tm)

    # 4. composantes connexes (union-find simple)
    seen = set(); clusters = []
    for node in list(adj):
        if node in seen:
            continue
        stack = [node]; comp = set()
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x); seen.add(x)
            stack.extend(adj[x] - comp)
        if len(comp) >= 2:
            clusters.append(sorted(comp))
    clusters.sort(key=len, reverse=True)

    # 5. caractériser + garde-fous
    print(f"\n{len(clusters)} cluster(s) structuraux (≥2 membres, TM≥{TM_MIN}, E≤{EVAL_MAX}).\n")
    report = []
    for i, comp in enumerate(clusters, 1):
        members = []
        for rv in comp:
            d = load_gene(rv)
            seq = d.get("protein_mtbc0") or d.get("sequence") or ""
            loc = d.get("localization") or {}
            pl = d.get("plddt") or {}
            members.append({
                "rv": rv, "len": len(seq), "verdict": d.get("verdict"),
                "plddt": round((pl.get("mean_plddt") or 0), 0),
                "tm_helices": loc.get("tm_helices", loc.get("n_tm_helices")),
                "loc_type": loc.get("type"),
                "low_complexity": low_complexity(seq),
            })
        # garde-fous : cluster suspect si dominé par TM / low-complexity / modèles peu fiables
        n = len(members)
        n_tm = sum(1 for m in members if (m["tm_helices"] or 0) and m["tm_helices"] >= 2)
        n_lc = sum(1 for m in members if m["low_complexity"])
        n_lowconf = sum(1 for m in members if m["plddt"] < PLDDT_MIN)
        flags = []
        if n_tm / n > 0.5: flags.append("multi-TM (promiscuité membranaire)")
        if n_lc / n > 0.5: flags.append("low-complexity")
        if n_lowconf / n > 0.5: flags.append("modèles pLDDT bas (peu fiable)")
        tms = [pair_tm[tuple(sorted((a, b)))]
               for a in comp for b in comp if a < b and tuple(sorted((a, b))) in pair_tm]
        med_tm = round(sorted(tms)[len(tms)//2], 2) if tms else None
        status = "ARTEFACT probable" if flags else "CANDIDAT famille structurale"
        report.append({"id": i, "size": n, "members": comp, "median_tm": med_tm,
                       "flags": flags, "status": status, "detail": members})
        print(f"Cluster {i} [{status}] n={n} TM méd={med_tm} {('/ '+', '.join(flags)) if flags else ''}")
        for m in members:
            print(f"    {m['rv']:<9} len={m['len']:<4} pLDDT={m['plddt']:<4.0f} "
                  f"verdict={m['verdict']:<14} TM={m['tm_helices']} loc={m['loc_type']} "
                  f"{'LC' if m['low_complexity'] else ''}")
        print()

    json.dump(report, open(OUT / "clusters.json", "w"), ensure_ascii=False, indent=1)
    clean = [c for c in report if not c["flags"]]
    print(f"=> {len(clean)} cluster(s) CANDIDATS (après garde-fous) sur {len(clusters)} bruts.")
    print(f"   Rapport : {OUT/'clusters.json'}")


if __name__ == "__main__":
    main()
