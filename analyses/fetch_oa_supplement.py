#!/usr/bin/env python3
"""fetch_oa_supplement.py — reusable fetcher for open-access supplementary files.

Consolidates the download routes re-derived by hand across the P16 sessions
(springer static-content for Nature/Sci Rep, els-cdn for Elsevier/Cell Press,
PMC bin/ + OA package) into ONE tool, so future work does not re-write ad-hoc
curl chains. Encodes the KB lesson (~/.claude/knowledge/bioinformatics.md,
"Télécharger les supplémentaires d'articles OA ..."): the article HTML often
bounces on an auth interstitial and PMC bin/ serves HTML, but the publisher CDN
serves the static supplementary files in the open — so try the CDN routes and
ALWAYS validate the download with `file` (a ~1 KB "XHTML" .xlsx is an error page).

Routes, by DOI prefix:
  10.1038 (Nature / Sci Rep / Nat Commun) -> scrape static-content.springer.com links
  10.1016 (Elsevier / Cell Press)         -> resolve PII from doi.org, build els-cdn mmc URLs
  any + --pmcid                            -> PMC bin/<file> and PMC OA package tar.gz
  10.1093 (Oxford / NAR)                   -> best-effort article page with a realistic UA
                                              (often anti-bot; may fail -> report honestly)

Usage:
  python fetch_oa_supplement.py --doi 10.1038/s41598-024-82465-w --out DIR
  python fetch_oa_supplement.py --doi 10.1016/j.chembiol.2021.09.002 --out DIR
  python fetch_oa_supplement.py --doi 10.1093/nar/gkag252 --pmcid PMC13096801 --out DIR

Prints, per candidate URL, whether a REAL file was obtained (via `file`), and
exits 0 iff at least one real supplementary file was downloaded. stdlib + curl
+ file only (no extra deps), matching the project environment.
"""
from __future__ import annotations
import argparse
import os
import re
import subprocess
import sys

UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120 Safari/537.36"
# a real supplementary file is a binary office/archive/pdf; an error page is (X)HTML/text
GOOD_FILE = re.compile(r"Microsoft (Excel|Word|OOXML)|Composite Document|Zip archive|"
                       r"PDF document|OpenDocument|CSV text|Excel 2007", re.I)
BAD_FILE = re.compile(r"HTML|XML 1\.0|ASCII text|empty", re.I)


def curl(url: str, out: str | None = None, max_time: int = 60) -> tuple[int, str]:
    """GET url with a realistic UA, following redirects. Returns (http_code, effective_url)."""
    cmd = ["curl", "-sL", "--max-time", str(max_time), "-A", UA]
    if out:
        cmd += ["-o", out, "-w", "%{http_code} %{url_effective}"]
    else:
        cmd += ["-w", "%{http_code} %{url_effective}", "-o", "/dev/null"]
    try:
        r = subprocess.run(cmd + [url], capture_output=True, text=True, timeout=max_time + 15)
        parts = (r.stdout or "").strip().rsplit(" ", 1)
        code = parts[0].split()[-1] if parts and parts[0] else "000"
        eff = parts[1] if len(parts) > 1 else url
        return (int(code) if code.isdigit() else 0), eff
    except Exception:
        return 0, url


def curl_text(url: str, max_time: int = 40) -> str:
    try:
        r = subprocess.run(["curl", "-sL", "--max-time", str(max_time), "-A", UA, url],
                           capture_output=True, text=True, timeout=max_time + 15)
        return r.stdout or ""
    except Exception:
        return ""


def filetype(path: str) -> str:
    try:
        return subprocess.run(["file", "-b", path], capture_output=True, text=True).stdout.strip()
    except Exception:
        return "?"


def is_real(path: str) -> bool:
    if not os.path.exists(path) or os.path.getsize(path) < 1500:
        return False
    ft = filetype(path)
    return bool(GOOD_FILE.search(ft)) and not BAD_FILE.search(ft)


def springer_links(doi: str) -> list[str]:
    html = curl_text(f"https://doi.org/{doi}")
    if len(html) < 2000:  # bounced; try the direct nature.com article URL
        html = curl_text(f"https://www.nature.com/articles/{doi.split('/', 1)[1]}")
    return sorted(set(re.findall(r"https://static-content\.springer\.com[^\"'\s]+", html)))


def elsevier_links(doi: str) -> list[str]:
    html = curl_text(f"https://doi.org/{doi}")
    m = re.search(r"S\d{16}", html)          # PII, e.g. S2451945621004335
    if not m:
        return []
    pii = m.group(0)
    urls = []
    for n in range(1, 9):
        for ext in ("xlsx", "xls", "docx", "pdf", "csv"):
            urls.append(f"https://ars.els-cdn.com/content/image/1-s2.0-{pii}-mmc{n}.{ext}")
    return urls


def pmc_links(pmcid: str) -> list[str]:
    pmcid = pmcid if pmcid.upper().startswith("PMC") else f"PMC{pmcid}"
    urls = []
    # OA web service gives the canonical package href (ftp -> https mirror)
    oa = curl_text(f"https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmcid}")
    for href in re.findall(r'href="(ftp://[^"]+\.tar\.gz)"', oa):
        urls.append(href.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov"))
    # common in-article supplementary filenames under bin/
    for name in ("mmc1.xlsx", "mmc2.xlsx", "media-1.xlsx", "supplementary.xlsx"):
        urls.append(f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/bin/{name}")
    return urls


def oup_links(doi: str) -> list[str]:
    # OUP (NAR) supplementary lives on silverchair; the article page is often anti-bot.
    html = curl_text(f"https://doi.org/{doi}")
    return sorted(set(re.findall(r"https://oup\.silverchair-cdn\.com[^\"'\s]+", html)))


def elife_links(doi: str) -> list[str]:
    # eLife serves supplements on cdn.elifesciences.org, versioned: elife-<id>-supp<N>-v<V>.xlsx.
    # DOI is 10.7554/eLife.<id>; try a few supp indices and file versions.
    m = re.search(r"elife\.(\d+)", doi, re.I)
    if not m:
        return []
    aid = m.group(1)
    urls = []
    for n in range(1, 6):
        for v in (2, 1, 3):
            for ext in ("xlsx", "zip", "docx", "pdf"):
                urls.append(f"https://cdn.elifesciences.org/articles/{aid}/elife-{aid}-supp{n}-v{v}.{ext}")
    return urls


def candidates(doi: str, pmcid: str | None) -> list[str]:
    urls: list[str] = []
    if doi.startswith("10.1038"):
        urls += springer_links(doi)
    elif doi.startswith("10.1016"):
        urls += elsevier_links(doi)
    elif doi.startswith("10.1093"):
        urls += oup_links(doi)
    elif doi.startswith("10.7554"):
        urls += elife_links(doi)
    if pmcid:
        urls += pmc_links(pmcid)
    # de-dup preserving order
    seen, out = set(), []
    for u in urls:
        if u not in seen:
            seen.add(u); out.append(u)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--doi", required=True)
    ap.add_argument("--pmcid", default="")
    ap.add_argument("--out", required=True, help="output directory")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    urls = candidates(args.doi, args.pmcid or None)
    if not urls:
        print(f"[fetch_oa_supplement] no candidate URLs derived for {args.doi} "
              f"(publisher route unknown or article page bounced).")
        return 2
    got = 0
    for i, url in enumerate(urls):
        name = url.rstrip("/").split("/")[-1] or f"suppl_{i}"
        path = os.path.join(args.out, name)
        code, _ = curl(url, out=path)
        if code == 200 and is_real(path):
            print(f"  OK   {name}  ({filetype(path)})  <- {url}")
            got += 1
        else:
            if os.path.exists(path) and not is_real(path):
                os.remove(path)                 # drop error pages
            print(f"  fail [{code}] {name}  <- {url}")
    print(f"[fetch_oa_supplement] {got} real supplementary file(s) in {args.out}")
    if got == 0:
        print("  All routes walled. If the article is very recent or behind anti-bot "
              "(OUP), download the supplementary manually from the article page and "
              "drop it in the output directory.")
    return 0 if got else 1


if __name__ == "__main__":
    sys.exit(main())
