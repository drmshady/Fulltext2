#!/usr/bin/env python3
"""
fetch_fulltext.py - download legally open full texts, one file per study_id.

    pip install requests pypdf
    python fetch_fulltext.py fulltext_retrieval_plan.csv --email you@example.com

Works on fulltext_retrieval_plan.csv or on the original reference_coverage_audit.csv
(needs columns: study_id, title, year, doi, pmid; an optional "fetch" column with
"no" skips that row).

Output
    fulltext/<study_id>.pdf      the paper
    fulltext/<study_id>.xml      Europe PMC full-text XML, when --xml and available
    fulltext_log.csv             one row per study_id: what happened and why

Order of attempts per record
    0. plan  -> oa_url column: first when oa_search_result says "open", otherwise
               tried last (never when it says "no open copy", a trial registration
               or an abstract)
    1. PMCID -> pmcid column, Europe PMC / PMC PDF
       PMID  -> Europe PMC   (authoritative DOI + PMCID; PDF/XML if in PMC)
    2. no id -> Crossref, then OpenAlex title search (accepted only on a close title match)
    3. DOI   -> Unpaywall, OpenAlex, Semantic Scholar OA locations
    4. DOI   -> publisher landing page, <meta name="citation_pdf_url">
               (this is what catches EKB, J-STAGE, SciELO and OJS journals)

Only open-access locations are used. Paywalled papers end up as "no_oa_found"
in the log - that list is your library / ILL request list.
"""
import argparse
import csv
import io
import re
import sys
import time
import unicodedata
from difflib import SequenceMatcher
from pathlib import Path
from urllib.parse import quote, urljoin

import requests

try:
    from pypdf import PdfReader
except ImportError:  # title check is skipped without pypdf
    PdfReader = None

TIMEOUT = 40
BROWSER_UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
              "(KHTML, like Gecko) Chrome/126.0 Safari/537.36")
LOG_FIELDS = ["study_id", "status", "file", "source", "url", "doi_in_csv", "doi_used",
              "doi_mismatch", "pmid", "pmcid", "resolved_by_title", "title_match",
              "pdf_title_check", "note"]


# ----------------------------------------------------------------- helpers ---
def norm(s):
    s = unicodedata.normalize("NFKD", s or "").lower()
    s = re.sub(r"<[^>]+>", " ", s)
    s = re.sub(r"[^a-z0-9\u0600-\u06ff\u3040-\u30ff\u4e00-\u9fff]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def title_similarity(a, b):
    a, b = norm(a), norm(b)
    if not a or not b:
        return 0.0
    return SequenceMatcher(None, a, b).ratio()


def clean_doi(d):
    d = (d or "").strip()
    d = re.sub(r"^https?://(dx\.)?doi\.org/", "", d, flags=re.I)
    return d.lower()


def clean_pmid(p):
    p = (p or "").strip()
    if not p or p.lower() == "nan":
        return ""
    try:
        return str(int(float(p)))
    except ValueError:
        return ""


def clean_title(t):
    t = (t or "").strip()
    if t.lower() == "nan":
        return ""
    if t.startswith("[") and t.endswith("]") and "', '" in t:      # "['arabic', 'english']"
        parts = re.findall(r"'([^']+)'", t)
        latin = [p for p in parts if re.search(r"[A-Za-z]{4,}", p)]
        t = (latin or parts or [t])[-1]
    return t.strip("[]. ")


def clean_pmcid(p):
    m = re.search(r"(?:PMC)?(\d+)", (p or "").strip(), re.I)
    return f"PMC{m.group(1)}" if m else ""


# oa_search_result wording that means the hand-found oa_url is not the paper itself
NOT_THE_PAPER = re.compile(r"^no open copy|trial registration|poster abstract|abstract only", re.I)


def plan_url(row):
    """-> (url, where): where is "first" (marked open), "last" (unconfirmed) or "" (unusable)."""
    url = (row.get("oa_url") or "").strip()
    result = (row.get("oa_search_result") or "").strip()
    if not url.lower().startswith("http") or NOT_THE_PAPER.search(result):
        return url, ""
    if re.match(r"(probable match, )?open\b", result, re.I):
        return url, "first"
    return url, "last"


def pmc_candidates(pmcid):
    return [("europepmc", "https://europepmc.org/backend/ptpmcrender.fcgi"
                          f"?accid={pmcid}&blobtype=pdf"),
            ("europepmc", f"https://europepmc.org/articles/{pmcid}?pdf=render"),
            ("pmc", f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/pdf/")]


class Net:
    def __init__(self, email, delay):
        self.s = requests.Session()
        self.email = email
        self.delay = delay
        self.api_headers = {"User-Agent": f"fulltext-fetch/1.0 (mailto:{email})"}

    def json(self, url, params=None):
        time.sleep(self.delay)
        try:
            r = self.s.get(url, params=params, headers=self.api_headers, timeout=TIMEOUT)
            if r.status_code == 200:
                return r.json()
        except (requests.RequestException, ValueError):
            pass
        return None

    def get(self, url):
        time.sleep(self.delay)
        try:
            return self.s.get(url, headers={"User-Agent": BROWSER_UA, "Accept": "*/*"},
                              timeout=TIMEOUT, allow_redirects=True)
        except requests.RequestException:
            return None


# ---------------------------------------------------------- metadata lookups ---
def europepmc(net, pmid):
    """-> dict(doi, pmcid, is_oa, title, urls[]) or None"""
    j = net.json("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                 {"query": f"EXT_ID:{pmid} AND SRC:MED", "resultType": "core", "format": "json"})
    try:
        hit = j["resultList"]["result"][0]
    except (TypeError, KeyError, IndexError):
        return None
    urls = []
    for u in (hit.get("fullTextUrlList") or {}).get("fullTextUrl", []):
        if u.get("availabilityCode") in ("OA", "F"):          # open access / free
            urls.append((u.get("documentStyle", ""), u.get("url", "")))
    urls.sort(key=lambda x: 0 if x[0] == "pdf" else 1)
    return {"doi": clean_doi(hit.get("doi")), "pmcid": hit.get("pmcid", ""),
            "is_oa": hit.get("isOpenAccess") == "Y", "title": hit.get("title", ""),
            "urls": [u for _, u in urls if u]}


def resolve_by_title(net, title, year, openalex_key):
    """-> (doi, oa_urls, score, source). Accepts only a close title match within +-1 year."""
    best = ("", [], 0.0, "")
    j = net.json("https://api.crossref.org/works",
                 {"query.bibliographic": title, "rows": 5, "mailto": net.email})
    for it in ((j or {}).get("message", {}).get("items", []) or []):
        sc = title_similarity(title, (it.get("title") or [""])[0])
        yr = None
        for k in ("published-print", "published-online", "issued"):
            try:
                yr = it[k]["date-parts"][0][0]; break
            except (KeyError, IndexError, TypeError):
                continue
        if year and yr and abs(int(yr) - int(year)) > 1:
            continue
        if sc > best[2]:
            best = (clean_doi(it.get("DOI")), [], sc, "crossref")
    params = {"search": title, "per-page": 5, "mailto": net.email}
    if openalex_key:
        params["api_key"] = openalex_key
    j = net.json("https://api.openalex.org/works", params)
    for it in ((j or {}).get("results", []) or []):
        sc = title_similarity(title, it.get("title") or it.get("display_name") or "")
        yr = it.get("publication_year")
        if year and yr and abs(int(yr) - int(year)) > 1:
            continue
        if sc > best[2]:
            best = (clean_doi(it.get("doi")), openalex_urls(it), sc, "openalex")
    return best


def openalex_urls(w):
    out = []
    for loc in [w.get("best_oa_location")] + (w.get("locations") or []):
        if loc and loc.get("is_oa"):
            out += [loc.get("pdf_url"), loc.get("landing_page_url")]
    out.append((w.get("open_access") or {}).get("oa_url"))
    return [u for u in out if u]


def oa_candidates(net, doi, openalex_key):
    """All open locations for a DOI, PDF links first. -> list of (source, url)"""
    cands = []
    j = net.json(f"https://api.unpaywall.org/v2/{quote(doi, safe='/')}", {"email": net.email})
    if j:
        locs = [j.get("best_oa_location")] + (j.get("oa_locations") or [])
        for loc in locs:
            if loc:
                cands += [("unpaywall", loc.get("url_for_pdf")), ("unpaywall", loc.get("url"))]
    params = {"mailto": net.email}
    if openalex_key:
        params["api_key"] = openalex_key
    j = net.json(f"https://api.openalex.org/works/https://doi.org/{quote(doi, safe='/')}", params)
    if j:
        cands += [("openalex", u) for u in openalex_urls(j)]
    j = net.json(f"https://api.semanticscholar.org/graph/v1/paper/DOI:{quote(doi, safe='/')}",
                 {"fields": "openAccessPdf"})
    if j and (j.get("openAccessPdf") or {}).get("url"):
        cands.append(("semanticscholar", j["openAccessPdf"]["url"]))
    seen, out = set(), []
    for src, u in cands:
        if u and u not in seen:
            seen.add(u); out.append((src, u))
    out.sort(key=lambda x: 0 if ".pdf" in x[1].lower() or "/pdf" in x[1].lower() else 1)
    return out


META_PDF = re.compile(
    r"""<meta[^>]+(?:name|property)=["']citation_pdf_url["'][^>]*content=["']([^"']+)["']"""
    r"""|<meta[^>]+content=["']([^"']+)["'][^>]*(?:name|property)=["']citation_pdf_url["']""", re.I)


def pdf_from_response(net, resp, depth=0):
    """Return PDF bytes from a response, following citation_pdf_url once if it is HTML."""
    if resp is None or resp.status_code != 200:
        return None, ""
    if b"%PDF-" in resp.content[:1024]:
        return resp.content, resp.url
    if depth >= 2:
        return None, ""
    m = META_PDF.search(resp.text[:400000])
    if m:
        pdf_url = urljoin(resp.url, (m.group(1) or m.group(2)).replace("&amp;", "&"))
        return pdf_from_response(net, net.get(pdf_url), depth + 1)
    return None, ""


def pdf_title_check(pdf_bytes, title):
    """Share of the title's long words found on pages 1-2. '' if it cannot be checked."""
    if PdfReader is None or not re.search(r"[A-Za-z]{4,}", title or ""):
        return ""
    try:
        reader = PdfReader(io.BytesIO(pdf_bytes))
        text = " ".join((p.extract_text() or "") for p in reader.pages[:2])
    except Exception:
        return ""
    text = norm(text)
    if len(text) < 200:
        return "scan_no_text"
    words = [w for w in norm(title).split() if len(w) >= 5]
    if not words:
        return ""
    return f"{sum(w in text for w in words) / len(words):.2f}"


# -------------------------------------------------------------------- main ---
def process(row, net, outdir, args):
    sid = row["study_id"].strip()
    title = clean_title(row.get("title"))
    year = (row.get("year") or "").strip()
    doi_csv = clean_doi(row.get("doi"))
    pmid = clean_pmid(row.get("pmid"))
    log = dict.fromkeys(LOG_FIELDS, "")
    log.update(study_id=sid, doi_in_csv=doi_csv, pmid=pmid)
    target = outdir / f"{sid}.pdf"

    if (row.get("fetch") or "yes").strip().lower() == "no":
        log.update(status="skipped", note="fetch=no in plan")
        return log
    if target.exists() and not args.overwrite:
        log.update(status="already_have", file=target.name)
        return log
    if args.skip_got and (row.get("got") or "").strip().lower() == "yes":
        log.update(status="skipped", note="got=yes in plan (--skip-got)")
        return log

    doi, cands, last = doi_csv, [], []

    # 0. an open copy located by hand in the plan: first if marked open, else as a last resort
    hand_url, where = plan_url(row)
    if where == "first":
        cands.append(("plan_oa_url", hand_url))
    elif where == "last":
        last.append(("plan_oa_url_unconfirmed", hand_url))

    # 1. PMCID from the plan, then PMID -> Europe PMC (authoritative DOI + PMC copy)
    pmcid = clean_pmcid(row.get("pmcid"))
    if pmcid:
        cands += pmc_candidates(pmcid)
    if pmid:
        e = europepmc(net, pmid)
        if e:
            if e["doi"] and doi_csv and e["doi"] != doi_csv:
                log["doi_mismatch"] = f"csv DOI differs from PubMed DOI {e['doi']}"
            if e["doi"]:
                doi = e["doi"]
            if e["pmcid"] and e["pmcid"] != pmcid:
                pmcid = e["pmcid"]
                cands += pmc_candidates(pmcid)
            cands += [("europepmc", u) for u in e["urls"]]
        elif doi_csv:
            log["note"] = "PMID not found in Europe PMC; csv DOI used unverified. "
    log["pmcid"] = pmcid
    if pmcid and args.xml:          # the endpoint answers 404 when the paper is not OA
        r = net.get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML")
        if r is not None and r.status_code == 200 and r.content.lstrip()[:1] == b"<":
            (outdir / f"{sid}.xml").write_bytes(r.content)

    # 2. no DOI yet: resolve by title
    if not doi:
        if not title and not pmid:
            log.update(status="unresolved", note="no DOI, PMID or title")
            return log
        if title:
            rdoi, rurls, score, src = resolve_by_title(net, title, year, args.openalex_key)
            log["title_match"] = f"{score:.2f}"
            if score >= args.min_title_match:
                doi = rdoi
                log["resolved_by_title"] = src
                cands += [(src, u) for u in rurls]
            elif not pmid and not cands and not last:
                log.update(status="unresolved",
                           note=f"no close title match (best {score:.2f}); search manually")
                return log

    # 3 + 4. OA indexes, then the publisher landing page
    if doi:
        log["doi_used"] = doi
        cands += oa_candidates(net, doi, args.openalex_key)
        cands.append(("landing_page", f"https://doi.org/{doi}"))
    cands += last

    tried, seen = 0, set()
    for src, url in cands:
        if url in seen:
            continue
        seen.add(url)
        tried += 1
        pdf, final_url = pdf_from_response(net, net.get(url))
        if pdf:
            target.write_bytes(pdf)
            check = pdf_title_check(pdf, title)
            note = log["note"]
            try:
                if check and check != "scan_no_text" and float(check) < 0.5:
                    note += "PDF first pages do not match the title - open it and check. "
            except ValueError:
                pass
            if src == "plan_oa_url_unconfirmed":
                note += "Came from an oa_url the plan did not confirm as open - check it is the paper. "
            if log["doi_mismatch"]:
                note += "DOI in your CSV points to another paper - fix the record. "
            log.update(status="downloaded", file=target.name, source=src, url=final_url,
                       pdf_title_check=check, note=note.strip())
            return log

    log.update(status="no_oa_found",
               note=(log["note"] + f"{tried} location(s) tried, none returned a PDF").strip())
    return log


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv_path")
    ap.add_argument("--email", required=True, help="real address; Unpaywall/Crossref require one")
    ap.add_argument("--out", default="fulltext")
    ap.add_argument("--log", default="fulltext_log.csv")
    ap.add_argument("--openalex-key", default="", help="OpenAlex API key, if your account needs one")
    ap.add_argument("--xml", action="store_true", help="also save Europe PMC full-text XML")
    ap.add_argument("--only", default="", help="comma-separated study_ids to process")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--skip-got", action="store_true",
                    help='skip rows whose "got" column is "yes" (off by default: those were read '
                         "through a connector, not saved as PDFs)")
    ap.add_argument("--delay", type=float, default=1.0, help="seconds between requests")
    ap.add_argument("--min-title-match", type=float, default=0.88)
    args = ap.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    with open(args.csv_path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    only = {s.strip() for s in args.only.split(",") if s.strip()}
    if only:
        rows = [r for r in rows if r["study_id"].strip() in only]

    net = Net(args.email, args.delay)
    counts = {}
    with open(args.log, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=LOG_FIELDS)
        w.writeheader()
        for i, row in enumerate(rows, 1):
            try:
                log = process(row, net, outdir, args)
            except Exception as exc:            # one bad record must not stop the run
                log = dict.fromkeys(LOG_FIELDS, "")
                log.update(study_id=row.get("study_id", ""), status="error", note=repr(exc))
            w.writerow(log)
            f.flush()
            counts[log["status"]] = counts.get(log["status"], 0) + 1
            print(f"[{i:>3}/{len(rows)}] {log['study_id']:<22} {log['status']:<14} {log['source']}")

    print("\nSummary:", ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    print(f"PDFs in ./{outdir}/, details in {args.log}")
    if PdfReader is None:
        print("pypdf not installed: downloaded PDFs were not checked against their titles.")


if __name__ == "__main__":
    sys.exit(main())
