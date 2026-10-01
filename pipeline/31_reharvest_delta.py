#!/usr/bin/env python3
"""Harvest GEO series deposited since the atlas's last cutoff, and classify them.

The full sweep in stages 01-04 queries ~450,000 records and takes hours. Bringing the
atlas up to date does not need that: GEO series carry a release date, so the same term
sweep bounded to dates after the last series already held reaches only what is new.

The delta is classified by importing `classify` from `04_classify`, not by reimplementing
it, so a new sample is annotated by exactly the rules the rest of the atlas was built
with. Every correction applied downstream -- the RMS subtype rule, the purity pass, the
RNA-assay guard -- is re-run over the merged table afterwards by stage 32.

Idempotent: series already present in T4 are skipped, so re-running only picks up what has
appeared since.
"""
import csv, gzip, json, os, sys, re, time, urllib.parse, urllib.request, collections
import importlib.util, concurrent.futures as cf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import DATA, WORK, resolve, twrite
csv.field_size_limit(10**7)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "sarcoma-epigenomic-atlas/1.0 (research; delta harvest)"}
GCOLS = ["gsm", "gse", "families", "title", "source_name", "characteristics", "organism",
         "library_strategy", "library_source", "library_selection", "molecule",
         "instrument", "platform", "sample_type", "suppl_files", "status_date",
         "contact_institute"]


def _mod(name, path):
    spec = importlib.util.spec_from_file_location(name, os.path.join(
        os.path.dirname(os.path.abspath(__file__)), path))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def eget(url, data=None, tries=5):
    for a in range(tries):
        try:
            rq = urllib.request.Request(url, headers=UA,
                                        data=data.encode() if data else None)
            with urllib.request.urlopen(rq, timeout=90) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            time.sleep(1.5 * (a + 1))
    return ""


def esearch(term, lo):
    q = urllib.parse.quote(f"({term}) AND GSE[Entry Type] AND {lo}:3000[PDAT]")
    js = eget(f"{EUTILS}/esearch.fcgi?db=gds&term={q}&retmax=5000&retmode=json")
    try:
        return json.loads(js)["esearchresult"].get("idlist", [])
    except Exception:
        return []


def main():
    # ---- where the atlas currently ends
    t4 = resolve("T4_samples_atomic.tsv")
    have, last = set(), ""
    op = gzip.open if t4.endswith(".gz") else open
    with op(t4, "rt", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("gse"):
                have.add(r["gse"])
            d = r.get("gse_date") or ""
            if d > last:
                last = d
    lo = last.replace("-", "/")[:10] or "2026/01/01"
    print(f"atlas holds {len(have):,} series; latest release date {lo}")

    # ---- the same term sweep as stage 01, plus the curated model names stage 02 uses
    terms = []
    for g, ts in _mod("s01", "01_sweep_disease_terms.py").TERMS.items():
        terms += ts
    models = set()
    for r in csv.DictReader(open(resolve("T2_models.tsv"), newline=""), delimiter="\t"):
        n = (r.get("model_name") or "").strip()
        if len(n) >= 4 and not n.lower().startswith("none"):
            models.add(f'"{n}"')
    terms += sorted(models)
    print(f"sweeping {len(terms):,} queries ({len(models):,} of them model names) "
          f"bounded to {lo} onward")

    uids = set()
    with cf.ThreadPoolExecutor(max_workers=3) as ex:
        for i, got in enumerate(ex.map(lambda t: esearch(t, lo), terms)):
            uids.update(got)
            if (i + 1) % 200 == 0:
                print(f"    {i+1}/{len(terms)} queries, {len(uids):,} uids", flush=True)
    print(f"  {len(uids):,} distinct GEO uids released since {lo}")

    # ---- summaries, then keep only series the atlas does not already hold
    uids = sorted(uids)
    ctx, new = {}, []
    for i in range(0, len(uids), 300):
        js = eget(f"{EUTILS}/esummary.fcgi",
                  data="db=gds&retmode=json&id=" + ",".join(uids[i:i+300]))
        try:
            res = json.loads(js).get("result", {})
        except Exception:
            continue
        for uid in res.get("uids", []):
            d = res[uid]
            acc = d.get("accession", "")
            if not acc.startswith("GSE") or acc in have:
                continue
            ctx[acc] = (f"{d.get('title','')} || {d.get('summary','')}",
                        ";".join(str(x) for x in (d.get("pubmedids") or [])),
                        d.get("PDAT", "") or d.get("pdat", ""),
                        d.get("taxon", ""), "")
            new.append(acc)
    new = sorted(set(new))
    print(f"  {len(new):,} series are new to the atlas")
    if not new:
        print("nothing to add"); return

    # ---- sample records for the new series
    parse_gsm = _mod("s03", "03_harvest_geo.py").parse_gsm
    rows, done = [], 0

    def fetch(gse):
        txt = eget("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                   f"?acc={gse}&targ=gsm&form=text&view=brief")
        return gse, txt

    with cf.ThreadPoolExecutor(max_workers=4) as ex:
        for gse, txt in ex.map(fetch, new):
            done += 1
            if txt:
                rows += parse_gsm(txt, gse, "")
            if done % 50 == 0:
                print(f"    {done}/{len(new)} series, {len(rows):,} samples", flush=True)
    print(f"  {len(rows):,} sample records fetched")

    # ---- classify with the atlas's own rules
    classify = _mod("c04", "04_classify.py").classify
    out, seen = [], set()
    for r in rows:
        g = r.get("gsm")
        if not g or g in seen:
            continue
        seen.add(g)
        full = {c: r.get(c, "") for c in GCOLS}
        out.append(classify(full, ctx))

    keep = [r for r in out if r["in_scope"] == "Y"]
    with twrite("T30_delta_samples.tsv.gz") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()), delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader(); [w.writerow(r) for r in out]

    print(f"\nclassified {len(out):,} samples from {len(new):,} new series")
    print(f"  in scope            : {len(keep):,}")
    print(f"  epigenomic & in scope: "
          f"{sum(1 for r in keep if r['is_epigenomic']=='Y'):,}")
    ac = collections.Counter(r["assay_class"] for r in keep if r["is_epigenomic"] == "Y")
    print("  assay classes:", dict(ac.most_common(10)))
    dc = collections.Counter(r["disease"] for r in keep if r["is_epigenomic"] == "Y")
    print("  entities:", dict(dc.most_common(10)))
    print("\n  wrote T30_delta_samples.tsv.gz")


if __name__ == "__main__":
    main()
