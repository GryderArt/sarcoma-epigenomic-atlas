#!/usr/bin/env python3
"""Stop the chondroblastoma oncohistone from matching the H3K36 methylation mark.

Correction #22. Chondroblastoma is recognised by the H3.3 K36M substitution, and giant
cell tumour of bone by G34W, so the entity pattern carried `k36m`. Written without a
lookahead it also matches the "K36m" inside H3K36me1, H3K36me2 and H3K36me3 -- the name
of an antibody on a large share of all ChIP-seq and CUT&RUN in GEO. Every such record
read as chondroblastoma: ENCODE histone tracks, neuroblastoma, prostate cancer, and
sarcomas of other entities whose own identity was overwritten. `K27M_GLIOMA` carried the
same flaw against H3K27me1/2/3.

The patterns in `04_classify` now require the mutation form (`k36m(?!e)`, `k27m(?!e)`).
This stage repairs the shipped table, and only where the call rested on that collision
alone: a record keeps its entity if anything else in its text supports it -- "giant cell
tumour of bone", "chondroblastoma", G34W, or a real K36M. Otherwise the disease is
re-derived from the record's own words with the fixed patterns, by the same precedence
`classify` uses, so a sarcoma that was overwritten returns to its own entity rather than
being thrown away.

  22a  the call came only from a methyl mark, and the record names another entity
       -> re-filed under that entity
  22b  ... and it names none                 -> EXCLUDED - not a sarcoma

Run after 41, before 08. Idempotent.
"""
import csv, gzip, json, os, re, shutil, sys, time, collections, importlib.util
import urllib.parse, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import resolve, topen, twrite
csv.field_size_limit(10**7)


def _mod(name, fn):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(os.path.dirname(os.path.abspath(__file__)), fn))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


C4 = _mod("c04", "04_classify.py")
LOG = "T36_oncohistone_corrections.tsv"
COLS = ["gsm", "gse", "rule", "disease_from", "disease_to", "disease_source_to",
        "epi_target", "assay_class", "is_epigenomic", "basis", "title", "gse_title"]
EXCLUDED = "EXCLUDED - not a sarcoma"
# the patterns as they stood before the fix
OLD = {"GCTB/Chondroblastoma":
       re.compile(r"giant\s+cell\s+tumou?r\s+of\s+bone|\bgctb\b|chondroblastoma|g34w|k36m", re.I)}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "sarcoma-epigenomic-atlas/1.0 (research; oncohistone recheck)"}
SUMM = "T36b_series_summaries.tsv"


def summaries(accs):
    """{GSE -> title || summary}. T4 keeps only the series title, so a call made on the
    summary cannot be re-tested from the table; these are fetched once and cached."""
    have = {}
    try:
        for r in csv.DictReader(topen(SUMM), delimiter="\t"):
            have[r["gse"]] = r["context"]
    except FileNotFoundError:
        pass
    todo = sorted(set(accs) - set(have))
    for i in range(0, len(todo), 40):          # the accession query is long; keep the URI short
        chunk = todo[i:i + 40]
        q = urllib.parse.quote(" OR ".join(f"{a}[Accession]" for a in chunk))
        try:
            js = urllib.request.urlopen(urllib.request.Request(
                f"{EUTILS}/esearch.fcgi?db=gds&retmax=500&retmode=json&term={q}",
                headers=UA), timeout=90).read().decode()
            uids = json.loads(js)["esearchresult"].get("idlist", [])
            js = urllib.request.urlopen(urllib.request.Request(
                f"{EUTILS}/esummary.fcgi", headers=UA,
                data=("db=gds&retmode=json&id=" + ",".join(uids)).encode()),
                timeout=90).read().decode()
            res = json.loads(js).get("result", {})
        except Exception as e:
            print(f"  !! summary fetch failed for {len(chunk)} series: {e}")
            continue
        for uid in res.get("uids", []):
            d = res[uid]
            if d.get("accession") in chunk:
                have[d["accession"]] = (f"{d.get('title','')} || {d.get('summary','')}"
                                        ).replace("\t", " ").replace("\n", " ")
        time.sleep(0.4)
    with twrite(SUMM, gz=False) as f:
        w = csv.DictWriter(f, fieldnames=["gse", "context"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for k in sorted(have):
            w.writerow({"gse": k, "context": have[k]})
    return have


def main():
    path = resolve("T4_samples_atomic.tsv")
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        hdr, rows = list(rd.fieldnames), list(rd)

    need = {r["gse"] for r in rows if r["disease"] in OLD
            and r["relevance"] != "off_target" and r["disease_source"] == "series"}
    ctx = summaries(need) if need else {}
    print(f"series summaries available for {len(ctx)} of {len(need)} series-sourced calls")

    log = []
    for r in rows:
        dis = r["disease"]
        if dis not in OLD or r["relevance"] == "off_target":
            continue
        sblob = " ;; ".join([r.get("title", ""), r.get("source_name", ""),
                             r.get("characteristics", ""), r.get("antibody_raw", ""),
                             r.get("epi_target_norm", "")])
        gctx = ctx.get(r["gse"]) or r.get("gse_title", "")
        blob = sblob + " ;; " + gctx
        new = dict(C4.DIS_PATS)[dis]
        if new.search(blob) or not OLD[dis].search(blob):
            continue                      # real evidence, or the call never came from here
        # re-derive by the same precedence classify uses
        dh_s = [k for k, p in C4.DIS_PATS if p.search(sblob)]
        dh_c = [k for k, p in C4.DIS_PATS if p.search(gctx)]
        sp_s = [d for d in dh_s if d != "Sarcoma NOS"]
        sp_c = [d for d in dh_c if d != "Sarcoma NOS"]
        pick = (sp_s or sp_c or dh_s or dh_c or [""])[0]
        src = ("sample" if sp_s else "series" if sp_c else
               "sample" if dh_s else "series" if dh_c else "none")
        frm = dis
        if pick:
            rule, basis = "22a", (f"the chondroblastoma call matched only the K36 methyl "
                                  f"mark; the record's own words name {pick}")
            r["disease"], r["disease_source"] = pick, src
            r["relevance"] = ("sample_explicit" if sp_s else "series_explicit" if sp_c
                              else "sarcoma_nos")
            r["age_class"] = C4.AGE_CLASS.get(pick, r.get("age_class", ""))
        else:
            rule, basis = "22b", ("the chondroblastoma call matched only the K36 methyl "
                                  "mark; nothing else in the record names a sarcoma")
            r["disease"], r["entity_kind"] = EXCLUDED, "nonsarcoma"
            r["disease_source"] = f"correction #22: {basis}"
            r["relevance"] = "off_target"
        log.append({"gsm": r["gsm"], "gse": r["gse"], "rule": rule, "disease_from": frm,
                    "disease_to": r["disease"], "disease_source_to": r["disease_source"][:40],
                    "epi_target": r.get("epi_target_norm", ""),
                    "assay_class": r["assay_class"], "is_epigenomic": r["is_epigenomic"],
                    "basis": basis, "title": r["title"][:90],
                    "gse_title": r["gse_title"][:90]})

    state = {r["gsm"]: r["disease"] for r in rows}
    done = {x["gsm"] for x in log}
    try:
        for x in csv.DictReader(topen(LOG), delimiter="\t"):
            if x["gsm"] not in done and state.get(x["gsm"]) == x["disease_to"]:
                log.append(x)
    except FileNotFoundError:
        pass
    log.sort(key=lambda x: (x["rule"], x["gse"], x["gsm"]))

    tmp = path + ".tmp"
    with op(tmp, "wt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader(); [w.writerow(r) for r in rows]
    shutil.move(tmp, path)
    with twrite(LOG, gz=False) as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader(); [w.writerow(x) for x in log]

    by = collections.Counter(x["rule"] for x in log)
    epi = collections.Counter(x["rule"] for x in log if x["is_epigenomic"] == "Y")
    print(f"corrected {len(done)} samples this run; {len(log)} logged -> {LOG}")
    for k in sorted(by):
        print(f"   {k:5s} {by[k]:5d}   epigenomic {epi.get(k, 0):5d}")
    moved = collections.Counter(x["disease_to"] for x in log if x["rule"] == "22a")
    print("\n  re-filed to:")
    for k, v in moved.most_common():
        print(f"    {v:5}  {k}")


if __name__ == "__main__":
    main()
