#!/usr/bin/env python3
"""Apply the species gate to the shipped sample table, and record what each sample is.

A collaborator opened the GEO link for ESS-1, an endometrial stromal sarcoma line, and
found Saccharomyces cerevisiae: GSE147927 is a budding-yeast ChIP-exo atlas, and "Ess1",
the yeast prolyl isomerase, had been read as the cell line. The organism was in the
harvest all along -- nothing had ever asked what it was.

This stage adds the answer to the table and acts on it:

  species        the resolved scientific name(s), from stage 40
  host_class     human | vertebrate | fungus | plant | bacterium | archaeon |
                 invertebrate | other | unknown

and applies `04_classify.species_material()` -- correction #21 -- to every record,
whatever its sample_type. Each change is logged with its basis in T35.

Run after 40, before 08; re-run whenever the table is rebuilt. Idempotent: a second run
changes nothing and keeps the log.
"""
import csv, gzip, os, sys, shutil, collections, importlib.util
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import resolve, topen, twrite
csv.field_size_limit(10**7)


def _mod(name, fn):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(os.path.dirname(os.path.abspath(__file__)), fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


C4 = _mod("c04", "04_classify.py")
LOG = "T35_species_corrections.tsv"
COLS = ["gsm", "gse", "rule", "organism", "species", "host_class", "disease_from",
        "disease_to", "model_from", "model_to", "assay_class", "is_epigenomic",
        "basis", "title", "gse_title"]
EXCLUDED = "EXCLUDED - not a sarcoma"
RMS = {"FP-RMS", "FN-RMS", "RMS-MYOD1"}


def counted(rows):
    """Epigenomic samples the atlas counts for an entity, as the figures count them."""
    c = collections.Counter()
    for r in rows:
        if r.get("entity_kind") == "sarcoma" and r["is_epigenomic"] == "Y" \
                and r["relevance"] != "off_target" and r["disease"]:
            c[r["disease"]] += 1
    return c


def main():
    tax = {r["organism"]: r for r in
           csv.DictReader(topen("T34_organism_taxonomy.tsv"), delimiter="\t")}
    path = resolve("T4_samples_atomic.tsv")
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        hdr, rows = list(rd.fieldnames), list(rd)
    for c in ("species", "host_class"):
        if c not in hdr:
            hdr.insert(hdr.index("organism") + 1, c)
    before = counted(rows)

    log = []
    for r in rows:
        names = [p.strip() for p in r.get("organism", "").replace("/", ";").split(";")
                 if p.strip()]
        r["species"] = ";".join(dict.fromkeys(
            tax[n]["scientific_name"] for n in names if n in tax))
        r["host_class"] = C4.host_class(r.get("organism", ""))
        if r["relevance"] == "off_target" and r["disease"].startswith(EXCLUDED):
            continue          # already out of the counts; re-logging it every run is noise
        frm = (r["disease"], r["model_matched"])
        rules, why = [], []
        for rule, ch, basis in C4.species_material(r):
            if ch.get("clear_model"):
                m = r["model_matched"]
                if r.get("donor_key") in (m, "donor:" + m):
                    r["donor_key"] = ""
                for k in ("model_matched", "model_rrid", "model_subtype",
                          "model_problematic"):
                    r[k] = ""
                if not ch.get("exclude") and r["disease"] in RMS \
                        and m and m in (r.get("rms_call_basis") or ""):
                    r["disease"] = "RMS-NOS"
                    r["subtype_evidence"] = ""
                    r["rms_call_basis"] = (f"reverted to RMS-NOS: the subtype came from "
                                           f"{m}, matched in {r['host_class']} material "
                                           f"(correction #21)")
                if r["relevance"] == "model_confirmed":
                    r["relevance"] = {"sample": "sample_explicit",
                                      "series": "series_explicit"}.get(
                                          r["disease_source"], "sarcoma_nos")
            if ch.get("exclude"):
                r["disease"] = EXCLUDED
                r["entity_kind"] = "nonsarcoma"
                r["disease_source"] = f"correction #21: {basis}"
                r["relevance"] = "off_target"
            rules.append(rule); why.append(basis)
        if not rules:
            continue
        log.append({"gsm": r["gsm"], "gse": r["gse"], "rule": "+".join(rules),
                    "organism": r["organism"], "species": r["species"],
                    "host_class": r["host_class"], "disease_from": frm[0],
                    "disease_to": r["disease"], "model_from": frm[1],
                    "model_to": r["model_matched"], "assay_class": r["assay_class"],
                    "is_epigenomic": r["is_epigenomic"], "basis": " | ".join(why),
                    "title": r["title"][:90], "gse_title": r["gse_title"][:90]})

    # keep earlier entries whose correction still stands, so a re-run does not erase the log
    state = {r["gsm"]: r["disease"] for r in rows}
    done = {x["gsm"] for x in log}
    try:
        for x in csv.DictReader(topen(LOG), delimiter="\t"):
            if x["gsm"] not in done and state.get(x["gsm"]) == x["disease_to"]:
                log.append(x)
    except FileNotFoundError:
        pass
    fresh = len(done)
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

    after = counted(rows)
    by = collections.Counter(x["rule"] for x in log)
    epi = collections.Counter(x["rule"] for x in log if x["is_epigenomic"] == "Y")
    print(f"corrected this run: {fresh} samples; logged in total: {len(log)} "
          f"({sum(epi.values())} epigenomic) -> {LOG}")
    for k in sorted(by):
        print(f"   {k:6s} {by[k]:6d}   epigenomic {epi.get(k, 0):5d}")
    print("\nhost_class across the table:",
          dict(collections.Counter(r["host_class"] for r in rows).most_common()))
    print(f"\nepigenomic samples counted for an entity: {sum(before.values()):,} -> "
          f"{sum(after.values()):,}")
    for d in sorted(set(before) | set(after), key=lambda d: after[d] - before[d]):
        if before[d] != after[d]:
            print(f"   {d:34s} {before[d]:5d} -> {after[d]:5d}")


if __name__ == "__main__":
    main()
