#!/usr/bin/env python3
"""Recover epigenomic samples the assay classifier typed as something else.

Correction #23, found while checking whether two Cancer Discovery CDK8 papers
(Zhang et al. and McGinnis et al., Cancer Discov 2026;16(11)) had reached the atlas.
Both had; both were under-counted, by two distinct routes:

  GSE302426  22 H3K27ac ChIP-seq samples read as CRISPR-screen. Their own titles name
             the mark but no assay ("A673_EWSFLI-AID_NChIP-H3K27ac_DMSO_rep1"), so the
             assay word came from the series summary, which describes the paper's pooled
             CRISPR screen. The depositor's own library_strategy said ChIP-Seq.
  GSE303028  50 CUT&RUN samples read as other/unknown. GEO lets CUT&RUN declare
             library_strategy OTHER, so the title has to carry the assay, and these
             abbreviate it "C&R_Rh30_P3F_rep1". The pattern knew "CnR" but not "C&R".

Both rules now live in `04_classify`, and this stage applies them to the shipped table.

  23a  the sample's own words name no assay, the context's answer is not an epigenomic
       assay, and library_strategy names one   -> the strategy's assay
  23b  the sample's own title names CUT&RUN or CUT&Tag and the record is typed as
       something non-epigenomic                 -> that assay

Targets are re-derived for every re-typed row, since a record that was not a ChIP-family
assay never had one read. Run after 43, before 08. Idempotent.
"""
import csv, gzip, os, shutil, sys, collections, importlib.util
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import resolve, topen, twrite
csv.field_size_limit(10**7)


def _mod(name, fn):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(os.path.dirname(os.path.abspath(__file__)), fn))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


C4 = _mod("c04", "04_classify.py")
LOG = "T37_assay_strategy_corrections.tsv"
COLS = ["gsm", "gse", "rule", "library_strategy", "assay_from", "assay_to",
        "target_from", "target_to", "disease", "epi_from", "epi_to", "basis",
        "title", "gse_title"]
CHIPFAM = ("ChIP-seq", "CUT&RUN", "CUT&Tag", "ChIP-exo", "ChIP-chip", "HiChIP", "ChIA-PET")


def retarget(r, assay):
    """The target block of `classify`, for a record that has just become ChIP-family."""
    sblob = " ;; ".join([r.get("title", ""), r.get("source_name", ""),
                         r.get("characteristics", "")])
    if assay in ("WGBS", "RRBS", "Methyl-array", "MeDIP/hMeDIP", "Bisulfite-PCR"):
        return ("5hmC" if "5hmC" in C4.first(C4.MARK_PATS, sblob) else "5mC"), ""
    if assay not in CHIPFAM:
        return "", ""
    m = C4.ENRICH.search(r.get("characteristics", ""))
    if m and C4.INPUTY.search(m.group(1)):
        return "input/none", ""
    if C4.INPUTY.search(r.get("title", "")):
        return "input/none", ""
    a = C4.AB_RE.search(r.get("characteristics", ""))
    ab = (a.group(1).strip() if a else "")[:40]
    combo = ab + " ;; " + sblob
    mk = C4.first(C4.MARK_PATS, combo); tf = C4.first(C4.TF_PATS, combo)
    t = mk[0] if mk else (tf[0] if tf else (ab or "unspecified"))
    return ("input/none" if t in ("input", "IgG") else t), ab


def main():
    path = resolve("T4_samples_atomic.tsv")
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        hdr, rows = list(rd.fieldnames), list(rd)

    log = []
    for r in rows:
        a, ls = r["assay_class"], r.get("library_strategy", "")
        if a in C4.EPI:
            continue
        sblob = " ;; ".join([r.get("title", ""), r.get("source_name", ""),
                             r.get("characteristics", "")])
        own = C4.first(C4.ASSAY_PATS, sblob)
        if own and own[0] in ("CUT&RUN", "CUT&Tag"):
            new, rule = own[0], "23b"
            basis = (f"the sample's own title names {own[0]}; GEO carried it as "
                     f"library_strategy {ls or 'OTHER'}, which names no assay")
        elif not own and C4.LIBSTRAT.get(ls) in C4.EPI:
            new, rule = C4.LIBSTRAT[ls], "23a"
            basis = (f"the sample's own words name no assay and '{a}' came from the "
                     f"series; the depositor's library_strategy is {ls}")
        else:
            continue
        frm = (a, r.get("epi_target_norm", ""))
        tgt, ab = retarget(r, new)
        r["assay_class"] = new
        r["epi_target_norm"] = tgt
        if ab and not r.get("antibody_raw"):
            r["antibody_raw"] = ab
        r["is_epigenomic"] = "Y" if new in C4.EPI else "N"
        log.append({"gsm": r["gsm"], "gse": r["gse"], "rule": rule,
                    "library_strategy": ls, "assay_from": frm[0], "assay_to": new,
                    "target_from": frm[1], "target_to": tgt, "disease": r["disease"],
                    "epi_from": "N", "epi_to": r["is_epigenomic"], "basis": basis,
                    "title": r["title"][:90], "gse_title": r["gse_title"][:90]})

    state = {r["gsm"]: r["assay_class"] for r in rows}
    done = {x["gsm"] for x in log}
    try:
        for x in csv.DictReader(topen(LOG), delimiter="\t"):
            if x["gsm"] not in done and state.get(x["gsm"]) == x["assay_to"]:
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
    print(f"recovered this run: {len(done)} samples; {len(log)} logged -> {LOG}")
    for k in sorted(by):
        print(f"   {k:5s} {by[k]:5d}")
    print("\n  assays recovered:",
          dict(collections.Counter(x["assay_to"] for x in log).most_common()))
    print("  entities:",
          dict(collections.Counter(x["disease"] for x in log).most_common(8)))


if __name__ == "__main__":
    main()
