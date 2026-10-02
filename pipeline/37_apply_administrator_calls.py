#!/usr/bin/env python3
"""Let the St Jude data administrator's model list govern GEO samples too.

Correction #11 obtained a definitive per-model diagnosis list from the St Jude data
administrator, and applied it to the CSTN browser tracks it was requested for. It was
never applied to the GEO samples generated on those same models -- and those samples carry
no diagnosis of their own. GSE174376's records say `tissue: xenograft` and nothing else,
so the subtype came from matching the model name against the curated catalogue, which
disagrees with the administrator on most of these O-PDX lines.

The result: 255 samples across six series filed under the wrong RMS subtype, 174 of them
belonging to FN-RMS and sitting under RMS-MYOD1, RMS-NOS or FP-RMS. A reader looking up
FN-RMS saw none of the St Jude PDX work, which is how this was found.

The precedence rule follows correction #3, one level further out. There, curated model
identity outranked series context. Here, the institution that derived, holds and diagnosed
the model outranks our catalogue's reading of it. Matching is on the model base
(SJRHB013758 covers _X1 and _X2), which is safe because no base carries conflicting calls
across its passages -- asserted below rather than assumed.

Also corrects a second fault in the same series. Fifteen samples titled `*_barcode`, with
`protocol: Barcode dialout PCR` and `library_strategy: OTHER`, were classified scATAC-seq
because the series is a single-cell ATAC study. A barcode dial-out is a PCR readout of
lineage labels, not a chromatin profile. This is the series-context leak of corrections #4
and #15 once more, in a series whose own words were clear enough to prevent it: the real
scATAC count for GSE174376 is 7, not 22.

Run after 07, before 08.
"""
import csv, gzip, os, re, sys, shutil, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import resolve, topen, twrite
csv.field_size_limit(10**7)

BASE = re.compile(r"(SJ[A-Z]{2,4}[0-9]+)")
# A barcode dial-out names itself three ways at once; require two so an ordinary sample
# that merely mentions barcodes is never demoted.
BC_TITLE = re.compile(r"_barcode\b", re.I)
BC_PROTO = re.compile(r"barcode\s+dial", re.I)


def admin_map():
    """model base -> entity, from the administrator's list only."""
    by_base = collections.defaultdict(set)
    for r in csv.DictReader(topen("T28_stjude_viz_tracks.tsv"), delimiter="\t"):
        m, e, b = r.get("model_id", ""), r.get("atlas_entity", ""), r.get("entity_call_basis", "")
        if not (m and e) or e in ("normal/reference", "Sarcoma NOS"):
            continue
        if not b.startswith("St Jude data administrator"):
            continue
        g = BASE.match(m)
        if g:
            by_base[g.group(1)].add(e)
    bad = {k: v for k, v in by_base.items() if len(v) > 1}
    if bad:
        sys.exit(f"model bases carry conflicting administrator calls, cannot match on "
                 f"base: {bad}")
    return {k: v.pop() for k, v in by_base.items()}


def main():
    admin = admin_map()
    print(f"administrator list resolves to {len(admin)} model bases")

    path = resolve("T4_samples_atomic.tsv")
    rows = list(csv.DictReader(topen("T4_samples_atomic.tsv"), delimiter="\t"))
    hdr = list(rows[0].keys())

    log, n_sub, n_bc = [], 0, 0
    for r in rows:
        # ---- subtype, from the institution that holds the model
        base = ""
        for f in (r["model_matched"], r["title"], r["source_name"]):
            g = BASE.search(f or "")
            if g:
                base = g.group(1); break
        e = admin.get(base)
        if e and r["disease"] != e:
            log.append([r["gsm"], r["gse"], r["disease"], e, base,
                        "St Jude data administrator per-model list (correction #11)",
                        r["title"][:80]])
            r["disease"] = e
            r["disease_source"] = ("St Jude data administrator per-model list; "
                                   "outranks the curated catalogue for models the "
                                   "institution derived and diagnosed")
            r["subtype_evidence"] = "institutional_diagnosis"
            r["rms_call_basis"] = f"administrator list for {base}"
            n_sub += 1

        # ---- barcode dial-out is not a chromatin profile
        text = " ".join((r["title"], r["source_name"], r["characteristics"]))
        if r["is_epigenomic"] == "Y" and BC_TITLE.search(r["title"] or "") \
                and BC_PROTO.search(text):
            log.append([r["gsm"], r["gse"], r["assay_class"], "not epigenomic", "",
                        "barcode dial-out PCR, not a chromatin assay", r["title"][:80]])
            r["assay_class"] = "barcode-dialout"
            r["epi_target_norm"] = r["antibody_raw"] = ""
            r["is_epigenomic"] = "N"
            n_bc += 1

    tmp = path + ".tmp"
    op = gzip.open if path.endswith(".gz") else open
    with op(tmp, "wt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader(); [w.writerow(r) for r in rows]
    shutil.move(tmp, path)

    with twrite("T32_administrator_corrections.tsv", gz=False) as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["gsm", "gse", "from", "to", "model_base", "basis", "title"])
        [w.writerow(x) for x in log]

    print(f"\nsubtype reassigned from the administrator's list: {n_sub} samples")
    c = collections.Counter((x[2], x[3]) for x in log if x[4])
    for (a, b), n in c.most_common():
        print(f"   {n:4d}  {a:14s} -> {b}")
    print(f"\nbarcode dial-out demoted from epigenomic: {n_bc} samples")
    print("  wrote T32_administrator_corrections.tsv")


if __name__ == "__main__":
    main()
