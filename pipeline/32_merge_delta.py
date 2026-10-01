#!/usr/bin/env python3
"""Merge the delta harvest (stage 31) into the atomic sample table.

Design rule: the delta is appended, never reconciled by hand. Every per-sample annotation
on a new row is produced by the same code that produced it for the rest of the atlas --
`classify` (stage 04, imported by stage 31), then the RMS subtype rule (06), the purity
pass (07) and the RNA-assay guard (30), each re-run over the merged table by the caller.
This stage only does what those stages cannot: align the delta's columns to T4's, and
carry across the two annotations that are table-level rather than row-level.

  donor_key      Donor-level grouping, so that two cell lines derived from the same
                 patient (RH4 and RH41) are not counted as two donors. The mapping is
                 curated, so it is read back out of the existing table -- model_matched
                 -> the donor_key already assigned to it -- rather than re-derived. A
                 model the atlas has never seen takes its own name, which is what the
                 existing table does for every unpaired model.

  entity_kind    Defaults to "sarcoma" exactly as stage 05 defaults it; stage 07 demotes
                 the ones that are not.

Not carried across: `is_duplicate` / `duplicate_of_gsm` / `series_redundant`. The
cross-series duplicate pass that produced them is not re-run here, so new rows are
unflagged. The bound on that is small and known: 197 of 85,698 existing rows carry the
flag (0.23%), and the de-duplicated denominators used in every headline figure therefore
move by less than their rounding. This is recorded in docs/methods.md rather than hidden.

Idempotent: GSMs already in T4 are skipped, so re-running adds nothing twice.
"""
import csv, collections, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import resolve, topen, twrite
csv.field_size_limit(10**7)


def main():
    t4_path = resolve("T4_samples_atomic.tsv")
    rows = list(csv.DictReader(topen("T4_samples_atomic.tsv"), delimiter="\t"))
    hdr = list(rows[0].keys())
    have = {r["gsm"] for r in rows}
    print(f"T4 holds {len(rows):,} samples in {len({r['gse'] for r in rows}):,} series")

    # Curated donor grouping, read back out of the table that carries it.
    donor = {}
    for r in rows:
        if r.get("model_matched") and r.get("donor_key"):
            donor.setdefault(r["model_matched"], collections.Counter())[r["donor_key"]] += 1
    donor = {m: c.most_common(1)[0][0] for m, c in donor.items()}
    print(f"donor map covers {len(donor):,} models")

    delta = list(csv.DictReader(topen("T30_delta_samples.tsv.gz"), delimiter="\t"))
    print(f"delta holds {len(delta):,} classified samples")
    new = [r for r in delta if r["in_scope"] == "Y" and r["gsm"] not in have]
    print(f"  in scope and not already held: {len(new):,} "
          f"in {len({r['gse'] for r in new}):,} series")
    if not new:
        print("nothing to merge"); return

    out = []
    for r in new:
        m = r.get("model_matched", "")
        row = {c: r.get(c, "") for c in hdr}
        row["donor_key"] = donor.get(m, m) if m else ""
        row["entity_kind"] = "sarcoma"
        row["is_duplicate"] = row["duplicate_of_gsm"] = row["series_redundant"] = ""
        row["subtype_evidence"] = row["rms_call_basis"] = ""
        out.append(row)

    with twrite(t4_path, gz=True) as f:
        w = csv.DictWriter(f, fieldnames=hdr, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows + out:
            w.writerow(r)
    print(f"\nwrote {len(rows) + len(out):,} samples -> {os.path.basename(t4_path)} "
          f"(+{len(out):,})")

    epi = [r for r in out if r["is_epigenomic"] == "Y"]
    print(f"  of the added rows, {len(epi):,} are epigenomic")
    print("  assay classes:",
          dict(collections.Counter(r["assay_class"] for r in epi).most_common(8)))
    print("  entities:",
          dict(collections.Counter(r["disease"] for r in epi).most_common(8)))
    yrs = collections.Counter(r["gse_date"][:4] for r in out if r["gse_date"])
    print("  release years:", dict(sorted(yrs.items())))


if __name__ == "__main__":
    main()
