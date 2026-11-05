#!/usr/bin/env python3
"""What GEO has added since the atlas's last scrape, and what landed in the last month.

Reads the delta table stage 31 wrote and reports it two ways: everything new since the
cutoff, and the trailing month on its own, which is the number worth watching because it
says whether the field is still depositing into these entities.

Counts only what the atlas counts: a sample in scope, flagged epigenomic, on material that
can carry a sarcoma.
"""
import csv, collections, datetime, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import topen

WINDOW = 30


def main():
    # The delta is reported from the merged table, not from T30: every correction stage
    # runs over T4 after the merge, so T30 still holds the pre-correction call.
    delta = {r["gse"] for r in csv.DictReader(topen("T30_delta_samples.tsv"), delimiter="\t")}
    rows = [r for r in csv.DictReader(topen("T4_samples_atomic.tsv"), delimiter="\t")
            if r["gse"] in delta]
    for r in rows:
        r.setdefault("in_scope", "N" if r["relevance"] == "off_target" else "Y")
        r["in_scope"] = "N" if r["relevance"] == "off_target" else "Y"
    if not rows:
        print("no delta series in T4 -- run 31 then 32"); return
    dates = sorted({r["gse_date"] for r in rows if r.get("gse_date")})
    newest = datetime.date.fromisoformat(dates[-1].replace("/", "-")[:10])
    cut = newest - datetime.timedelta(days=WINDOW)
    print(f"delta covers {dates[0]} .. {dates[-1]}   "
          f"({len(rows):,} samples in {len({r['gse'] for r in rows}):,} series)\n")

    def window(r):
        d = (r.get("gse_date") or "").replace("/", "-")[:10]
        try:
            return datetime.date.fromisoformat(d) > cut
        except ValueError:
            return False

    for label, sel in (("since the last scrape", lambda r: True),
                       (f"last {WINDOW} days (after {cut})", window)):
        sub = [r for r in rows if sel(r)]
        ins = [r for r in sub if r.get("in_scope") == "Y"]
        epi = [r for r in ins if r["is_epigenomic"] == "Y"
               and not r["disease"].startswith("EXCLUDED")]
        print(f"=== {label} ===")
        print(f"  series seen                 {len({r['gse'] for r in sub}):>7,}")
        print(f"  samples seen                {len(sub):>7,}")
        print(f"  in scope                    {len(ins):>7,}")
        print(f"  epigenomic and in scope     {len(epi):>7,}"
              f"   in {len({r['gse'] for r in epi}):,} series")
        if epi:
            print("\n  assay classes")
            for k, v in collections.Counter(r["assay_class"] for r in epi).most_common():
                print(f"    {v:5}  {k}")
            print("\n  entities")
            for k, v in collections.Counter(r["disease"] for r in epi).most_common():
                print(f"    {v:5}  {k}")
            print("\n  sample types")
            for k, v in collections.Counter(r["sample_type"] for r in epi).most_common():
                print(f"    {v:5}  {k}")
            print("\n  series, largest first")
            by = collections.defaultdict(list)
            for r in epi:
                by[(r["gse"], r["gse_title"][:62], r["gse_date"])].append(r)
            for (g, t, d), rs in sorted(by.items(), key=lambda kv: -len(kv[1]))[:18]:
                ent = collections.Counter(x["disease"] for x in rs).most_common(1)[0][0]
                print(f"    {len(rs):4}  {g:<11} {d}  {ent[:22]:<22} {t}")
        print()

    blocked = [r for r in rows if r["disease"].startswith("EXCLUDED")
               and "correction #2" in (r.get("disease_source") or "")]
    print(f"the gates turned away {len(blocked):,} of the new samples "
          f"({len({r['gse'] for r in blocked}):,} series) before they could be counted")
    for k, v in collections.Counter(r["host_class"] for r in blocked).most_common():
        print(f"    {v:5}  {k}")


if __name__ == "__main__":
    main()
