#!/usr/bin/env python3
"""Correct samples whose epigenomic assay was inherited from their series, not their own
record.

The assay detector in `04_classify.py` searched the sample's own text first and fell back
to series-level context. For a multi-assay series -- "ChIP-seq and RNA-seq of X", or a 10x
multiome whose ATAC and GEX halves sit side by side -- that fallback stamps the series'
epigenomic assay onto its RNA samples. 658 sarcoma samples carried a DNA-based assay class
while their own GEO `library_strategy` read RNA-Seq, and 623 of them had no epigenomic
assay word anywhere in their own title, source name or characteristics.

`library_strategy` is a structured field set by the depositor and is treated here as
authoritative: a sample GEO calls RNA-Seq is not ChIP-seq, whatever the series is titled.
The 35 that do mention an assay in their own text are almost all the GEX half of a
multiome, where the series legitimately contains ATAC as well.

Stage 04 is fixed for any future full re-run; this stage repairs the shipped table, which
cannot be regenerated without re-running the multi-hour GEO sweep. It is idempotent and
writes `T29_rna_assay_corrections.tsv` so every reassignment is auditable.
"""
import csv, gzip, os, sys, re, collections, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import DATA, resolve, twrite
csv.field_size_limit(10**7)

RNA_STRATEGIES = {"RNA-Seq", "RNA-seq", "miRNA-Seq", "ncRNA-Seq", "RIP-Seq", "ssRNA-seq"}
DNA_ASSAYS = {"ChIP-seq", "ChIP-chip", "ChIP-exo", "CUT&RUN", "CUT&Tag", "ATAC-seq",
              "scATAC-seq", "DNase-seq", "FAIRE-seq", "MNase-seq", "Hi-C", "HiChIP",
              "Micro-C", "Capture-HiC", "ChIA-PET", "4C-seq", "Repli-seq",
              "WGBS", "RRBS", "MeDIP/hMeDIP", "Bisulfite-PCR"}
SINGLE = re.compile(r"single[-\s]cell|single[-\s]nucle|\bscrna|\bsnrna|10x|multiome|"
                    r"smart[-\s]?seq", re.I)
SELF = re.compile(r"chip[-_\s]?seq|cut[-_&\s]?(run|tag)|atac|dnase|faire|mnase|"
                  r"hi[-_]?c|hichip|micro[-_]?c|chia[-_]pet|4c[-_]seq|wgbs|rrbs", re.I)


def main():
    path = resolve("T4_samples_atomic.tsv")
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    cols = list(rows[0].keys())

    log, by_entity = [], collections.Counter()
    for r in rows:
        if r["library_strategy"] not in RNA_STRATEGIES or r["assay_class"] not in DNA_ASSAYS:
            continue
        own = bool(SELF.search(" ".join([r.get("title", ""), r.get("source_name", ""),
                                         r.get("characteristics", "")[:400]])))
        new = "scRNA-seq" if SINGLE.search(" ".join(
            [r.get("title", ""), r.get("characteristics", "")[:300],
             r.get("gse_title", "")])) else "RNA-seq"
        log.append({"gsm": r["gsm"], "gse": r["gse"], "disease": r["disease"],
                    "was_assay": r["assay_class"], "was_target": r["epi_target_norm"],
                    "library_strategy": r["library_strategy"], "now_assay": new,
                    "assay_word_in_own_record": "Y" if own else "",
                    "basis": "series context only" if not own else
                             "own record mentions an assay; library_strategy overrules",
                    "gse_title": r.get("gse_title", "")[:120]})
        if r["entity_kind"] == "sarcoma" and r["is_duplicate"] != "Y":
            by_entity[r["disease"]] += 1
        r["assay_class"] = new
        r["epi_target_norm"] = ""
        r["antibody_raw"] = ""
        r["is_epigenomic"] = "N"

    if not log:
        print("nothing to correct (already applied)")
        return

    tmp = path + ".tmp"
    with (gzip.open if path.endswith(".gz") else open)(tmp, "wt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); [w.writerow(r) for r in rows]
    shutil.move(tmp, path)

    with twrite("T29_rna_assay_corrections.tsv", gz=False) as f:
        w = csv.DictWriter(f, fieldnames=list(log[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); [w.writerow(x) for x in log]

    print(f"reassigned {len(log):,} samples from a DNA-based assay class to RNA")
    print(f"  from series context alone : "
          f"{sum(1 for x in log if not x['assay_word_in_own_record']):,}")
    print(f"  own record mentions assay : "
          f"{sum(1 for x in log if x['assay_word_in_own_record']):,}")
    print(f"  across {len({x['gse'] for x in log}):,} series")
    print("\n  entities most affected (de-duplicated sarcoma samples):")
    for d, n in by_entity.most_common(8):
        print(f"     {d[:36]:38s} {n:5d}")
    print("\n  wrote T29_rna_assay_corrections.tsv")


if __name__ == "__main__":
    main()
