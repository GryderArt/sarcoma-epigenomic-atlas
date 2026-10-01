#!/usr/bin/env python3
"""Purity pass on the v3 atlas.

The bone-tumor bin is populated partly by the H3.3 driver mutations (H3F3A G34W/L in
giant cell tumor of bone, H3F3B K36M in chondroblastoma). The SAME gene carries the
pediatric glioma residues K27M and G34R/V, and a large H3.3 chromatin literature uses
those alleles in cells and in plants. Those records must not sit in a sarcoma bin.
"""
import csv, re, collections
import os
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
csv.field_size_limit(10**7)
O = DATA

GLIOMA = re.compile(r"k27m(?!e)|\bg34r\b|\bg34v\b|\bdipg\b|diffuse midline|glioma|"
                    r"glioblastoma|astrocytom|oligodendroglio|pontine", re.I)
BONE_OK = re.compile(r"giant cell tumou?r|\bgctb\b|chondroblastoma|\bg34w\b|\bg34l\b|"
                     r"\bk36m\b|osteoclast|bone tumou?r", re.I)
NONHUMAN = re.compile(r"arabidopsis|thaliana|drosophila|yeast|cerevisiae|zebrafish|"
                      r"xenopus|c\.? ?elegans", re.I)

def main():
    rows = list(csv.DictReader(topen("T4_samples_atomic.tsv"), delimiter="\t"))
    hdr = list(rows[0].keys())
    moved = collections.Counter(); log = []
    for r in rows:
        if r["disease"] != "GCTB/Chondroblastoma": continue
        b = " ".join((r["title"], r["source_name"], r["characteristics"], r["gse_title"]))
        if BONE_OK.search(b): continue
        why = None
        if NONHUMAN.search(b) or r["organism"] not in ("Homo sapiens", ""):
            why = "non-human H3.3 model system, not a bone tumor"
        elif GLIOMA.search(b):
            why = "H3.3 glioma residue (K27M / G34R / G34V), not the GCTB G34W/L or CB K36M allele"
        if why:
            log.append((r["gsm"], r["gse"], r["disease"], "EXCLUDED", why, r["title"][:52],
                        r["gse_title"][:70]))
            r["disease"] = "EXCLUDED - not a sarcoma"; r["entity_kind"] = "nonsarcoma"
            r["disease_source"] = why
            moved[why.split(",")[0]] += 1
    print("moved out of the bone-tumor bin:", dict(moved), f"({sum(moved.values())} rows)")
    with twrite("T9c_purity_log.tsv") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gsm","gse","from","to","reason","title","gse_title"]); [w.writerow(x) for x in log]
    with twrite("T4_samples_atomic.tsv", gz=True) as f:
        w = csv.DictWriter(f, fieldnames=hdr, delimiter="\t", extrasaction="ignore")
        w.writeheader(); [w.writerow(r) for r in rows]
    g = [r for r in rows if r["disease"] == "GCTB/Chondroblastoma"]
    left = [r for r in g if GLIOMA.search(" ".join((r["title"], r["source_name"],
                                                    r["characteristics"], r["gse_title"])))]
    print(f"GCTB/chondroblastoma bin now {len(g)} rows; residual glioma-wording rows "
          f"(kept because the record also names GCTB/chondroblastoma): {len(left)}")

if __name__ == "__main__":
    main()
