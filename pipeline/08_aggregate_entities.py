#!/usr/bin/env python3
"""Re-derive the per-disease model / patient / assay counts from the v3 atomic table
(v3 = v2 + the GSE140686 entity recovery). Writes T3_entity_counts.tsv."""
import csv, collections, os
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
csv.field_size_limit(10**7)
O = DATA

EPI = {"ChIP-seq","CUT&RUN","CUT&Tag","ChIP-exo","ChIP-chip","ATAC-seq","scATAC-seq",
       "DNase-seq","FAIRE-seq","MNase-seq","Hi-C","HiChIP","Micro-C","Capture-HiC",
       "ChIA-PET","4C-seq","Repli-seq","WGBS","RRBS","Methyl-array","MeDIP/hMeDIP",
       "Bisulfite-PCR"}
ACC = {"ATAC-seq","scATAC-seq","DNase-seq","FAIRE-seq"}
D3  = {"Hi-C","HiChIP","Micro-C","Capture-HiC","ChIA-PET","4C-seq"}
DNAME = {"WGBS","RRBS","Methyl-array","MeDIP/hMeDIP","Bisulfite-PCR"}
MODEL_ST  = {"cell_line","organoid","xenograft(CDX)","mouse_model"}
PDX_ST    = {"PDX"}
PATIENT_ST = {"primary_tumor","metastasis","recurrence"}

def main():
    rows = list(csv.DictReader(topen("T4_samples_atomic.tsv"), delimiter="\t"))
    epi  = [r for r in rows if r["assay_class"] in EPI]
    ded  = [r for r in epi if r["is_duplicate"] != "Y"]

    acls = {}
    for r in rows:
        if r["disease"] and r["age_class"]:
            acls.setdefault(r["disease"], collections.Counter())[r["age_class"]] += 1

    dis = sorted({r["disease"] for r in rows if r["disease"]})
    out = []
    for d in dis:
        dr  = [r for r in ded if r["disease"] == d]
        dall= [r for r in epi if r["disease"] == d]
        drow= [r for r in rows if r["disease"] == d]
        def units(pred):
            return {r["donor_key"] or r["model_matched"] for r in drow
                    if pred(r) and (r["donor_key"] or r["model_matched"])}
        mu  = units(lambda r: r["sample_type"] in MODEL_ST)
        px  = units(lambda r: r["sample_type"] in PDX_ST)
        pat = {(r["gse"], r["gsm"]) for r in drow if r["sample_type"] in PATIENT_ST
               and r["assay_class"] in EPI and r["is_duplicate"] != "Y"}
        def cnt(pred): return sum(1 for r in dr if pred(r))
        ismodel = lambda r: r["sample_type"] in MODEL_ST
        ispdx   = lambda r: r["sample_type"] in PDX_ST
        ispat   = lambda r: r["sample_type"] in PATIENT_ST
        k27 = lambda r: r["epi_target_norm"] == "H3K27ac"
        out.append({
            "disease": d,
            "entity_kind": collections.Counter(r["entity_kind"] for r in drow).most_common(1)[0][0]
                            if drow else "sarcoma",
            "age_class": acls.get(d, collections.Counter()).most_common(1)[0][0] if d in acls else "",
            "model_units": len(mu), "pdx_units": len(px),
            "patient_epigenomes": len(pat),
            "model_H3K27ac": cnt(lambda r: ismodel(r) and k27(r)),
            "patient_H3K27ac": cnt(lambda r: ispat(r) and k27(r)),
            "pdx_H3K27ac": cnt(lambda r: ispdx(r) and k27(r)),
            "model_accessibility": cnt(lambda r: ismodel(r) and r["assay_class"] in ACC),
            "patient_accessibility": cnt(lambda r: ispat(r) and r["assay_class"] in ACC),
            "model_3D": cnt(lambda r: ismodel(r) and r["assay_class"] in D3),
            "patient_3D": cnt(lambda r: ispat(r) and r["assay_class"] in D3),
            "model_DNAme": cnt(lambda r: ismodel(r) and r["assay_class"] in DNAME),
            "patient_DNAme": cnt(lambda r: ispat(r) and r["assay_class"] in DNAME),
            "epigenomic_samples_raw": len(dall),
            "epigenomic_samples_dedup": len(dr),
            "epi_ChIPseq": cnt(lambda r: r["assay_class"] in ("ChIP-seq","CUT&RUN","CUT&Tag",
                                                             "ChIP-exo","ChIP-chip")),
            "epi_accessibility": cnt(lambda r: r["assay_class"] in ACC),
            "epi_3D": cnt(lambda r: r["assay_class"] in D3),
            "epi_DNAme": cnt(lambda r: r["assay_class"] in DNAME),
            "total_samples_any_assay": len(drow),
        })
    out.sort(key=lambda r: -r["epigenomic_samples_dedup"])
    with twrite("T3_entity_counts.tsv") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()), delimiter="\t")
        w.writeheader(); [w.writerow(r) for r in out]
    print(f"{len(out)} disease bins -> T3_entity_counts.tsv")
    print(f"{len(ded):,} deduplicated epigenomic samples ({len(epi):,} raw)")
    print(f"\n{'disease':34s} {'kind':10s} {'epi':>6s} {'K27ac':>6s} {'ATAC':>5s} {'3D':>4s} "
          f"{'DNAme':>6s} {'models':>7s} {'PDX':>4s} {'pt-epi':>7s}")
    for r in out:
        print(f"{r['disease'][:33]:34s} {r['entity_kind'][:9]:10s} {r['epigenomic_samples_dedup']:6d} "
              f"{r['model_H3K27ac']+r['patient_H3K27ac']+r['pdx_H3K27ac']:6d} {r['epi_accessibility']:5d} "
              f"{r['epi_3D']:4d} {r['epi_DNAme']:6d} {r['model_units']:7d} {r['pdx_units']:4d} "
              f"{r['patient_epigenomes']:7d}")

if __name__ == "__main__":
    main()
