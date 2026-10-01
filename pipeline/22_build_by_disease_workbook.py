#!/usr/bin/env python3
"""SarcomaEpigenomics_by_disease_v3.xlsx -- one tab per sarcoma entity, at the resolution
the working group asked for: 'H3K27ac ChIP-seq for RH4', one row at a time.

Columns follow the agreed layout: model type | model name | datatype | GEO source,
with the target, sample counts, links and quality flags to the right of them."""
import csv, collections, os, re
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
csv.field_size_limit(10**7)
O = DATA
OUT = os.path.join(WORKBOOKS, "SarcomaEpigenomics_by_disease.xlsx")

ARIAL = "Arial"
HDR_FIL = PatternFill("solid", fgColor="0D366B")
HDR_FNT = Font(name=ARIAL, size=9, bold=True, color="FFFFFF")
BODY = Font(name=ARIAL, size=9)
LINK = Font(name=ARIAL, size=9, color="0563C1", underline="single")
WARN = PatternFill("solid", fgColor="FDE9E4")
TITLE = Font(name=ARIAL, size=11, bold=True, color="0D366B")
SUB = Font(name=ARIAL, size=9, italic=True, color="4A5560")
BOLD = Font(name=ARIAL, size=9, bold=True)

EPI = {"ChIP-seq","CUT&RUN","CUT&Tag","ChIP-exo","ChIP-chip","ATAC-seq","scATAC-seq",
       "DNase-seq","FAIRE-seq","MNase-seq","Hi-C","HiChIP","Micro-C","Capture-HiC",
       "ChIA-PET","4C-seq","Repli-seq","WGBS","RRBS","Methyl-array","MeDIP/hMeDIP",
       "Bisulfite-PCR"}
REG = {"ChIP-seq","CUT&RUN","CUT&Tag","ChIP-exo","ChIP-chip","ATAC-seq","scATAC-seq",
       "DNase-seq","FAIRE-seq","MNase-seq","Hi-C","HiChIP","Micro-C","Capture-HiC",
       "ChIA-PET","4C-seq"}
DNAME = {"WGBS","RRBS","Methyl-array","MeDIP/hMeDIP","Bisulfite-PCR"}
PATIENT = {"primary_tumor","metastasis","recurrence"}
MODELISH = {"cell_line","organoid","xenograft(CDX)","mouse_model","PDX"}

MODEL_TYPE = {"cell_line":"cell line","organoid":"organoid","xenograft(CDX)":"CDX",
              "mouse_model":"GEM / mouse model","PDX":"PDX","primary_tumor":"primary tumor",
              "metastasis":"metastasis","recurrence":"recurrence",
              "normal/reference":"normal / reference","unspecified":"unspecified"}

COLS = ["model type","model name","datatype","GEO source","epigenetic target",
        "n samples","sample type","organism","assay class","subtype evidence",
        "series title","PubMed","year","GEO link","PubMed link","model flag",
        "duplicate note"]

def tabname(d):
    s = re.sub(r"[\[\]\:\*\?/\\]", "-", d)[:31]
    return s or "unnamed"

def main():
    rows = list(csv.DictReader(topen("T4_samples_atomic.tsv"), delimiter="\t"))
    epi = [r for r in rows if r["assay_class"] in EPI and r["is_duplicate"] != "Y"]
    prob = {}
    for m in csv.DictReader(topen("T2_models.tsv"), delimiter="\t"):
        if (m.get("problematic_flag") or "").strip() and m["problematic_flag"].lower() not in ("no","none",""):
            prob[m["model_name"].strip().lower()] = m["problematic_flag"][:120]

    k686 = {r["gsm"]: r for r in csv.DictReader(topen("GSE140686_recovered_diagnoses.tsv"),
                                                delimiter="\t")}
    by = collections.defaultdict(list)
    for r in epi:
        if r["entity_kind"] != "sarcoma": continue
        by[r["disease"]].append(r)
    order = sorted(by, key=lambda d: -len(by[d]))

    wb = openpyxl.Workbook(); wb.remove(wb.active)

    # ---- index tab
    ws = wb.create_sheet("INDEX")
    ws["A1"] = "Sarcoma epigenomics by disease — atlas v3 (all ages)"; ws["A1"].font = TITLE
    ws["A2"] = ("One tab per entity. Each row is one unique combination of model, assay and "
                "epigenetic target within one GEO series — the 'H3K27ac ChIP-seq for RH4' "
                "resolution. Duplicated samples and ChIP input controls are labeled, not "
                "silently dropped.")
    ws["A2"].font = SUB; ws["A2"].alignment = Alignment(wrap_text=True, vertical="top")
    ws.row_dimensions[2].height = 42
    hdr = ["tab","entity","epigenomic samples","regulatory","DNA methylation",
           "H3K27ac","accessibility","3D genome","primary-tumor regulatory",
           "distinct models","distinct series"]
    for j, h in enumerate(hdr, 1):
        c = ws.cell(row=4, column=j, value=h); c.font = HDR_FNT; c.fill = HDR_FIL
    r0 = 5
    for d in order:
        g = by[d]
        isreg = lambda x: x["assay_class"] in REG and x["epi_target_norm"] not in ("input/none","none")
        vals = [tabname(d), d, len(g), sum(1 for x in g if isreg(x)),
                sum(1 for x in g if x["assay_class"] in DNAME),
                sum(1 for x in g if x["epi_target_norm"] == "H3K27ac"),
                sum(1 for x in g if x["assay_class"] in ("ATAC-seq","scATAC-seq","DNase-seq","FAIRE-seq")),
                sum(1 for x in g if x["assay_class"] in ("Hi-C","HiChIP","Micro-C","Capture-HiC","ChIA-PET","4C-seq")),
                sum(1 for x in g if x["sample_type"] in PATIENT and isreg(x)),
                len({x["model_matched"] for x in g if x["model_matched"]}),
                len({x["gse"] for x in g})]
        for j, v in enumerate(vals, 1):
            c = ws.cell(row=r0, column=j, value=v); c.font = BODY
            if j in (4, 9) and v == 0: c.fill = WARN; c.font = Font(name=ARIAL, size=9, color="C0392B", bold=True)
        r0 += 1
    ws.freeze_panes = "A5"
    ws.auto_filter.ref = f"A4:{get_column_letter(len(hdr))}{r0-1}"
    for j, w in enumerate([26,34,19,11,17,10,14,11,25,15,14], 1):
        ws.column_dimensions[get_column_letter(j)].width = w

    # ---- one tab per entity
    for d in order:
        g = by[d]
        agg = collections.OrderedDict()
        for r in g:
            if r["gse"] == "GSE140686":
                k = k686.get(r["gsm"], {})
                name = ("DKFZ classifier " + (k.get("set") or "") + " cases — "
                        + (k.get("institutional_diagnosis") or r["disease"]))[:70]
            elif r["model_matched"]:
                name = r["model_matched"]
            elif r["sample_type"] in PATIENT:
                name = "patient tumor: " + (r["source_name"][:40] or r["title"][:40]
                                             or "unlabeled")
            else:
                name = r["source_name"][:44] or r["title"][:44] or "unnamed"
            key = (r["sample_type"], name, r["assay_class"], r["epi_target_norm"], r["gse"])
            a = agg.setdefault(key, {"n": 0, "row": r})
            a["n"] += 1
        ws = wb.create_sheet(tabname(d))
        for j, h in enumerate(COLS, 1):
            c = ws.cell(row=1, column=j, value=h); c.font = HDR_FNT; c.fill = HDR_FIL
        i = 2
        for (st, name, assay, target, gse), a in sorted(
                agg.items(), key=lambda kv: (kv[0][0] not in MODELISH, kv[0][0], kv[0][1],
                                             kv[0][2], kv[0][3])):
            r = a["row"]
            pm = r["gse_pubmed"]
            flag = prob.get((r["model_matched"] or "").strip().lower(), "")
            dup = "input / IgG control" if target in ("input/none","none") else ""
            vals = [MODEL_TYPE.get(st, st), name, assay, gse,
                    target or "—", a["n"], st, r["organism"], assay,
                    r.get("subtype_evidence") or r.get("rms_call_basis") or "",
                    r["gse_title"][:110], pm, (r["gse_date"] or "")[:4],
                    f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}",
                    f"https://pubmed.ncbi.nlm.nih.gov/{pm}/" if pm else "",
                    flag, dup]
            for j, v in enumerate(vals, 1):
                c = ws.cell(row=i, column=j, value=v)
                c.font = LINK if (j in (14, 15) and v) else BODY
                if j in (14, 15) and v: c.hyperlink = v
                if j == 16 and v: c.fill = WARN
            i += 1
        ws.freeze_panes = "A2"
        if i > 2: ws.auto_filter.ref = f"A1:{get_column_letter(len(COLS))}{i-1}"
        for j, w in enumerate([17,30,15,13,26,10,17,15,15,34,52,11,7,17,17,34,20], 1):
            ws.column_dimensions[get_column_letter(j)].width = w

    wb.save(OUT)
    print(f"{len(order)} entity tabs + INDEX -> {OUT} ({os.path.getsize(OUT)/1e6:.1f} MB)")
    for d in order[:12]:
        print(f"   {tabname(d):32s} {len(by[d]):6d} epigenomic samples")

if __name__ == "__main__":
    main()
