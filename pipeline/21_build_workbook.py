#!/usr/bin/env python3
"""SarcomaEpigenomicDataAtlas_v3.xlsx -- all-ages sarcoma epigenomic data atlas.
Built for clean Google Sheets import: one flat header row per tab, no merged cells,
no formulas that Sheets would have to re-resolve, frozen header, autofilter."""
import csv, os, collections
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
csv.field_size_limit(10**7)
O = DATA
OUT = os.path.join(WORKBOOKS, "SarcomaEpigenomicDataAtlas.xlsx")

ARIAL   = "Arial"
HDR_FIL = PatternFill("solid", fgColor="0D366B")
HDR_FNT = Font(name=ARIAL, size=9, bold=True, color="FFFFFF")
BODY    = Font(name=ARIAL, size=9)
TITLE   = Font(name=ARIAL, size=11, bold=True, color="0D366B")
NOTE    = Font(name=ARIAL, size=9, italic=True, color="4A5560")

SHEETS = [
    ("README",            None),
    ("T13_Incidence_vs_Data", "T13_incidence_vs_data.tsv"),
    ("G3_EntityCounts",   "T3_entity_counts.tsv"),
    ("T2_Models",         "T2_models.tsv"),
    ("T1_Taxonomy",       "T1_taxonomy.tsv"),
    ("T4_Samples_atomic", "T4_samples_atomic.tsv"),
    ("T14_EGA_sarcoma",   "T14_ega_sarcoma_datasets.tsv"),
    ("T15_Controlled_epigenomics", "T15_controlled_access.tsv"),
    ("T16_StJude_CSTN",   "T16_stjude_cstn_inventory.tsv"),
    ("T17_StJude_OPDX_models", "T17_stjude_opdx_models.tsv"),
    ("T12c_EBI_unique",   "T12c_ebi_unique.tsv"),
    ("T12b_EBI_sarcoma_epi", "T12b_ebi_sarcoma_epigenomic.tsv"),
    ("T12_EBI_all",       "T12_ebi_all.tsv"),
    ("GSE140686_recovered", "GSE140686_recovered_diagnoses.tsv"),
    ("T9b_RMS_recall_log", "T9b_rms_subtype_log.tsv"),
    ("T9c_Purity_log",    "T9c_purity_log.tsv"),
    ("T11_H3K27ac_fidelity", "T11_h3k27ac_fidelity.tsv"),
    ("T11b_fidelity_by_study", "T11b_fidelity_by_study.tsv"),
]

README = [
 ("Sarcoma Epigenomic Data Atlas — v3 (all ages)", "title"),
 ("", ""),
 ("Built for the SASS Epigenetics Working Group white paper, section: "
  "Critical assessment of existing datasets.", ""),
 ("Scope: paediatric AND adult sarcoma. Sources: NCBI GEO (sample-level), "
  "EBI ArrayExpress / BioStudies, ENA, EGA, EpiRR (IHEC).", ""),
 ("", ""),
 ("WHAT CHANGED IN v3", "head"),
 ("1. GSE140686 recovered. The DKFZ sarcoma methylation classifier (Koelsche et al., "
  "Nat Commun 2021, PMID 33479225) deposited 1,505 arrays in GEO labelled only "
  "'sarcoma classifier reference case N'; every sample reads 'tissue: sarcoma'. The "
  "per-case diagnosis exists only in the paper's Supplementary Data, joined through an "
  "id GEO hides in !Sample_description. Recovered here: 1,315 sarcoma, 154 benign or "
  "control, 36 non-sarcoma, across 46 entities.", ""),
 ("2. RMS subtypes re-derived on the user's rule: fusion-negative = neither a PAX3/7 "
  "fusion nor a mutant MYOD1 (a RAS mutation is NOT the criterion). Curated cell line "
  "identity now outranks series context, so PAX-fusion lines (RH4/RH30/RH41/CW9019) "
  "can no longer be swept into FN-RMS by their TP53 or other annotations.", ""),
 ("3. Bone-tumour bin purified. H3F3A/H3F3B capture pulled in H3.3 K27M and G34R/V "
  "records (paediatric glioma residues, plus plant and mouse model systems); 288 rows "
  "removed. GCTB uses G34W/L, chondroblastoma uses K36M.", ""),
 ("4. ChIP/CUT&RUN input and IgG controls are excluded from 'regulatory epigenomic "
  "samples' — they are controls, not profiles.", ""),
 ("", ""),
 ("KEY DEFINITIONS", "head"),
 ("regulatory epigenomics = ChIP-seq, ChIP-exo, ChIP-chip, CUT&RUN, CUT&Tag, ATAC-seq, "
  "scATAC-seq, DNase-seq, FAIRE-seq, MNase-seq, Hi-C, HiChIP, Micro-C, Capture-HiC, "
  "ChIA-PET, 4C-seq — excluding input/IgG controls.", ""),
 ("DNA methylation = WGBS, RRBS, methylation array, MeDIP/hMeDIP, bisulfite PCR.", ""),
 ("epigenomic = regulatory + DNA methylation + input controls.", ""),
 ("age_class = which age group the entity predominantly affects (paediatric / both / "
  "adult); it is a property of the DISEASE, not of the individual sample.", ""),
 ("", ""),
 ("HOW TO READ T13", "head"),
 ("US_cases_per_year_all_ages is a curated anchor with its source in incidence_basis. "
  "US_low / US_high bracket the disagreement between US registry coding and European "
  "central-pathology-review series (NETSARC), which for several entities differ 2-fold. "
  "Entities flagged no_population_rate=Y have NO published population rate anywhere; "
  "they are shown but excluded from every per-case ratio.", ""),
 ("", ""),
 ("NEW IN THIS BUILD: CONTROLLED-ACCESS AND ST JUDE", "head"),
 ("T14 is the full EGA sarcoma catalogue: 560 datasets, 22,832 samples, 100% controlled "
  "access. T15 maps the epigenomic subset onto atlas entities. T16/T17 are the St Jude "
  "Childhood Solid Tumor Network inventory and its 389 O-PDX / cell models.", ""),
 ("The distinction these tabs add: a GEO-only survey cannot tell 'never measured' from "
  "'measured but locked'. 1,207 regulatory epigenomes and 7,591 methylation samples exist "
  "for sarcoma behind data-access committees. The entire St Jude CSTN epigenome — 756 "
  "ChIP-seq libraries plus WGBS — is EGA-only; none of it is in GEO or in St Jude Cloud.", ""),
 ("Crucially, none of the controlled-access holdings cover the nine entities that have no "
  "regulatory epigenomics at all. Those zeros survive the addition of EGA.", ""),
 ("", ""),
 ("CAVEATS THAT MATTER", "head"),
 ("· A sample count is not an information count. 4,511 osteosarcoma epigenomes are not "
  "4,511 independent observations; see G3 for distinct models and distinct patients.", ""),
 ("· Absence in this atlas means absence of a public, entity-labelled deposit — not "
  "proof that no experiment was ever done. Data in EGA under controlled access, or "
  "deposited without an entity label, is invisible to any search of this kind. "
  "GSE140686 is the proof that this failure mode is real and large.", ""),
 ("· RMS subtype calls carry their evidence in rms_call_basis and, for GSE140686, in "
  "subtype_evidence. Calls from histology or methylation class are weaker than calls "
  "from a stated fusion.", ""),
]

def main():
    wb = openpyxl.Workbook(); wb.remove(wb.active)
    for name, fn in SHEETS:
        ws = wb.create_sheet(name[:31])
        if fn is None:
            r = 1
            for text, kind in README:
                c = ws.cell(row=r, column=1, value=text)
                c.font = TITLE if kind == "title" else (
                         Font(name=ARIAL, size=10, bold=True) if kind == "head" else BODY)
                c.alignment = Alignment(wrap_text=True, vertical="top")
                ws.row_dimensions[r].height = None if len(text) < 90 else 14*(len(text)//95 + 1)
                r += 1
            ws.column_dimensions["A"].width = 132
            continue
        try:
            fh = topen(fn)
        except FileNotFoundError:
            ws.cell(row=1, column=1,
                    value=f"(source table {fn} not present - run the pipeline stage "
                          f"that produces it)").font = NOTE
            print(f"  !! MISSING {fn} - tab {name} left empty")
            continue
        with fh as f:
            rd = csv.reader(f, delimiter="\t")
            for i, row in enumerate(rd, start=1):
                for j, v in enumerate(row, start=1):
                    if v == "": continue
                    val = v
                    if i > 1:
                        s = v.replace(",", "")
                        try:
                            val = int(s) if s.lstrip("-").isdigit() else (
                                  float(s) if s.replace(".","",1).replace("-","",1).isdigit() else v)
                        except Exception: val = v
                    c = ws.cell(row=i, column=j, value=val)
                    c.font = HDR_FNT if i == 1 else BODY
                    if i == 1:
                        c.fill = HDR_FIL
                        c.alignment = Alignment(vertical="center")
        ws.freeze_panes = "A2"
        if ws.max_row > 1 and ws.max_column > 0:
            ws.auto_filter.ref = f"A1:{get_column_letter(ws.max_column)}{ws.max_row}"
        for j in range(1, min(ws.max_column, 40) + 1):
            w = max((len(str(ws.cell(row=i, column=j).value or ""))
                     for i in range(1, min(ws.max_row, 250) + 1)), default=10)
            ws.column_dimensions[get_column_letter(j)].width = min(max(w + 2, 9), 52)
        print(f"  {name:26s} {ws.max_row-1:7,d} rows x {ws.max_column} cols")
    wb.save(OUT)
    print(f"\nwrote {OUT}  ({os.path.getsize(OUT)/1e6:.1f} MB)")

if __name__ == "__main__":
    main()
