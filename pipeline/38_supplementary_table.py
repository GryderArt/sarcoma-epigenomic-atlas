#!/usr/bin/env python3
"""Supplementary Table 1 -- the entity reference.

Matt Hemming's third review ask: "a supplemental table with your list of entities ...
and provide abbreviations not listed elsewhere (e.g. ASPS)."

One row per named entity. Everything numeric is read from T13 so the table cannot
disagree with the figures; the abbreviation, the expansion and the defining lesion are
curated here because no upstream table carries a figure-caption short form.

Writes  data/S1_entity_reference.tsv
        workbooks/SupplementaryTable1_entities.xlsx
        docs/supplementary_table1.md
"""
import csv, os, importlib.util
from _paths import DATA, WORKBOOKS, DOCS, topen, twrite

# atlas_disease -> (abbreviation as it appears in a figure, full name, defining lesion)
# "" as the abbreviation means the figures spell the entity out and it needs no gloss.
REF = {
 "ASPS": ("ASPS", "Alveolar soft part sarcoma", "ASPSCR1::TFE3 t(X;17)(p11;q25)"),
 "Angiomatoid fibrous histiocytoma": ("AFH", "Angiomatoid fibrous histiocytoma",
    "EWSR1::CREB1, EWSR1::ATF1 or FUS::ATF1"),
 "Angiosarcoma": ("", "Angiosarcoma", "MYC amplification in radiation-associated and "
    "lymphoedema-associated disease; TP53, KDR, PLCG1"),
 "AFX/PDS": ("AFX / PDS", "Atypical fibroxanthoma and pleomorphic dermal sarcoma",
    "UV-signature TP53 mutation; no recurrent fusion"),
 "BCOR-sarcoma": ("", "BCOR::CCNB3 sarcoma", "BCOR::CCNB3 inv(X)(p11.4p11.22); "
    "BCOR internal tandem duplication in the infantile form"),
 "CIC-DUX4": ("", "CIC::DUX4 sarcoma", "CIC::DUX4 t(4;19) or t(10;19); "
    "CIC::FOXO4 and other variant partners"),
 "Chondrosarcoma": ("", "Chondrosarcoma", "IDH1 / IDH2 in conventional and "
    "dedifferentiated disease; COL2A1"),
 "Chordoma": ("", "Chordoma", "TBXT (brachyury) duplication and super-enhancer "
    "dependency; SMARCB1 loss in the poorly differentiated form"),
 "Clear cell sarcoma": ("CCS", "Clear cell sarcoma of soft tissue",
    "EWSR1::ATF1 t(12;22)(q13;q12); EWSR1::CREB1"),
 "CCSK": ("CCSK", "Clear cell sarcoma of the kidney",
    "BCOR internal tandem duplication; YWHAE::NUTM2"),
 "DSRCT": ("DSRCT", "Desmoplastic small round cell tumor",
    "EWSR1::WT1 t(11;22)(p13;q12)"),
 "DFSP": ("DFSP", "Dermatofibrosarcoma protuberans",
    "COL1A1::PDGFB t(17;22)(q22;q13), usually as a supernumerary ring"),
 "Desmoid": ("", "Desmoid fibromatosis (aggressive fibromatosis)",
    "CTNNB1 exon 3 mutation; APC loss in the familial form"),
 "Endometrial stromal sarcoma": ("ESS", "Endometrial stromal sarcoma",
    "JAZF1::SUZ12 in the low-grade form; YWHAE::NUTM2 and BCOR alterations "
    "in the high-grade form"),
 "EHE": ("EHE", "Epithelioid hemangioendothelioma",
    "WWTR1::CAMTA1 t(1;3)(p36;q25); YAP1::TFE3"),
 "Epithelioid sarcoma": ("", "Epithelioid sarcoma",
    "SMARCB1 (INI1) loss -- a core chromatin-remodelling lesion"),
 "Ewing": ("", "Ewing sarcoma",
    "EWSR1::FLI1 t(11;22)(q24;q12) in ~85%; EWSR1::ERG and other ETS partners"),
 "EMC": ("EMC", "Extraskeletal myxoid chondrosarcoma",
    "NR4A3::EWSR1 t(9;22)(q22;q12); NR4A3::TAF15"),
 "FN-RMS": ("FN-RMS", "Fusion-negative rhabdomyosarcoma",
    "RAS-pathway mutation (NRAS, KRAS, HRAS), FGFR4, TP53; no FOXO1 fusion"),
 "FP-RMS": ("FP-RMS", "Fusion-positive rhabdomyosarcoma",
    "PAX3::FOXO1 t(2;13)(q35;q14) or PAX7::FOXO1 t(1;13)(p36;q14)"),
 "Fibrosarcoma NOS": ("", "Fibrosarcoma of adults",
    "a diagnosis of exclusion; no defining lesion"),
 "GCTB/Chondroblastoma": ("GCTB", "Giant cell tumor of bone, with chondroblastoma",
    "H3-3A p.G34W in GCTB and H3-3B p.K36M in chondroblastoma -- oncohistones"),
 "GIST": ("GIST", "Gastrointestinal stromal tumor",
    "KIT or PDGFRA mutation; SDH deficiency in the wild-type group, which is "
    "a DNA-hypermethylation phenotype"),
 "Infantile fibrosarcoma": ("", "Infantile fibrosarcoma",
    "ETV6::NTRK3 t(12;15)(p13;q25)"),
 "IMT": ("IMT", "Inflammatory myofibroblastic tumor",
    "ALK fusion in ~50%; ROS1, NTRK3 and RET fusions"),
 "Intimal sarcoma": ("", "Intimal sarcoma", "MDM2 amplification with CDK4 and PDGFRA"),
 "Kaposi sarcoma": ("KS", "Kaposi sarcoma",
    "human herpesvirus 8 (HHV-8 / KSHV) latent infection"),
 "LGFMS/SEF": ("LGFMS / SEF", "Low-grade fibromyxoid sarcoma and sclerosing "
    "epithelioid fibrosarcoma", "FUS::CREB3L2 t(7;16)(q33;p11); FUS::CREB3L1; "
    "EWSR1::CREB3L1 in SEF"),
 "Leiomyosarcoma": ("LMS", "Leiomyosarcoma",
    "TP53, RB1 and PTEN loss with widespread copy-number change; "
    "alternative lengthening of telomeres"),
 "Liposarcoma-NOS": ("LPS, NOS", "Liposarcoma, not otherwise specified or mixed",
    "heterogeneous; see the subtype rows"),
 "Liposarcoma-dediff": ("LPS, DD", "Dedifferentiated liposarcoma",
    "MDM2 and CDK4 amplification on a 12q13-15 neochromosome"),
 "Liposarcoma-myxoid": ("LPS, myxoid", "Myxoid and round cell liposarcoma",
    "FUS::DDIT3 t(12;16)(q13;p11); EWSR1::DDIT3"),
 "Liposarcoma-pleomorphic": ("LPS, pleo", "Pleomorphic liposarcoma",
    "TP53 and RB1 loss; complex karyotype without MDM2 amplification"),
 "Liposarcoma-WD": ("LPS, WD", "Well-differentiated liposarcoma / atypical "
    "lipomatous tumor", "MDM2 and CDK4 amplification on a 12q13-15 "
    "supernumerary ring or giant marker chromosome"),
 "MPNST": ("MPNST", "Malignant peripheral nerve sheath tumor",
    "NF1 and CDKN2A loss with PRC2 inactivation (EED or SUZ12), which erases "
    "H3K27me3 -- a defining chromatin lesion"),
 "RMS-MYOD1": ("MYOD1-RMS", "MYOD1-mutant rhabdomyosarcoma, usually spindle cell "
    "or sclerosing", "MYOD1 p.L122R, frequently with PIK3CA"),
 "Myoepithelial carcinoma": ("", "Myoepithelial carcinoma of soft tissue",
    "EWSR1::POU5F1, EWSR1::PBX1, EWSR1::ZNF444; SMARCB1 loss in a subset"),
 "Myxofibrosarcoma": ("MFS", "Myxofibrosarcoma",
    "complex karyotype; TP53 and CDKN2A loss; no recurrent fusion"),
 "Ossifying fibromyxoid tumor": ("OFMT", "Ossifying fibromyxoid tumor",
    "PHF1 rearrangement (EP400::PHF1, MEAF6::PHF1) -- a PRC2-component fusion"),
 "Osteosarcoma": ("OS", "Osteosarcoma",
    "TP53 and RB1 loss with chromothripsis and kataegis; no defining fusion"),
 "PEComa": ("PEComa", "Perivascular epithelioid cell tumor, malignant",
    "TSC1 or TSC2 loss; TFE3 rearrangement in a distinct subset"),
 "Rhabdoid tumor/ATRT": ("ATRT", "Malignant rhabdoid tumor and atypical "
    "teratoid/rhabdoid tumor", "biallelic SMARCB1 (INI1) loss; SMARCA4 loss in a "
    "small subset -- the canonical SWI/SNF sarcoma"),
 "Solitary fibrous tumor": ("SFT", "Solitary fibrous tumor",
    "NAB2::STAT6 inv(12)(q13q13)"),
 "Synovial sarcoma": ("SS", "Synovial sarcoma",
    "SS18::SSX1, SS18::SSX2 or SS18::SSX4 t(X;18)(p11;q11), which rewires "
    "the BAF complex"),
 "UPS/MFH": ("UPS / MFH", "Undifferentiated pleomorphic sarcoma, formerly "
    "malignant fibrous histiocytoma",
    "complex karyotype; TP53, RB1 and ATRX loss; no defining lesion"),
}

# The label a reader actually sees on a patch or an axis comes from stage 35. Import it
# rather than retyping it: a hand-copied label is a label that drifts.
_spec = importlib.util.spec_from_file_location(
    "quilt", os.path.join(os.path.dirname(os.path.abspath(__file__)), "35_burden_quilt.py"))
_q = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_q)
fig_label = _q.short

AGE = {"pediatric": "Pediatric-predominant", "adult": "Adult-predominant",
       "both": "Both age ranges"}

COLS = [
  ("abbrev",      "Abbreviation (as printed in the figures)"),
  ("full_name",   "Entity"),
  ("also_known",  "Other common abbreviation"),
  ("lesion",      "Defining molecular lesion"),
  ("age",         "Age class"),
  ("cases",       "US cases per year (all ages)"),
  ("cases_range", "Plausible range"),
  ("basis",       "Basis for the incidence estimate"),
  ("epi",         "Epigenomic samples"),
  ("reg",         "Regulatory samples"),
  ("meth",        "DNA methylation samples"),
  ("pd_reg",      "Regulatory samples on patient-derived material"),
  ("pt_epi",      "Epigenomic samples on primary tumor"),
  ("models",      "Distinct models"),
  ("studies",     "Studies"),
  ("per_case",    "Epigenomic samples per new case per year"),
  ("figs",        "Appears in"),
  ("note",        "Note"),
]

def I(r, k):
    try: return int(float(r[k] or 0))
    except (ValueError, TypeError): return 0

def build():
    T13 = list(csv.DictReader(topen("T13_incidence_vs_data.tsv"), delimiter="\t"))
    rows = []
    for r in T13:
        d = r["atlas_disease"]
        if d not in REF:
            raise SystemExit("no reference entry for " + d)   # a new entity must be curated
        ab, full, lesion = REF[d]
        label = fig_label(r["display_name"])
        cases = I(r, "US_cases_per_year_all_ages")
        reg, pd_reg = I(r, "regulatory_epigenomic_samples"), I(r, "patient_derived_regulatory")
        epi = I(r, "epigenomic_samples")

        # Which panels the entity reaches, and -- when it is absent from one -- why.
        # Stating the reason is the point of the table: Matt's confusion was over
        # entities appearing in one figure and not another.
        figs, note = [], ""
        figs.append("Fig. 2A-C")
        if cases > 0:
            figs.append("Fig. 2D")
        else:
            note = ("No population incidence rate is published anywhere and none can be "
                    "derived, so the entity has no x-position in Fig. 2D and no area in "
                    "Fig. 5; it is named beneath both.")
        figs.append("Fig. 3B-C")
        line = reg - pd_reg
        if d == "LGFMS/SEF":
            note = ("Its regulatory data is from engineered non-tumor cells, so the "
                    "cell-line-versus-patient comparison in Fig. 4 would be misleading "
                    "and the entity is withheld from that figure only.")
        elif reg > 0 and (line > 0 or pd_reg > 0):
            figs.append("Fig. 4")
        elif reg == 0:
            # An entity can fail two tests at once -- no published rate AND no
            # regulatory data -- and the reader needs both reasons, not the first one.
            z = ("No regulatory epigenomic data of any kind, so the entity is absent "
                 "from Fig. 4" + (" as well." if not cases else
                 " and sits on the floor of Fig. 2D."))
            note = (note + " " + z) if note else z
        if cases > 0:
            figs.append("Fig. 5")

        rows.append({
          "abbrev": label if label != full else "—",
          "full_name": full,
          "also_known": (ab if ab and ab != label else ""),
          "lesion": lesion,
          "age": AGE[r["age_class"]],
          "cases": cases if cases else "not established",
          "cases_range": (f"{I(r,'US_low'):,}–{I(r,'US_high'):,}"
                          if cases else "—"),
          "basis": r["incidence_basis"],
          "epi": epi, "reg": reg,
          "meth": I(r, "DNA_methylation_samples"),
          "pd_reg": pd_reg,
          "pt_epi": I(r, "primary_tumor_epigenomic"),
          "models": I(r, "distinct_models"),
          "studies": I(r, "n_studies"),
          "per_case": (round(epi / cases, 3) if cases else "—"),
          "figs": ", ".join(figs),
          "note": note,
        })
    rows.sort(key=lambda x: x["full_name"].lower())
    return rows

def write_tsv(rows):
    out = os.path.join(DATA, "S1_entity_reference.tsv")
    with twrite(out) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([h for _, h in COLS])
        for r in rows: w.writerow([r[k] for k, _ in COLS])
    return out

def write_xlsx(rows):
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
        from openpyxl.utils import get_column_letter
    except ImportError:
        print("openpyxl missing -- skipping the workbook"); return None
    wb = Workbook(); ws = wb.active; ws.title = "Supplementary Table 1"
    head = Font(name="Arial", size=9, bold=True, color="FFFFFF")
    body = Font(name="Arial", size=9)
    fill = PatternFill("solid", fgColor="1F3A5F")
    thin = Side(style="thin", color="D9D9D9")
    ws.append([h for _, h in COLS])
    for c in ws[1]:
        c.font = head; c.fill = fill
        c.alignment = Alignment(vertical="center", wrap_text=True)
    for r in rows:
        ws.append([r[k] for k, _ in COLS])
    widths = [13, 44, 30, 46, 20, 13, 14, 60, 11, 11, 13, 16, 14, 10, 9, 13, 34, 62]
    for i, wdt in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(i)].width = wdt
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.font = body
            c.border = Border(bottom=thin)
            c.alignment = Alignment(vertical="top",
                                    wrap_text=c.column in (4, 8, 17, 18),
                                    horizontal="right" if c.column in
                                        (6, 9, 10, 11, 12, 13, 14, 15, 16) else "left")
    ws.freeze_panes = "C2"
    ws.row_dimensions[1].height = 42
    ws.auto_filter.ref = ws.dimensions
    out = os.path.join(WORKBOOKS, "SupplementaryTable1_entities.xlsx")
    os.makedirs(WORKBOOKS, exist_ok=True)
    wb.save(out)
    return out

def write_md(rows):
    keep = ["abbrev", "full_name", "lesion", "age", "cases", "epi", "reg", "pd_reg", "figs"]
    hdr = {k: h for k, h in COLS}
    lines = ["# Supplementary Table 1. Entity reference", "",
             "Every named entity in the atlas, the abbreviation used for it in the "
             "figures, its defining lesion, its annual US incidence, the epigenomic "
             "data held for it, and the figures it appears in.", "",
             "| " + " | ".join(hdr[k] for k in keep) + " |",
             "|" + "|".join("---" for _ in keep) + "|"]
    for r in rows:
        lines.append("| " + " | ".join(str(r[k]).replace("|", "\\|") for k in keep) + " |")
    out = os.path.join(DOCS, "supplementary_table1.md")
    with twrite(out) as fh:
        fh.write("\n".join(lines) + "\n")
    return out

if __name__ == "__main__":
    rows = build()
    for f in (write_tsv(rows), write_xlsx(rows), write_md(rows)):
        if f: print("wrote", f)
    print(f"{len(rows)} entities; "
          f"{sum(1 for r in rows if 'Fig. 4' in r['figs'])} reach Fig. 4; "
          f"{sum(1 for r in rows if r['note'])} carry a reason for absence")
