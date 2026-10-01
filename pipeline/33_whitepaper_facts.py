#!/usr/bin/env python3
"""Every quantity the white-paper section asserts, computed from the atlas tables.

Why this exists. The first draft carried its numbers as typed literals. When the RNA
assay leak was corrected (methods §7, #15) three figures silently kept the old totals,
and an earlier figure had a hardcoded "1,207" that no longer matched anything. A
document whose numbers are typed cannot be re-derived, and in an atlas that is re-
harvested it will drift the moment the data moves.

So the prose is built from this file. `33` writes docs/whitepaper_facts.json; the
document build reads it and interpolates. Re-harvest, re-run, rebuild: the text follows
the data with no hand-editing, and anything the data no longer supports fails loudly
rather than quietly persisting.

Run after 32 / 06 / 07 / 30 / 08 / 13 / 12.
"""
import csv, json, os, sys, datetime, collections, importlib.util
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import DOCS, topen
csv.field_size_limit(10**7)

# The figure script is the authority for every quantity that also appears on a panel.
# Importing it -- rather than restating its definitions here -- is what keeps the prose
# and the figures from drifting apart; the two disagreed on six quantities when this
# file carried its own copies.
_spec = importlib.util.spec_from_file_location(
    "wpfig", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "28_whitepaper_figures.py"))
WP = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(WP)

REG, PATIENT, PDERIVED = WP.REG, WP.PATIENT, WP.PDERIVED
MODALITY, project, HORIZON = WP.MODALITY, WP.project, WP.HORIZON
FAM, KEY = {}, ["chip", "meth", "acc", "cnr", "td", "sc"]
for _k, (_lab, _cls, _c) in zip(KEY, MODALITY):
    for _a in _cls:
        FAM[_a] = _k
MODERN = WP.MODERN
I, F = WP.I, WP.F
# Entity names that begin with an eponym keep their capital inside a sentence; the rest
# are lowercased when they appear mid-list. Doing this here, once, is why the document
# build needs no case logic of its own -- and why "Kaposi sarcoma" survived a pass that
# had been lowercasing it.
EPONYM = ("Kaposi", "Ewing", "Askin", "Wilms", "Merkel", "Hodgkin", "Langerhans")


def sentence_case(name):
    if name.startswith(EPONYM) or (len(name) > 1 and name[1].isupper()) \
            or name[0].isdigit() or not name[0].isalpha():
        return name
    return name[0].lower() + name[1:]


# Display names are written for a table column; a few read badly inside a sentence.
PROSE = {"Liposarcoma, pleomorphic": "pleomorphic liposarcoma",
         "Liposarcoma, well-differentiated": "well-differentiated liposarcoma",
         "Liposarcoma, NOS / mixed": "liposarcoma NOS",
         "Chondrosarcoma, primary central": "primary central chondrosarcoma",
         "Atypical fibroxanthoma / pleomorphic dermal sarcoma":
             "atypical fibroxanthoma / pleomorphic dermal sarcoma",
         "PEComa (malignant)": "malignant PEComa",
         "Desmoid / aggressive fibromatosis": "desmoid fibromatosis",
         "Fibrosarcoma (adult)": "adult fibrosarcoma",
         "FN-RMS (fusion-negative)": "fusion-negative rhabdomyosarcoma",
         "GCTB / chondroblastoma": "giant cell tumor of bone",
         "Myoepithelial carcinoma of soft tissue": "myoepithelial carcinoma of soft tissue"}
WORD = {0: "zero", 16: "sixteen", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six",
        7: "seven", 8: "eight", 9: "nine", 10: "ten", 11: "eleven", 12: "twelve",
        13: "thirteen", 14: "fourteen", 15: "fifteen", 16: "sixteen", 17: "seventeen",
        18: "eighteen", 19: "nineteen", 20: "twenty", 21: "twenty-one",
        22: "twenty-two", 23: "twenty-three", 24: "twenty-four", 25: "twenty-five",
        26: "twenty-six", 27: "twenty-seven", 28: "twenty-eight", 29: "twenty-nine",
        30: "thirty", 31: "thirty-one", 40: "forty", 44: "forty-four", 45: "forty-five"}
ORD = {1: "first", 2: "second", 3: "third", 4: "fourth", 5: "fifth",
       13: "thirteenth", 14: "fourteenth", 15: "fifteenth", 20: "twentieth",
       21: "twenty-first", 22: "twenty-second", 23: "twenty-third",
       24: "twenty-fourth", 25: "twenty-fifth", 26: "twenty-sixth"}


def main():
    f = {}
    d = WP.load()
    zcol = lambda k: sum(1 for r in d["T13"] if I(r, k) == 0)
    T13, T15, epi, reg = d["T13"], d["T15"], d["epi"], d["reg"]
    t4 = list(csv.DictReader(topen("T4_samples_atomic.tsv"), delimiter="\t"))
    T2 = list(csv.DictReader(topen("T2_models.tsv"), delimiter="\t"))
    dated = [r for r in epi if r["gse_date"][:4].isdigit()]

    # ---------- scale
    # Two counts, not one blended "studies": a GEO series and an EGA dataset are not
    # the same object, and an earlier draft summed them into a number that could not be
    # reproduced from any table.
    f["studies"] = len({r["gse"] for r in epi if r["gse"]})
    f["controlled_datasets"] = len({r["accession"] for r in T15 if r.get("accession")})
    f["samples"] = len(dated)
    f["entities"] = len(T13)
    f["reg_samples"] = len(reg)

    yr = lambda r: int(r["gse_date"][:4])
    now = max(yr(r) for r in dated)
    f["latest_year"] = now
    f["pct_last10"] = round(100 * sum(1 for r in dated if yr(r) > now - 10) / len(dated))
    f["pct_last5"] = round(100 * sum(1 for r in dated if yr(r) > now - 5) / len(dated))

    # ---------- modality mix
    fam = collections.Counter(FAM.get(r["assay_class"], "") for r in dated)
    f["fam"] = {k: fam[k] for k in ("chip", "meth", "acc", "cnr", "td", "sc")}
    f["fam_total"] = sum(f["fam"].values())
    f["pct_chip"] = round(100 * fam["chip"] / len(dated))
    f["accessibility"] = fam["acc"]; f["cutrun"] = fam["cnr"]
    f["threed"] = fam["td"]; f["singlecell"] = fam["sc"]; f["methylation"] = fam["meth"]
    f["chip"] = fam["chip"]

    def window(lo, hi):
        # Formatted exactly as the panel formats it, so the sentence and the figure
        # cannot print different integers for the same share.
        w = [r for r in dated if lo <= yr(r) <= hi]
        m = sum(1 for r in w if r["assay_class"] in MODERN)
        return f"{100 * m / len(w):.0f}" if w else "0"
    f["modern_early"] = window(2011, 2015)
    f["modern_late"] = window(2021, now)
    f["early_lo"], f["early_hi"] = 2011, 2015
    f["late_lo"], f["late_hi"] = 2021, now

    # first year each family passed fifty cumulative deposits
    arrival, run = {}, collections.Counter()
    for r in sorted(dated, key=yr):
        k = FAM.get(r["assay_class"], "")
        if not k:
            continue
        run[k] += 1
        if run[k] >= 50 and k not in arrival:
            arrival[k] = yr(r)
    f["arrival"] = arrival

    # entity coverage by year
    f["coverage"] = {}
    names = {r["display_name"] for r in T13}
    for y in (2010, 2015, 2020, now):
        seen = {r["disease"] for r in dated if yr(r) <= y}
        f["coverage"][str(y)] = sum(1 for r in T13 if r["atlas_disease"] in seen)
    f["coverage_2010"] = f["coverage"]["2010"]
    f["coverage_2015"] = f["coverage"]["2015"]
    f["coverage_2020"] = f["coverage"]["2020"]

    # ---------- projection: bootstrap over whole months, 24-month window, 6-month horizon
    ym = lambda r: (int(r["gse_date"][:4]), int(r["gse_date"].replace("-", "/")[5:7]))
    last = max(datetime.date(*ym(r), int(r["gse_date"].replace("-", "/")[8:10]))
               for r in dated)
    cutoff = datetime.date(last.year, last.month, 1)
    permo = collections.Counter(ym(r) for r in dated
                                if datetime.date(*ym(r), 1) < cutoff)
    lo, md, hi, _w = project(permo, cutoff)
    f["cutoff"] = cutoff.isoformat()
    f["cutoff_human"] = cutoff.strftime("%-d %B %Y")
    f["cutoff_month"] = cutoff.strftime("%B %Y")
    f["horizon_word"] = WORD.get(HORIZON, str(HORIZON))
    f["horizon_months"] = HORIZON
    f["proj6"] = [lo, md, hi]

    # ---------- burden
    agg = collections.defaultdict(collections.Counter)
    for r in T13:
        a = agg[r["age_class"]]
        a["cases"] += F(r, "US_cases_per_year_all_ages")
        a["reg"] += I(r, "regulatory_epigenomic_samples")
        a["epi"] += I(r, "epigenomic_samples")
    f["cases_total"] = round(sum(agg[g]["cases"] for g in agg))
    f["cases_total_round"] = int(round(f["cases_total"], -3))
    f["cases_ped"] = round(agg["pediatric"]["cases"])
    f["cases_adult"] = round(agg["adult"]["cases"])
    f["cases_both"] = round(agg["both"]["cases"])
    f["pct_cases_adult"] = round(100 * (agg["adult"]["cases"] + agg["both"]["cases"])
                                 / f["cases_total"])
    f["reg_ped"] = agg["pediatric"]["reg"]; f["reg_adult"] = agg["adult"]["reg"]
    f["epi_ped"] = agg["pediatric"]["epi"]; f["epi_adult"] = agg["adult"]["epi"]
    f["epi_rated_total"] = sum(agg[g]["epi"] for g in agg)
    # The ACS topography-coded total, which the atlas's own entity sum exceeds because
    # the atlas counts entities the registry codes to the organ of origin.
    f["cases_topography"] = 18020
    # NETSARC central-pathology-review incidence, extrapolated to the US population.
    # Typed as two inputs and a multiplication rather than as the 24,000 it produces.
    f["netsarc_rate"] = 70.7                 # per million per year [PUBMEDID: ####]
    f["us_population_m"] = 342               # US Census 2025 estimate
    f["netsarc_us_cases"] = int(round(f["netsarc_rate"] * f["us_population_m"], -3))
    f["reg_rated_total"] = sum(agg[g]["reg"] for g in agg)
    f["rate_ped"] = agg["pediatric"]["reg"] / agg["pediatric"]["cases"]
    f["rate_adult"] = agg["adult"]["reg"] / agg["adult"]["cases"]
    f["gap"] = round(f["rate_ped"] / f["rate_adult"])
    f["gap_word"] = ORD.get(f["gap"], f"{f['gap']}th")
    f["burden_x"] = round(agg["adult"]["cases"] / agg["pediatric"]["cases"])
    f["burden_x_word"] = WORD.get(f["burden_x"], str(f["burden_x"]))

    # entities with no regulatory data, and what they cost
    zero = [r for r in T13 if I(r, "regulatory_epigenomic_samples") == 0]
    f["n_zero_reg"] = len(zero)
    assert f["n_zero_reg"] == zcol("regulatory_epigenomic_samples")
    f["n_zero_reg_word"] = WORD.get(len(zero), str(len(zero)))
    f["zero_reg_cases"] = int(round(sum(F(r, "US_cases_per_year_all_ages")
                                        for r in zero), -2))
    f["zero_reg_named"] = [r["display_name"] for r in
                           sorted(zero, key=lambda r: -F(r, "US_cases_per_year_all_ages"))]
    f["n_zero_norate"] = sum(1 for r in zero if not F(r, "US_cases_per_year_all_ages"))
    f["n_zero_rated"] = f["n_zero_reg"] - f["n_zero_norate"]
    f["n_zero_rated_word"] = WORD.get(f["n_zero_rated"], str(f["n_zero_rated"]))
    f["n_zero_reg_plus1_word"] = ORD.get(f["n_zero_reg"] + 1,
                                         str(f["n_zero_reg"] + 1) + "th")

    # the two labelled points the legend disambiguates
    for key, dis in (("ratrt", "Rhabdoid tumor/ATRT"), ("fprms", "FP-RMS")):
        row = next((r for r in T13 if r["atlas_disease"] == dis), None)
        if row:
            f[key] = [round(F(row, "US_cases_per_year_all_ages")),
                      I(row, "regulatory_epigenomic_samples")]

    # ---------- composition
    def split(rows):
        cl = sum(1 for r in rows if r["sample_type"] == "cell_line")
        px = sum(1 for r in rows if r["sample_type"] in ("PDX", "organoid"))
        pt = sum(1 for r in rows if r["sample_type"] in PATIENT)
        return cl, px, pt
    for tag, rows_ in (("", reg), ("_all", epi)):
        cl, px, pt = split(rows_); tot = cl + px + pt
        f["pct_cell_line" + tag] = round(100 * cl / tot)
        f["pct_pdx" + tag] = round(100 * px / tot)
        f["pct_tumor" + tag] = round(100 * pt / tot)

    f["n_models"] = len(T2)
    f["n_problematic"] = sum(1 for r in T2 if (r.get("problematic_flag") or "").strip() == "Y")
    import re as _re
    f["n_no_rrid"] = sum(1 for r in T2
                         if not _re.search(r"CVCL_\w+", r.get("RRID_or_Cellosaurus") or ""))

    have_pd = {r["disease"] for r in reg if r["sample_type"] in PDERIVED}
    f["n_have_pd"] = sum(1 for r in T13 if r["atlas_disease"] in have_pd)
    f["n_have_pd"] = f["entities"] - zcol("patient_derived_regulatory")
    f["n_zero_pd"] = zcol("patient_derived_regulatory")
    f["n_have_pd_word"] = WORD.get(f["n_have_pd"], str(f["n_have_pd"]))
    f["n_zero_pd_word"] = WORD.get(f["n_zero_pd"], str(f["n_zero_pd"]))
    pd_growth = {}
    for y in (2015, 2020, now):
        s = {r["disease"] for r in reg
             if r["sample_type"] in PDERIVED and r["gse_date"][:4].isdigit() and yr(r) <= y}
        pd_growth[str(y)] = sum(1 for r in T13 if r["atlas_disease"] in s)
    f["pd_growth"] = pd_growth
    f["pd_2015_word"] = WORD.get(pd_growth["2015"], str(pd_growth["2015"]))
    # The entities a burden-led allocation would start with: no patient-derived
    # regulatory data, ordered by the cases each represents.
    f["zero_pd_by_burden"] = [r["display_name"] for r in
                              sorted((r for r in T13 if I(r, "patient_derived_regulatory") == 0),
                                     key=lambda r: -F(r, "US_cases_per_year_all_ages"))]
    f["pd_2020_word"] = WORD.get(pd_growth["2020"], str(pd_growth["2020"]))

    # zero-coverage counts by assay class
    f["zero_k27ac"] = zcol("H3K27ac_samples")
    f["zero_acc"] = zcol("accessibility_samples")
    f["zero_3d"] = zcol("samples_3D")
    for k in ("zero_k27ac", "zero_acc", "zero_3d"):
        f[k + "_word"] = WORD.get(f[k], str(f[k]))

    have_model = {r["disease"] for r in epi if r["sample_type"] in ("cell_line", "PDX", "organoid")}
    f["n_no_model"] = sum(1 for r in T13 if r["atlas_disease"] not in have_model)
    f["n_no_model_word"] = WORD.get(f["n_no_model"], str(f["n_no_model"]))

    # control / reference material as a share of regulatory
    ctrl = sum(1 for r in reg if r["sample_type"] == "normal/reference")
    f["pct_control"] = round(100 * ctrl / len(reg)) if reg else 0

    # ---------- validation: cell line vs patient-derived, per entity
    per = collections.defaultdict(collections.Counter)
    for r in reg:
        if r["sample_type"] == "cell_line":
            per[r["disease"]]["line"] += 1
        elif r["sample_type"] in PDERIVED:
            per[r["disease"]]["pd"] += 1
    rows = [(d2, c["line"], c["pd"]) for d2, c in per.items() if c["line"] or c["pd"]]
    name = {r["atlas_disease"]: r["display_name"] for r in T13}
    rows = [(name.get(d2, d2), ln, pd) for d2, ln, pd in rows if d2 in name]
    f["n_validation_rows"] = len(rows)
    f["n_validatable"] = sum(1 for _, ln, pd in rows if pd)
    f["n_no_patient"] = sum(1 for _, ln, pd in rows if not pd)
    for k in ("n_validation_rows", "n_validatable", "n_no_patient"):
        f[k + "_word"] = WORD.get(f[k], str(f[k]))
    ratio = sorted(((ln / pd, n) for n, ln, pd in rows if pd and ln), reverse=True)
    f["worst"] = [[n, round(v)] for v, n in ratio[:5]]
    f["best"] = [[n, round(v, 1)] for v, n in sorted(((v, n) for v, n in ratio))[:2]]
    burden = {r["display_name"]: F(r, "US_cases_per_year_all_ages") for r in T13}
    f["no_patient_named"] = [n for n, ln, pd in
                             sorted(((n, ln, pd) for n, ln, pd in rows if not pd),
                                    key=lambda t: -burden.get(t[0], 0))]

    # ---------- H3K27ac fidelity, restricted to within-study comparisons (methods #9)
    T11b = list(csv.DictReader(topen("T11b_fidelity_by_study.tsv"), delimiter="\t"))
    hi = sum(1 for r in T11b
             if float(r["model_model_reference"] or 0) > float(r["patient_patient_reference"] or 0))
    f["fidelity_hi"], f["fidelity_n"] = hi, len(T11b)
    f["fidelity_hi_word"] = WORD.get(hi, str(hi))
    f["fidelity_n_word"] = WORD.get(len(T11b), str(len(T11b)))

    # ---------- named entities the prose qualifies
    def ent(dis, col):
        r = next((x for x in T13 if x["atlas_disease"] == dis), None)
        return I(r, col) if r else 0
    f["lgfms_reg"] = ent("LGFMS/SEF", "regulatory_epigenomic_samples")
    f["lgfms_reg_word"] = WORD.get(f["lgfms_reg"], str(f["lgfms_reg"]))
    f["ess_reg"] = ent("Endometrial stromal sarcoma", "regulatory_epigenomic_samples")
    f["ess_normal"] = sum(1 for r in reg if r["disease"] == "Endometrial stromal sarcoma"
                          and r["sample_type"] == "normal/reference")
    f["ess_reg_word"] = WORD.get(f["ess_reg"], str(f["ess_reg"]))
    f["ess_normal_word"] = WORD.get(f["ess_normal"], str(f["ess_normal"]))
    f["n_corrections"] = sum(
        1 for ln in open(os.path.join(DOCS, "methods.md"), encoding="utf-8")
        if ln.startswith("| ") and ln[2:ln.find("|", 2)].strip().isdigit())
    f["n_corrections_word"] = WORD.get(f["n_corrections"], str(f["n_corrections"]))

    # ---------- access and provenance
    f["controlled_total"] = sum(I(r, "n_samples") for r in T15)
    f["controlled_reg"] = sum(I(r, "n_samples") for r in T15
                              if r["assay_family"] == "regulatory")

    g = [r for r in t4 if r["gse"] == "GSE140686"]
    f["gse140686_n"] = len(g)
    f["gse140686_labels"] = len({r["disease"] for r in g})
    onlyg = [d2 for d2 in {r["disease"] for r in epi}
             if {r["gse"] for r in epi if r["disease"] == d2} == {"GSE140686"}]
    f["gse140686_only"] = sum(1 for r in T13 if r["atlas_disease"] in onlyg)
    f["gse140686_only_word"] = WORD.get(f["gse140686_only"], str(f["gse140686_only"]))
    f["gse140686_only_named"] = sorted(name[d2] for d2 in onlyg if d2 in name)

    # Entities with no epigenomic data of any kind -- the "only named entity" claim
    seen_any = {r["disease"] for r in epi}
    f["no_epi_named"] = sorted(sentence_case(PROSE.get(r["display_name"], r["display_name"]))
                               for r in T13 if r["atlas_disease"] not in seen_any)
    f["entities_word"] = WORD.get(f["entities"], str(f["entities"]))
    for k in ("pct_last10", "pct_last5", "n_zero_norate"):
        f[k + "_word"] = WORD.get(f[k], str(f[k]))
    # Prose spellings of every named list, so the document never re-words a label.
    for k in ("zero_reg_named", "no_patient_named", "gse140686_only_named",
              "zero_pd_by_burden"):
        f[k] = [sentence_case(PROSE.get(x, x)) for x in f[k]]
    f["worst"] = [[PROSE.get(nm, nm), v] for nm, v in f["worst"]]
    out = os.path.join(DOCS, "whitepaper_facts.json")
    with open(out, "w") as fh:
        json.dump(f, fh, indent=1, sort_keys=True)
    print(f"wrote {out}  ({len(f)} facts)")
    for k in sorted(f):
        print(f"  {k:22s} {f[k]}")


if __name__ == "__main__":
    main()
