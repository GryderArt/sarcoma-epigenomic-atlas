#!/usr/bin/env python3
"""Adult-vs-paediatric burden and gap figures for the SASS white paper.
Arial throughout, TrueType-embedded (pdf.fonttype=42) so all text stays editable."""
import csv, collections, math, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Patch, FancyBboxPatch
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
FIG = FIGURES

fm._load_fontmanager(try_read_cache=False)
plt.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.04,
})
O = DATA
csv.field_size_limit(10**7)

PED, BOTH, ADULT = "#2a78d6", "#1baf7a", "#eb6834"
CLS = {"paediatric": PED, "both": BOTH, "adult": ADULT}
CLSLAB = {"paediatric": "Paediatric-predominant", "both": "Both ages",
          "adult": "Adult-predominant"}
GREY, DARK, RED = "#b9c2cc", "#2b3440", "#c0392b"
METH, REG = "#9dc3ee", "#0d366b"
K = FuncFormatter(lambda v, p: f"{int(v):,}")

T13 = [r for r in csv.DictReader(topen("T13_incidence_vs_data.tsv"), delimiter="\t")]
for r in T13:
    for k in ("US_cases_per_year_all_ages","US_low","US_high","epigenomic_samples",
              "regulatory_epigenomic_samples","DNA_methylation_samples","H3K27ac_samples",
              "accessibility_samples","samples_3D","primary_tumour_epigenomic",
              "primary_tumour_regulatory","patient_derived_regulatory","n_studies",
              "distinct_models","distinct_pdx",
              "models_with_H3K27ac","patients_with_H3K27ac"):
        r[k] = int(r[k] or 0)
    for k in ("pct_from_largest_study","pct_DNA_methylation"):
        r[k] = float(r[k]) if r[k] else 0.0
RATED = [r for r in T13 if r["US_cases_per_year_all_ages"] > 0]

def save(fig, name):
    fig.savefig(f"{FIG}/{name}.pdf"); plt.close(fig); print("  wrote", name)

def tidy(ax):
    ax.tick_params(length=3, pad=2)
    for s in ("left","bottom"):
        if s in ax.spines: ax.spines[s].set_color("#5a6572")

# ---------------------------------------------------------------- F18
def f18():
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.75),
                             gridspec_kw={"width_ratios":[1,1,0.8], "wspace":0.5})
    groups = ["paediatric","both","adult"]
    cases = {g: sum(r["US_cases_per_year_all_ages"] for r in RATED if r["age_class"]==g)
             for g in groups}
    reg   = {g: sum(r["regulatory_epigenomic_samples"] for r in RATED if r["age_class"]==g)
             for g in groups}
    allep = {g: sum(r["epigenomic_samples"] for r in RATED if r["age_class"]==g) for g in groups}
    y = list(range(3))

    ax = axes[0]
    ax.barh(y, [cases[g] for g in groups], color=[CLS[g] for g in groups], height=0.6,
            edgecolor="none")
    for i, g in enumerate(groups):
        ax.text(cases[g]*1.04, i, f"{cases[g]:,}", va="center", fontsize=7.5, color=DARK)
    ax.set_yticks(y); ax.set_yticklabels([CLSLAB[g] for g in groups], fontsize=7.5)
    ax.tick_params(axis="y", length=0)
    ax.invert_yaxis(); ax.set_xlim(0, max(cases.values())*1.30)
    ax.set_xticks([0, 5000, 10000]); ax.xaxis.set_major_formatter(K)
    ax.set_xlabel("US cases per year", fontsize=8)
    ax.set_title("Disease burden", fontsize=8.5, weight="bold", pad=14, loc="left")
    tidy(ax)

    ax = axes[1]
    ax.barh(y, [reg[g] for g in groups], color=[CLS[g] for g in groups], height=0.6,
            edgecolor="none")
    ax.barh(y, [allep[g]-reg[g] for g in groups], left=[reg[g] for g in groups],
            color=[CLS[g] for g in groups], height=0.6, alpha=0.30, edgecolor="none")
    for i, g in enumerate(groups):
        ax.text(allep[g]*1.04, i, f"{allep[g]:,}", va="center", fontsize=7.5, color=DARK)
    ax.set_yticks(y); ax.set_yticklabels([]); ax.tick_params(axis="y", length=0)
    ax.invert_yaxis(); ax.set_xlim(0, max(allep.values())*1.34)
    ax.set_xticks([0, 5000, 10000]); ax.xaxis.set_major_formatter(K)
    ax.set_xlabel("Epigenomic samples in public archives", fontsize=8)
    ax.set_title("Data generated", fontsize=8.5, weight="bold", pad=14, loc="left")
    ax.legend(handles=[Patch(facecolor=DARK, label="Regulatory (ChIP / ATAC / 3D)"),
                       Patch(facecolor=DARK, alpha=0.30, label="DNA methylation only")],
              fontsize=6.4, frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.0),
              ncol=2, handlelength=1.0, handleheight=0.75, columnspacing=1.0,
              borderpad=0.0, handletextpad=0.4)
    tidy(ax)

    ax = axes[2]
    rr = [reg[g]/cases[g] for g in groups]
    ax.bar(range(3), rr, color=[CLS[g] for g in groups], width=0.6, edgecolor="none")
    for i, v in enumerate(rr):
        ax.text(i, v*1.05, f"{v:.2f}", ha="center", va="bottom", fontsize=7.5, color=DARK)
    ax.set_xticks(range(3)); ax.set_xticklabels(["Paediatric","Both","Adult"], fontsize=7.5)
    ax.set_ylabel("Regulatory epigenomes per\nnew case per year", fontsize=8)
    ax.set_ylim(0, max(rr)*1.28)
    ax.set_title(f"{rr[0]/rr[2]:.0f}× gap", fontsize=8.5, weight="bold", pad=14, loc="left")
    tidy(ax)

    fig.suptitle("Adult sarcomas carry ~5× the disease burden and a small fraction of "
                 "the epigenomic data", fontsize=9.2, weight="bold", y=1.10, x=0.0, ha="left")
    save(fig, "F18_adult_vs_paediatric_burden")

# ---------------------------------------------------------------- F19
def f19():
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    FLOOR = 0.42
    zero = sorted([r for r in RATED if r["regulatory_epigenomic_samples"] == 0],
                  key=lambda r: -r["US_cases_per_year_all_ages"])
    have = [r for r in RATED if r["regulatory_epigenomic_samples"] > 0]
    ax.axhspan(0.30, 0.60, color="#fbe9e4", zorder=0)
    for r in have:
        ax.scatter(r["US_cases_per_year_all_ages"], r["regulatory_epigenomic_samples"],
                   s=40, color=CLS[r["age_class"]], edgecolor="white", linewidth=0.6, zorder=3)
    for r in zero:
        ax.scatter(r["US_cases_per_year_all_ages"], FLOOR, s=44, color=RED, marker="v",
                   edgecolor="white", linewidth=0.6, zorder=4)
    ax.text(26000, 0.42, "ZERO", fontsize=7, color=RED, weight="bold",
            va="center", ha="right")

    # de-collide labels for the notable points
    show = [r for r in have if r["US_cases_per_year_all_ages"] > 900
            or r["regulatory_epigenomic_samples"] > 900
            or r["US_cases_per_year_all_ages"] < 45]
    placed = []            # (log10 x, log10 y) of each label anchor already used
    CAND = [(0, 7, "center"), (0, -11, "center"), (-7, -2.5, "right"), (7, -2.5, "left"),
            (0, 15, "center"), (0, -19, "center")]
    for r in sorted(show, key=lambda r: -r["regulatory_epigenomic_samples"]):
        x, yv = r["US_cases_per_year_all_ages"], r["regulatory_epigenomic_samples"]
        lx, ly = math.log10(x), math.log10(yv)
        for dx, dy, ha in CAND:
            ax_ = lx + dx*0.011 + (0.13 if ha == "left" else -0.13 if ha == "right" else 0)
            ay_ = ly + dy*0.016
            if all(abs(ax_-px) > 0.30 or abs(ay_-py) > 0.115 for px, py in placed):
                placed.append((ax_, ay_)); break
        else:
            dx, dy, ha = CAND[0]; placed.append((lx, ly+0.11))
        ax.annotate(r["display_name"].split(" (")[0], (x, yv), fontsize=6.2, color=DARK,
                    xytext=(dx, dy), textcoords="offset points", ha=ha,
                    va="center" if ha != "center" else ("bottom" if dy > 0 else "top"),
                    zorder=5)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(9, 30000); ax.set_ylim(0.29, 200000)
    ax.set_xlabel("US cases per year (all ages)", fontsize=8.5)
    ax.set_ylabel("Regulatory epigenomic samples in public archives\n"
                  "(ChIP-seq / CUT&RUN / ATAC / DNase / Hi-C / HiChIP)", fontsize=8.5)
    ax.set_yticks([FLOOR,1,10,100,1000,10000])
    ax.set_yticklabels(["0","1","10","100","1,000","10,000"])
    ax.set_ylim(0.29, 200000)
    ax.set_xticks([10,100,1000,10000]); ax.set_xticklabels(["10","100","1,000","10,000"])
    ax.grid(True, which="major", color="#eef1f4", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(handles=[plt.Line2D([],[],marker="o",ls="",color=CLS[g],markeredgecolor="white",
                                  label=CLSLAB[g], markersize=6)
                       for g in ("paediatric","both","adult")],
              fontsize=6.8, frameon=False, loc="lower right", handletextpad=0.2)

    txt = ("No regulatory epigenomics whatsoever\n" +
           "\n".join(f"· {r['display_name'].split(' (')[0][:38]}  "
                     f"({r['US_cases_per_year_all_ages']:,}/yr)" for r in zero))
    ax.text(0.012, 0.988, txt, transform=ax.transAxes, fontsize=6.1, color=RED,
            va="top", ha="left", linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.42", facecolor="#fdf3f0", edgecolor="#f0c8bd",
                      linewidth=0.6))
    ax.set_title("Regulatory epigenomics tracks research history, not disease burden",
                 fontsize=9.2, weight="bold", loc="left", pad=8)
    tidy(ax)
    save(fig, "F19_burden_vs_regulatory_epigenomics")

# ---------------------------------------------------------------- F20
def f20():
    d = sorted(RATED, key=lambda r: -r["US_cases_per_year_all_ages"])
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 6.2),
                             gridspec_kw={"width_ratios":[0.80,1.12,1.20], "wspace":0.0})
    y = list(range(len(d)))

    ax = axes[0]
    ax.barh(y, [r["US_cases_per_year_all_ages"] for r in d],
            color=[CLS[r["age_class"]] for r in d], height=0.66, edgecolor="none")
    ax.invert_xaxis(); ax.invert_yaxis()
    ax.set_yticks([]); ax.set_xticks([0,2000,4000]); ax.xaxis.set_major_formatter(K)
    ax.set_xlabel("US cases per year", fontsize=8)
    ax.set_title("Burden", fontsize=8.5, weight="bold", loc="right", pad=6)
    ax.spines["left"].set_visible(False)
    tidy(ax)

    ax = axes[1]; ax.axis("off"); ax.set_ylim(len(d)-0.5, -0.5); ax.set_xlim(0, 1)
    for i, r in enumerate(d):
        ax.text(0.5, i, r["display_name"], fontsize=6.3, ha="center", va="center",
                color=CLS[r["age_class"]] if r["regulatory_epigenomic_samples"] else RED)

    ax = axes[2]
    reg  = [r["regulatory_epigenomic_samples"] for r in d]
    meth = [r["DNA_methylation_samples"] for r in d]
    ax.barh(y, reg, color=REG, height=0.66, edgecolor="none", label="Regulatory epigenomics")
    ax.barh(y, meth, left=reg, color=METH, height=0.66, edgecolor="none",
            label="DNA methylation only")
    ax.invert_yaxis(); ax.set_yticks([])
    for i, r in enumerate(d):
        tot = r["epigenomic_samples"]
        if r["regulatory_epigenomic_samples"] == 0:
            ax.text(tot*1.55 + 2, i, f"{tot:,}  ·  DNA methylation only",
                    fontsize=5.9, color=RED, va="center", weight="bold")
        elif tot:
            ax.text(tot*1.18, i, f"{tot:,}", fontsize=5.9, color="#6b7580", va="center")
    ax.set_xscale("log"); ax.set_xlim(1, 60000)
    ax.set_xticks([1,10,100,1000,10000]); ax.set_xticklabels(["1","10","100","1,000","10,000"])
    ax.set_xlabel("Epigenomic samples (log scale)", fontsize=8)
    ax.set_title("Data, split by what it can answer", fontsize=8.5, weight="bold",
                 loc="left", pad=6)
    ax.legend(fontsize=6.6, frameon=False, loc="lower right", handlelength=1.0,
              handleheight=0.75)
    ax.spines["left"].set_visible(False)
    tidy(ax)

    axes[0].legend(handles=[Patch(facecolor=CLS[g], label=CLSLAB[g])
                            for g in ("paediatric","both","adult")],
                   fontsize=6.6, frameon=False, loc="lower left", handlelength=1.0,
                   handleheight=0.75)
    fig.suptitle("Nine sarcoma entities with a published incidence rate — two of them above "
                 "1,300 US cases a year — have DNA methylation profiling and nothing else",
                 fontsize=9.0, weight="bold", y=0.925, x=0.0, ha="left")
    save(fig, "F20_methylation_only_entities")

# ---------------------------------------------------------------- F21
def f21():
    d = sorted(T13, key=lambda r: -r["US_cases_per_year_all_ages"])
    zero = [r for r in d if r["patient_derived_regulatory"] == 0]
    some = [r for r in d if r["patient_derived_regulatory"] > 0]
    GAP = 1.1
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    pos, lab, col = [], [], []
    for i, r in enumerate(zero):
        pos.append(i); lab.append(r["display_name"]); col.append(RED)
    for j, r in enumerate(some):
        pos.append(len(zero) + GAP + j); lab.append(r["display_name"])
        col.append(CLS[r["age_class"]])
    for p, r in zip(pos, zero + some):
        n = r["patient_derived_regulatory"]
        if n == 0:
            ax.scatter(1.0, p, marker="x", s=24, color=RED, linewidth=1.2, zorder=3)
            c = r["US_cases_per_year_all_ages"]
            ax.text(1.35, p, f"none, against ~{c:,} US cases per year" if c
                             else "none  (no population rate published anywhere)",
                    fontsize=6.2, color=RED, va="center")
        else:
            ax.barh(p, n, color=CLS[r["age_class"]], height=0.66, edgecolor="none")
            ax.text(n*1.14, p, f"{n:,}", fontsize=6.2, color="#6b7580", va="center")
    ax.set_yticks(pos); ax.set_yticklabels(lab, fontsize=6.4)
    for t, c in zip(ax.get_yticklabels(), col): t.set_color(c)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.invert_yaxis()
    ax.set_xscale("log"); ax.set_xlim(0.9, 4000)
    ax.set_xticks([1,10,100,1000]); ax.set_xticklabels(["1","10","100","1,000"])
    ax.set_xlabel("Regulatory epigenomic profiles of PATIENT-DERIVED MATERIAL — tumour "
                  "tissue, PDX or patient-derived organoid\n(ChIP-seq / CUT&RUN / ATAC / Hi-C; "
                  "input and IgG controls excluded)", fontsize=8)
    ax.axhline(len(zero) + GAP/2 - 0.5, color="#cfd6dd", linewidth=0.8, linestyle=(0,(3,3)))
    ax.spines["left"].set_visible(False)
    ax.set_title(f"{len(zero)} of {len(T13)} sarcoma entities have never had a single "
                 f"patient's tumour profiled for\nactive chromatin, chromatin accessibility "
                 f"or 3D genome architecture",
                 fontsize=9.2, weight="bold", loc="left", pad=8)
    tidy(ax)
    save(fig, "F21_zero_primary_tumour_regulatory")

# ---------------------------------------------------------------- F22
def f22():
    d = [r for r in T13 if r["epigenomic_samples"] > 0 and r["pct_from_largest_study"] >= 50]
    d = sorted(d, key=lambda r: (-r["pct_from_largest_study"], -r["epigenomic_samples"]))[:26]
    fig, ax = plt.subplots(figsize=(7.2, 4.9))
    y = list(range(len(d)))
    cols = [RED if r["largest_study"] == "GSE140686" else GREY for r in d]
    ax.barh(y, [r["pct_from_largest_study"] for r in d], color=cols, height=0.66,
            edgecolor="none")
    for i, r in enumerate(d):
        ax.text(r["pct_from_largest_study"] + 1.5, i,
                f"{r['largest_study']}  ·  {r['epigenomic_samples']:,} samples from "
                f"{r['n_studies']} stud{'y' if r['n_studies']==1 else 'ies'}",
                fontsize=6.0, va="center", color=DARK)
    ax.set_yticks(y); ax.set_yticklabels([r["display_name"] for r in d], fontsize=6.4)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.invert_yaxis(); ax.set_xlim(0, 178)
    ax.set_xticks([0,25,50,75,100]); ax.set_xticklabels(["0","25","50","75","100%"])
    ax.set_xlabel("Share of the entity's entire public epigenome contributed by "
                  "one single study", fontsize=8)
    ax.legend(handles=[Patch(facecolor=RED, label="GSE140686 — 1,505 arrays deposited with "
                                                  "the diagnosis stripped from every sample"),
                       Patch(facecolor=GREY, label="Some other single study")],
              fontsize=6.5, frameon=False, loc="upper left",
              bbox_to_anchor=(0.0, -0.10), ncol=1, handlelength=1.0, handleheight=0.75)
    ax.spines["left"].set_visible(False)
    ax.set_title("Most sarcoma entities rest on one deposit — and for many of them that "
                 "deposit\ncarries no entity label in GEO at all",
                 fontsize=9.2, weight="bold", loc="left", pad=8)
    tidy(ax)
    save(fig, "F22_single_study_dependence")

if __name__ == "__main__":
    print("figures ->", FIG)
    f18(); f19(); f20(); f21(); f22()
    print("done")

# ---------------------------------------------------------------- F23
def f23():
    """The national-level version of the point, built on primary registry sources
    rather than on this atlas's curated entity subset."""
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1),
                             gridspec_kw={"width_ratios":[1, 1.15], "wspace":0.42})

    ax = axes[0]
    grp = ["Children and\nadolescents (0-19)", "Adults (20+)"]
    share = [11.5, 0.8]          # sarcoma as % of all malignancies in that age group
    ax.bar([0,1], share, color=[PED, ADULT], width=0.55, edgecolor="none")
    for i, v in enumerate(share):
        ax.text(i, v+0.35, f"{v:.1f}%", ha="center", fontsize=8, weight="bold", color=DARK)
    ax.set_xticks([0,1]); ax.set_xticklabels(grp, fontsize=7.5)
    ax.set_ylim(0, 14); ax.set_ylabel("Sarcoma as a share of all\ncancers in that age group",
                                      fontsize=8)
    ax.set_title("Why sarcoma reads as a paediatric disease", fontsize=8.5, weight="bold",
                 loc="left", pad=6)
    tidy(ax)

    ax = axes[1]
    cases = {"paediatric": 1908, "adult": 16112}      # US, site-coded, ACS 2026
    reg   = {g: sum(r["regulatory_epigenomic_samples"] for r in RATED if r["age_class"]==g)
             for g in ("paediatric","both","adult")}
    reg["adult"] += reg.pop("both")                   # shared entities counted with adults
    tot_c = sum(cases.values()); tot_r = sum(reg.values())
    bars = [("US sarcoma cases\nper year", [100*cases["paediatric"]/tot_c,
                                            100*cases["adult"]/tot_c]),
            ("Regulatory epigenomic\nsamples ever generated", [100*reg["paediatric"]/tot_r,
                                                               100*reg["adult"]/tot_r])]
    for i, (lab, vals) in enumerate(bars):
        ax.barh(i, vals[0], color=PED, height=0.34, edgecolor="none")
        ax.barh(i, vals[1], left=vals[0], color=ADULT, height=0.34, edgecolor="none")
        ax.text(vals[0]/2, i, f"{vals[0]:.0f}%", ha="center", va="center", fontsize=7.5,
                color="white", weight="bold")
        ax.text(vals[0]+vals[1]/2, i, f"{vals[1]:.0f}%", ha="center", va="center",
                fontsize=7.5, color="white", weight="bold")
    ax.set_yticks([0,1]); ax.set_yticklabels([b[0] for b in bars], fontsize=7.5)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(1.75, -0.75)
    ax.set_xlim(0, 100); ax.set_xticks([0,25,50,75,100])
    ax.set_xticklabels(["0","25","50","75","100%"])
    ax.set_title("Where the cases are, and where the data is", fontsize=8.5, weight="bold",
                 loc="left", pad=6)
    ax.legend(handles=[Patch(facecolor=PED, label="Children and adolescents / "
                                                  "paediatric-predominant entities"),
                       Patch(facecolor=ADULT, label="Adults / adult- and both-age entities")],
              fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(0.0, -0.24))
    ax.spines["left"].set_visible(False)
    tidy(ax)

    fig.text(0.0, -0.30,
             "Sources: sarcoma share of all cancers and US case counts, Siegel et al., Cancer "
             "Statistics 2026 (PMID 41528114; 4,110 bone + 13,910 soft tissue = 18,020) and "
             "SEER 21 age distribution;\npaediatric share of childhood malignancy, Siegel DA "
             "et al., JNCI 2023 (PMID 37433078; ICCC groups VIII + IX = 11.5%). Sample counts "
             "from this atlas (T13).", fontsize=5.8, color="#6b7580", ha="left", va="top")
    fig.suptitle("Sarcoma is 1% of adult cancer and 11% of childhood cancer — and the "
                 "epigenomic literature followed the 11%",
                 fontsize=9.2, weight="bold", y=1.06, x=0.0, ha="left")
    save(fig, "F23_attention_follows_paediatric")

f23()

# ---------------------------------------------------------------- F24
def f24():
    """Three states a GEO-only survey conflates: open, locked behind a DAC, or absent."""
    ctrl = list(csv.DictReader(topen("T15_controlled_access.tsv"),
                               delimiter="\t"))
    # Entities whose EGA records do not state a subtype cannot be fanned out across the
    # atlas's subtype bins without inventing an attribution. They are reported as a
    # footnote instead.
    CMAP = {"UPS / undifferentiated sarcoma": "UPS/MFH"}
    UNRESOLVED = {"RMS (any subtype)", "Liposarcoma (any subtype)",
                  "Sarcoma, unspecified / mixed cohort"}
    creg, cdna = collections.Counter(), collections.Counter()
    unres = collections.Counter()
    for c in ctrl:
        n = int(c["n_samples"] or 0)
        if c["entity"] in UNRESOLVED:
            if c["assay_family"] == "regulatory": unres[c["entity"]] += n
            continue
        e = CMAP.get(c["entity"], c["entity"])
        (creg if c["assay_family"] == "regulatory" else cdna)[e] += n

    d = [r for r in T13 if r["US_cases_per_year_all_ages"] > 0]
    d.sort(key=lambda r: -r["US_cases_per_year_all_ages"])
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    y = list(range(len(d)))
    OPEN, LOCK = "#0d366b", "#eda100"
    for i, r in enumerate(d):
        pub = r["regulatory_epigenomic_samples"]
        lock = creg.get(r["atlas_disease"], 0)
        if pub: ax.barh(i, pub, color=OPEN, height=0.64, edgecolor="none")
        if lock: ax.barh(i, lock, left=max(pub, 0), color=LOCK, height=0.64, edgecolor="none")
        tot = pub + lock
        if tot == 0:
            ax.scatter(1.1, i, marker="x", s=22, color=RED, linewidth=1.2, zorder=3)
            ax.text(1.6, i, "nothing exists, open or closed", fontsize=6.0, color=RED,
                    va="center", weight="bold")
        else:
            bits = []
            if pub: bits.append(f"{pub:,} open")
            if lock: bits.append(f"{lock:,} controlled")
            ax.text(tot*1.13, i, "  ·  ".join(bits), fontsize=6.0, color="#6b7580", va="center")
    ax.set_yticks(y)
    ax.set_yticklabels([r["display_name"] for r in d], fontsize=6.4)
    for t, r in zip(ax.get_yticklabels(), d):
        t.set_color(RED if (r["regulatory_epigenomic_samples"]
                            + creg.get(r["atlas_disease"], 0)) == 0 else CLS[r["age_class"]])
    ax.tick_params(axis="y", length=0, pad=3)
    ax.invert_yaxis(); ax.set_xscale("log"); ax.set_xlim(1, 30000)
    ax.set_xticks([1,10,100,1000,10000]); ax.set_xticklabels(["1","10","100","1,000","10,000"])
    ax.set_xlabel("Regulatory epigenomic samples (log scale)", fontsize=8)
    ax.spines["left"].set_visible(False)
    ax.legend(handles=[Patch(facecolor=OPEN, label="Open — downloadable and reusable today"),
                       Patch(facecolor=LOCK, label="Controlled — exists in EGA, requires a "
                                                   "per-dataset data-access agreement"),
                       Patch(facecolor=RED, label="Absent — no record in any archive")],
              fontsize=6.6, frameon=False, loc="upper left", bbox_to_anchor=(0.0, -0.075),
              handlelength=1.0, handleheight=0.75)
    foot = ("A further " + f"{sum(unres.values()):,}" + " regulatory epigenomes sit in EGA "
            "under controlled access in cohorts whose subtype EGA does not state — "
            + "; ".join(f"{k}: {v:,}" for k, v in unres.most_common())
            + ".\nThe largest is the St Jude CSTN rhabdomyosarcoma ChIP-seq "
              "(EGAD00001004312 + EGAD00001006398), none of which is in GEO. "
              "They are excluded from the bars rather than attributed by guesswork.")
    fig.text(0.0, -0.135, foot, fontsize=5.9, color="#6b7580", ha="left", va="top")
    ax.set_title("Some of the gap is unmeasured biology; some of it is data that exists "
                 "and cannot be reused\n1,207 regulatory epigenomes sit behind data-access "
                 "committees — but not one covers the nine entities with nothing at all",
                 fontsize=8.8, weight="bold", loc="left", pad=8)
    tidy(ax)
    save(fig, "F24_open_vs_controlled_vs_absent")

f24()
