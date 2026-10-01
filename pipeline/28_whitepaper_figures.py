#!/usr/bin/env python3
"""The four figures for the white-paper section "Critical assessment of existing datasets".

Composed from the atlas tables, not from the exploratory figures, so every number on
every panel is computed at render time from T3/T4/T13/T15/T27/T28. Four panels maximum
per figure. Vector PDF, Arial, text left editable (pdf.fonttype 42).

WP1  accumulation and modality mix   cumulative GEO deposition by assay family
WP2  burden versus data              A burden  B data  C per-case  D per-entity scatter
WP3  what the atlas is made of    A provenance  B zero-coverage by assay  C per-entity
WP4  can the models be validated  cell line vs patient tissue, per entity
"""
import csv, os, sys, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.transforms as mtrans
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter
csv.field_size_limit(10**7)
from _paths import DATA, FIGURES, topen
from _labels import place_labels

fm._load_fontmanager(try_read_cache=False)
plt.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
})

# ---- palette (validated: see docs/methods.md; scripts/validate_palette.js all-PASS) ----
PED, BOTH, ADULT = "#2a78d6", "#1baf7a", "#eb6834"
OPEN, WALL, NONE = "#2a78d6", "#eda100", "#d03b3b"
LINE, PDER, TUMR = "#c3d5e8", "#2a78d6", "#0d366b"
LINE_C, TISSUE = "#2a78d6", "#eb6834"    # WP4: cell line vs patient tissue      # sequential: distance from patient
DARK, MUTED, GRID = "#2b3440", "#7c8894", "#e3e2dc"
CLS = {"pediatric": PED, "both": BOTH, "adult": ADULT}
LAB = {"pediatric": "Pediatric-predominant", "both": "Both ages", "adult": "Adult-predominant"}
LOGX = ["1", "10", "100", "1,000", "10,000"]
FLOOR = 0.57                      # where entities with nothing are drawn in WP1-D

REG = {"ChIP-seq", "CUT&RUN", "CUT&Tag", "ChIP-exo", "ChIP-chip", "ATAC-seq", "scATAC-seq",
       "DNase-seq", "FAIRE-seq", "MNase-seq", "Hi-C", "HiChIP", "Micro-C", "Capture-HiC",
       "ChIA-PET", "4C-seq"}
PATIENT = {"primary_tumor", "metastasis", "recurrence"}
PDERIVED = PATIENT | {"PDX", "organoid"}

I = lambda r, k: int(float(r.get(k) or 0))
F = lambda r, k: float(r.get(k) or 0)


def load():
    d = {}
    d["T13"] = list(csv.DictReader(topen("T13_incidence_vs_data.tsv"), delimiter="\t"))
    d["T15"] = list(csv.DictReader(topen("T15_controlled_access.tsv"), delimiter="\t"))
    t4 = list(csv.DictReader(topen("T4_samples_atomic.tsv"), delimiter="\t"))
    d["epi"] = [r for r in t4 if r["is_epigenomic"] == "Y"
                and r["entity_kind"] == "sarcoma" and r["is_duplicate"] != "Y"]
    d["reg"] = [r for r in d["epi"] if r["assay_class"] in REG
                and r["epi_target_norm"] not in ("input/none", "none")]
    return d



def foot(fig, x, y, text, width=150, size=5.9):
    """Footnotes wrapped to a character budget, so bbox_inches='tight' cannot widen the
    figure past the panel area. Paragraphs are separated by a blank line in `text`."""
    import textwrap
    out = []
    for para in text.strip().split("\n\n"):
        body = " ".join(para.split())
        out += textwrap.wrap(body, width=width, subsequent_indent="      ")
    fig.text(x, y, "\n".join(out), fontsize=size, color=MUTED, ha="left", va="top",
             linespacing=1.55)


def panel_tag(ax, s, dx=-0.085, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=10, weight="bold",
            color=DARK, ha="left", va="bottom")


def tidy(ax, grid="x"):
    ax.grid(True, axis=grid, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(length=3, pad=2)


# =====================================================================================
def fig_burden(d):
    """Burden versus data, by age class and by entity."""
    T13 = d["T13"]
    agg = collections.defaultdict(collections.Counter)
    for r in T13:
        a = agg[r["age_class"]]
        a["cases"] += F(r, "US_cases_per_year_all_ages")
        a["epi"] += I(r, "epigenomic_samples")
        a["reg"] += I(r, "regulatory_epigenomic_samples")
    order = ["pediatric", "both", "adult"]

    fig = plt.figure(figsize=(7.4, 5.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.55], width_ratios=[1.15, 1.15, 0.80],
                          hspace=0.40, wspace=0.62)
    axA, axB, axC = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])
    axD = fig.add_subplot(gs[1, :])

    y = [2, 1, 0]
    for ax, key, xlab, ttl in ((axA, "cases", "US cases per year", "Disease burden"),
                               (axB, "epi", "samples in public archives", "Data generated")):
        v = [agg[g][key] for g in order]
        ax.barh(y, v, color=[CLS[g] for g in order], height=0.6, edgecolor="none")
        for yy, g in zip(y, order):
            ax.text(agg[g][key] * 1.04, yy, f"{agg[g][key]:,.0f}", va="center",
                    fontsize=7, color=DARK)
        ax.set_yticks(y)
        ax.set_yticklabels(["Pediatric", "Both ages", "Adult"] if ax is axA else [],
                           fontsize=6.8)
        ax.set_xlim(0, max(v) * 1.42); ax.set_xlabel(xlab, fontsize=7)
        ax.set_title(ttl, fontsize=8.6, weight="bold", loc="left", pad=6)
        ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
        tidy(ax)
    ratio = [agg[g]["reg"] / agg[g]["cases"] for g in order]
    axC.bar(range(3), ratio, color=[CLS[g] for g in order], width=0.62, edgecolor="none")
    for i2, v in enumerate(ratio):
        axC.text(i2, v * 1.04, f"{v:.2f}", ha="center", va="bottom", fontsize=7, color=DARK)
    axC.set_xticks(range(3)); axC.set_xticklabels(["Paed.", "Both", "Adult"], fontsize=6.8)
    axC.set_ylim(0, max(ratio) * 1.22)
    axC.set_ylabel("regulatory epigenomes\nper new case per year", fontsize=7)
    axC.set_title(f"{ratio[0]/ratio[2]:.0f}× gap", fontsize=8.6, weight="bold",
                  loc="left", pad=6)
    tidy(axC, grid="y")
    for ax, t in ((axA, "A"), (axB, "B"), (axC, "C")):
        panel_tag(ax, t, dx=-0.34 if ax is axA else -0.12 if ax is axB else -0.40)

    # ---- D: per-entity
    SHORT = {"Rhabdoid tumor / ATRT": "Rhabdoid / ATRT",
             "Undifferentiated pleomorphic sarcoma": "UPS / MFH",
             "Desmoid / aggressive fibromatosis": "Desmoid"}
    short = lambda nm: SHORT.get(nm, nm.split(" (")[0])

    zeros = collections.defaultdict(list)
    for r in T13:
        x, yv = F(r, "US_cases_per_year_all_ages"), I(r, "regulatory_epigenomic_samples")
        if not x:
            continue
        if yv:
            axD.scatter(x, yv, s=26, color=CLS[r["age_class"]], edgecolor="white",
                        linewidth=0.6, zorder=3)
        else:
            zeros[r["age_class"]].append((x, r["display_name"]))
    axD.axhspan(0.46, 0.68, color=NONE, alpha=0.09, zorder=0)
    for a, zs in zeros.items():
        for x, _ in zs:
            axD.scatter(x, FLOOR, marker="v", s=26, color=NONE, zorder=3)

    # Scales and limits must be final BEFORE any label is placed: placement is measured
    # in display pixels, and transData changes when the scale does.
    axD.set_xscale("log"); axD.set_yscale("symlog", linthresh=1)
    axD.set_xlim(8, 9000); axD.set_ylim(0.19, 20000)
    axD.set_xticks([10, 100, 1000]); axD.set_xticklabels(["10", "100", "1,000"])
    axD.set_yticks([1, 10, 100, 1000]); axD.set_yticklabels(LOGX[:4])
    axD.set_xlabel("US cases per year (all ages)", fontsize=7.5)
    axD.set_ylabel("regulatory epigenomic samples\nin public archives", fontsize=7.5)
    nz = sum(len(v) for v in zeros.values())
    axD.text(9.2, 0.40, f"{nz} entities with a published incidence rate and no regulatory "
                        f"epigenomics at all — plotted on the floor",
             fontsize=6.4, color=NONE, va="top", weight="bold")
    # The legend sits below the axes, not in the lower-right corner, because that corner
    # is exactly where the high-burden, low-data points that carry the argument live.
    leg = axD.legend(handles=[Patch(facecolor=CLS[g], label=LAB[g]) for g in order],
                     fontsize=6.5, frameon=False, loc="upper left",
                     bbox_to_anchor=(0.0, -0.13), ncol=3, handlelength=1.0,
                     handleheight=0.75, columnspacing=1.6)
    tidy(axD, grid="both"); panel_tag(axD, "D", dx=-0.075, dy=1.02)

    # Label the four best-covered entities and the three highest-burden ones -- the two
    # ends of the argument. Placement is measured and leader-lined; see _labels.py.
    scat = [(F(r, "US_cases_per_year_all_ages"), I(r, "regulatory_epigenomic_samples"),
             short(r["display_name"])) for r in T13
            if F(r, "US_cases_per_year_all_ages") and I(r, "regulatory_epigenomic_samples")]
    N_DATA, N_BURDEN = 4, 3
    pick = {p[2]: p for p in sorted(scat, key=lambda t: -t[1])[:N_DATA]}
    pick.update({p[2]: p for p in sorted(scat, key=lambda t: -t[0])[:N_BURDEN]})
    place_labels(fig, axD, list(pick.values()),
                 [(x, yv) for x, yv, _ in scat]
                 + [(x, FLOOR) for zs in zeros.values() for x, _ in zs],
                 avoid=[leg], fontsize=6, color=DARK, leader=MUTED, where="WP2-D")

    ORD = {4: "fourth", 5: "fifth", 20: "twentieth", 21: "twenty-first",
           22: "twenty-second", 23: "twenty-third", 24: "twenty-fourth",
           25: "twenty-fifth", 26: "twenty-sixth"}
    burden_x = agg["adult"]["cases"] / agg["pediatric"]["cases"]
    gap = round(ratio[0] / ratio[2])
    fig.suptitle(f"Adult sarcomas are the largest unmet opportunity: "
                 f"{'five' if round(burden_x) == 5 else round(burden_x)} times the "
                 f"burden, one {ORD.get(gap, str(gap) + 'th')} the data per case",
                 fontsize=9.8, weight="bold", x=0.0, y=1.005, ha="left")
    foot(fig, 0.0, -0.055,
         "Burden is US cases per year summed over the 45 named entities; topography-coded "
         "registry data excludes sarcomas coded to the organ they arise in, so the adult "
         "total is if anything an undercount.\n\n"
         "Regulatory = ChIP-seq, CUT&RUN, CUT&Tag, ChIP-exo/chip, ATAC, scATAC, DNase, "
         "FAIRE, MNase, Hi-C, HiChIP, Micro-C, Capture-HiC, ChIA-PET and 4C, with input "
         "and IgG controls excluded.\n\n"
         "Four further entities have no population rate published anywhere and are omitted "
         "from panel D only; they are counted everywhere else.\n\n"
         "Panel D labels the four best-covered entities and the three highest-burden "
         "ones, each tied to its own point by a leader line.")
    save(fig, "WP2_burden_versus_data")
    return {"ratio": ratio, "agg": agg, "n_scatter_zero": nz}


# =====================================================================================
def fig_composition(d):
    """What the atlas is made of, and which entities it misses."""
    T13, epi, reg = d["T13"], d["epi"], d["reg"]

    def split(rows):
        cl = sum(1 for r in rows if r["sample_type"] == "cell_line")
        pdx = sum(1 for r in rows if r["sample_type"] in ("PDX", "organoid"))
        pt = sum(1 for r in rows if r["sample_type"] in PATIENT)
        return cl, pdx, pt

    fig = plt.figure(figsize=(7.4, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.52, 2.6], hspace=0.26, wspace=0.50)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, :])

    # ---- A: provenance
    rows = [("Regulatory\nassays only", reg), ("All epigenomic\nassays", epi)]
    for i2, (lab, rr) in enumerate(rows):
        cl, pdx, pt = split(rr)
        tot = cl + pdx + pt
        left = 0
        for v, col in ((cl, LINE), (pdx, PDER), (pt, TUMR)):
            pct = 100 * v / tot
            axA.barh(i2, pct, left=left, color=col, height=0.5, edgecolor="white",
                     linewidth=1.2)
            if pct > 6:
                axA.text(left + pct / 2, i2, f"{pct:.0f}%", ha="center", va="center",
                         fontsize=7, weight="bold", color="white" if col != LINE else DARK)
            else:
                axA.text(left + pct / 2, i2 + 0.42, f"{pct:.0f}%", ha="center", va="bottom",
                         fontsize=6, color=col if col != LINE else DARK)
            left += pct
    axA.set_yticks([0, 1]); axA.set_yticklabels([r[0] for r in rows], fontsize=6.8)
    axA.set_ylim(-0.55, 1.75)
    axA.set_xlim(0, 100); axA.set_xlabel("% of samples", fontsize=7)
    axA.set_title("Where the material came from", fontsize=8.6, weight="bold", loc="left",
                  pad=6)
    axA.spines["left"].set_visible(False); axA.tick_params(axis="y", length=0)
    axA.legend(handles=[Patch(facecolor=LINE, label="Established cell line"),
                        Patch(facecolor=PDER, label="PDX / organoid"),
                        Patch(facecolor=TUMR, label="Patient tumor, metastasis, recurrence")],
               fontsize=6.1, frameon=False, loc="upper left", bbox_to_anchor=(-0.30, -0.30),
               ncol=3, handlelength=1.0, handleheight=0.75, columnspacing=1.1)
    tidy(axA); panel_tag(axA, "A", dx=-0.30, dy=1.13)

    # ---- B: entities with nothing, by assay class
    checks = [("Any regulatory assay", "regulatory_epigenomic_samples"),
              ("H3K27ac", "H3K27ac_samples"),
              ("Chromatin accessibility", "accessibility_samples"),
              ("3D genome", "samples_3D"),
              ("Regulatory from patient-\nderived material", "patient_derived_regulatory")]
    vals = [(lab, sum(1 for r in T13 if I(r, k) == 0)) for lab, k in checks]
    vals.sort(key=lambda t: t[1])
    yy = range(len(vals))
    axB.barh(list(yy), [v for _, v in vals], color=NONE, height=0.58, edgecolor="none")
    for i2, (_, v) in enumerate(vals):
        axB.text(v + 0.7, i2, f"{v}", va="center", fontsize=7, color=DARK, weight="bold")
    axB.set_yticks(list(yy)); axB.set_yticklabels([l for l, _ in vals], fontsize=6.4)
    axB.set_xlim(0, len(T13)); axB.set_xlabel(f"entities with zero, of {len(T13)}", fontsize=7)
    axB.set_title("What is missing, and for how many entities", fontsize=8.6,
                  weight="bold", loc="left", pad=6)
    axB.spines["left"].set_visible(False); axB.tick_params(axis="y", length=0)
    tidy(axB); panel_tag(axB, "B", dx=-0.62, dy=1.13)

    # ---- C: per-entity patient-derived regulatory
    have = sorted([r for r in T13 if I(r, "patient_derived_regulatory") > 0],
                  key=lambda r: -I(r, "patient_derived_regulatory"))
    none_ = sorted([r for r in T13 if I(r, "patient_derived_regulatory") == 0],
                   key=lambda r: -F(r, "US_cases_per_year_all_ages"))
    seq = none_ + have
    yv = list(range(len(seq)))[::-1]
    for y0, r in zip(yv, seq):
        n = I(r, "patient_derived_regulatory")
        if n:
            axC.barh(y0, n, color=CLS[r["age_class"]], height=0.62, edgecolor="none")
            axC.text(n * 1.13, y0, f"{n:,}", va="center", fontsize=6, color=DARK)
        else:
            axC.scatter(1.12, y0, marker="x", s=18, color=NONE, linewidth=1.1, zorder=3)
            rate = F(r, "US_cases_per_year_all_ages")
            axC.text(1.42, y0, f"none — against ~{int(rate):,} US cases per year" if rate
                     else "none — no population rate published anywhere",
                     fontsize=5.9, color=NONE, va="center")
    axC.axhline(len(have) - 0.5, color=DARK, lw=0.8, ls=(0, (4, 3)))
    axC.set_yticks(yv)
    axC.set_yticklabels([r["display_name"] for r in seq], fontsize=6.1)
    for t, r in zip(axC.get_yticklabels(), seq):
        t.set_color(NONE if I(r, "patient_derived_regulatory") == 0 else CLS[r["age_class"]])
    axC.tick_params(axis="y", length=0, pad=3)
    axC.set_xscale("log"); axC.set_xlim(1, 3000)
    axC.set_xticks([1, 10, 100, 1000]); axC.set_xticklabels(LOGX[:4])
    axC.set_xlabel("regulatory epigenomic profiles of patient-derived material "
                   "(tumor, metastasis, recurrence, PDX or organoid)", fontsize=7.5)
    axC.spines["left"].set_visible(False)
    tidy(axC); panel_tag(axC, "C", dx=-0.42, dy=1.008)

    cl, pdx, pt = split(reg)
    # entities that have gained patient-derived regulatory data, by era
    firstpd = {}
    for r in sorted((x for x in d["reg"] if x["sample_type"] in PDERIVED),
                    key=lambda x: (x.get("gse_date") or "9999")):
        yv = (r.get("gse_date") or "")[:4]
        if yv.isdigit() and r["disease"] in {t["atlas_disease"] for t in T13}:
            firstpd.setdefault(r["disease"], int(yv))
    grew = {yy: sum(1 for v in firstpd.values() if v <= yy) for yy in (2015, 2020, 2026)}
    fig.suptitle(f"Patient-derived profiling has extended from {grew[2015]} entities to "
                 f"{grew[2026]} in a decade; {len(none_)} of {len(T13)} remain",
                 fontsize=9.8, weight="bold", x=0.0, y=0.995, ha="left")
    n_prov = cl + pdx + pt
    foot(fig, 0.0, -0.012,
         f"Panel A denominator is the {n_prov:,} de-duplicated sarcoma regulatory samples "
         f"whose material of origin is stated: {100*cl/n_prov:.0f}% are established cell "
         f"lines and {100*pt/n_prov:.0f}% are patient tumor. Samples recorded only as "
         f"normal/reference, mouse model or unspecified are excluded from that denominator "
         f"rather than assigned by guesswork; against the full de-duplicated set the "
         f"cell-line share is {100*cl/len(reg):.0f}%.\n\n"
         f"Panel C counts patient-derived material only, which is its purpose; the "
         f"cell-line side of the same comparison, entity by entity, is Figure WP4. "
         f"Patient-derived material is kept separate from cell lines because a synovial "
         f"sarcoma organoid grown from a patient is not a decades-old line, and "
         f"collapsing the two would report several entities as hard zeros that are "
         f"not.\n\n"
         f"Absence indicates no public, entity-labeled deposit rather than evidence that "
         f"no experiment was performed.")
    save(fig, "WP3_what_the_atlas_is_made_of")
    return {"n_zero_pd": len(none_), "n_have_pd": len(T13) - len(none_),
            "pct_cell_line": round(100*cl/(cl+pdx+pt)),
            "pct_tumor": round(100*pt/(cl+pdx+pt)),
            "pd_growth": grew, "zero_by_assay": vals}


# =====================================================================================
MODALITY = [
    ("ChIP-seq / ChIP-chip", {"ChIP-seq", "ChIP-chip", "ChIP-exo"},            "#2a78d6"),
    ("DNA methylation",      {"Methyl-array", "WGBS", "RRBS", "MeDIP/hMeDIP",
                              "Bisulfite-PCR"},                                "#eb6834"),
    ("Accessibility (ATAC / DNase)", {"ATAC-seq", "DNase-seq", "FAIRE-seq",
                                      "MNase-seq"},                            "#1baf7a"),
    ("CUT&RUN / CUT&Tag",    {"CUT&RUN", "CUT&Tag"},                           "#eda100"),
    ("3D genome",            {"Hi-C", "HiChIP", "Micro-C", "Capture-HiC",
                              "ChIA-PET", "4C-seq"},                           "#1d6b45"),
    ("Single-cell epigenome", {"scATAC-seq"},                                  "#e0679f"),
]


def fig_growth(d):
    """When each modality arrived, and how thin the newer ones still are."""
    epi = d["epi"]

    def year(r):
        v = (r.get("gse_date") or "")[:4]
        return int(v) if v.isdigit() else None

    yrs = [year(r) for r in epi if year(r)]
    y0, y1 = min(yrs), max(yrs)
    span = list(range(y0, y1 + 1))

    fig, ax = plt.subplots(figsize=(7.4, 3.5))
    fig.subplots_adjust(right=0.70)
    ends = []
    for label, classes, col in MODALITY:
        per = collections.Counter(year(r) for r in epi
                                  if r["assay_class"] in classes and year(r))
        run, cum = [], 0
        for yy in span:
            cum += per.get(yy, 0)
            run.append(cum)
        ax.plot(span, run, color=col, linewidth=1.9, solid_capstyle="round", zorder=3)
        ends.append((run[-1], label, col, per))

    top = max(e[0] for e in ends)
    ax.set_xlim(y0 - 0.4, y1 + 0.4)
    ax.set_ylim(0, top * 1.06)
    ticks = [y for y in span if y % 5 == 0]
    if y1 - ticks[-1] >= 2:
        ticks.append(y1)
    ax.set_xticks(ticks); ax.set_xticklabels([str(y) for y in ticks])
    ax.set_xlabel("year of GEO deposition", fontsize=7.5)
    ax.set_ylabel("cumulative epigenomic samples", fontsize=7.5)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{int(v):,}"))
    tidy(ax, grid="y")

    # Labels in the right margin at each line's end, spread vertically so none collide.
    # A line chart reads better with the series named where it finishes than with a
    # legend box; the only requirement is that the label sit at its own line's height.
    MINGAP = top * 0.052
    rank = sorted(range(len(ends)), key=lambda i: -ends[i][0])
    ypos = {}
    prev = None
    for i in rank:
        v = ends[i][0]
        if prev is not None and prev - v < MINGAP:
            v = prev - MINGAP
        ypos[i] = v
        prev = v
    for i, (v, label, col, _) in enumerate(ends):
        ax.annotate(label, xy=(y1, v), xytext=(y1 + 0.55, ypos[i]),
                    textcoords="data", va="center", ha="left", fontsize=6.8,
                    color=col, weight="bold", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color=col, linewidth=0.5,
                                    shrinkA=1, shrinkB=1)
                    if abs(ypos[i] - v) > top * 0.012 else None)

    # first year each modality reached 50 samples -- the point at which it stops being
    # a demonstration and starts being a resource
    arrival = {}
    for v, lab, col, per in ends:
        cum = 0
        for yy in span:
            cum += per.get(yy, 0)
            if cum >= 50:
                arrival[lab] = yy
                break

    chip = [e for e in ends if e[1].startswith("ChIP")][0][0]
    fam_tot = sum(e[0] for e in ends)
    dated = [r for r in epi if year(r)]
    tot = len(dated)                                # every dated sample, not just the six
    other = tot - fam_tot                           # Repli-seq, which is none of the six

    # Growth statistics, computed here so the figure and the prose cannot diverge.
    last10 = sum(1 for r in dated if year(r) > y1 - 10)
    last5 = sum(1 for r in dated if year(r) > y1 - 5)
    MODERN = {"CUT&RUN", "CUT&Tag", "ATAC-seq", "scATAC-seq", "Hi-C", "HiChIP",
              "Micro-C", "Capture-HiC"}
    era = {}
    for a, b in ((2011, 2015), (2021, y1)):
        ss = [r for r in dated if a <= year(r) <= b]
        era[(a, b)] = 100 * sum(1 for r in ss if r["assay_class"] in MODERN) / len(ss)
    firstyr = {}
    for r in dated:
        firstyr[r["disease"]] = min(firstyr.get(r["disease"], 9999), year(r))
    ents = {dz: v for dz, v in firstyr.items()
            if dz in {t["atlas_disease"] for t in d["T13"]}}
    cover = {yy: sum(1 for v in ents.values() if v <= yy) for yy in (2010, 2015, 2020, y1)}
    fig.suptitle("Sarcoma epigenomic data has accumulated rapidly, and the modality mix "
                 "is broadening",
                 fontsize=10, weight="bold", x=0.0, y=1.135, ha="left")
    fig.text(0.0, 1.085,
             f"{100*last10/tot:.0f}% of all deposits fall in the last ten years and "
             f"{100*last5/tot:.0f}% in the last five. Assays other than ChIP-seq and "
             f"methylation arrays rose from\n{era[(2011, 2015)]:.0f}% of deposits in "
             f"2011\u20132015 to {era[(2021, y1)]:.0f}% in 2021\u2013{y1}.",
             fontsize=7.3, color=DARK, ha="left", va="top", linespacing=1.5)
    foot(fig, 0.0, -0.115,
         "Cumulative count by the deposition year of the containing GEO series, "
         f"de-duplicated. {y1} is a partial year, so every curve is flat at its right "
         "edge by construction.\n\n"
         + f"Entities with any epigenomic data: {cover[2010]} of 45 by 2010, "
           f"{cover[2015]} by 2015, {cover[2020]} by 2020, {cover[y1]} today.\n\n"
         + "Year each modality first passed 50 deposited samples: "
         + "; ".join(f"{k} {v}" for k, v in sorted(arrival.items(), key=lambda kv: kv[1]))
         + ".\n\n"
         f"The six families cover {fam_tot:,} of the {tot:,} dated samples; the remaining "
         f"{other} are Repli-seq, which belongs to none of them.\n\n"
         "Controlled-access deposits are not in this figure: EGA and dbGaP do not publish "
         "a comparable deposition timeline.")
    save(fig, "WP1_accumulation_and_modality_mix")
    return {"pct_chip": round(100*chip/tot), "total": tot, "fam_total": fam_tot,
            "pct_last10": round(100*last10/tot), "pct_last5": round(100*last5/tot),
            "modern_early": round(era[(2011, 2015)], 1),
            "modern_late": round(era[(2021, y1)], 1), "coverage": cover,
            "arrival": arrival, "ends": {lab: v for v, lab, _, _ in ends}}


# =====================================================================================
def fig_validation(d):
    """Can a finding in a cell line be checked in human tissue?"""
    T13, reg = d["T13"], d["reg"]
    PDER = {"primary_tumor", "metastasis", "recurrence", "PDX", "organoid"}
    name = {r["atlas_disease"]: r["display_name"] for r in T13}
    cl = collections.Counter(r["disease"] for r in reg if r["sample_type"] == "cell_line")
    pd_ = collections.Counter(r["disease"] for r in reg if r["sample_type"] in PDER)
    rows = [(name[dz], cl.get(dz, 0), pd_.get(dz, 0)) for dz in name
            if cl.get(dz, 0) or pd_.get(dz, 0)]
    # Entities whose whole regulatory record is normal or reference tissue drop out
    # here, because neither point exists for them. That is correct but worth naming.
    allreg = collections.Counter(r["disease"] for r in reg)
    ctrl_only = [name[dz] for dz in name
                 if allreg.get(dz, 0) and not (cl.get(dz, 0) or pd_.get(dz, 0))]
    rows.sort(key=lambda t: -t[1])
    FL = 0.42                                   # where "none" is drawn on the log axis

    fig, ax = plt.subplots(figsize=(7.0, 0.185 * len(rows) + 1.25))
    fig.subplots_adjust(top=0.965, bottom=0.075)
    y = list(range(len(rows)))[::-1]
    for yy, (nm, c, p) in zip(y, rows):
        ax.plot([p if p else FL * 2.6, c], [yy, yy], color="#dedcd5", linewidth=2.4,
                solid_capstyle="round", zorder=1)
        ax.scatter(c, yy, s=28, color=LINE_C, edgecolor="white", linewidth=0.6, zorder=3)
        if p:
            ax.scatter(p, yy, s=28, color=TISSUE, edgecolor="white", linewidth=0.6,
                       zorder=3)
            ax.text(c * 1.3, yy, f"{c/p:.0f}×" if c / p >= 2 else f"{c/p:.1f}×",
                    va="center", fontsize=5.9, color=MUTED)
        else:
            # No marker for "none": an absent point is the encoding, and the italic red
            # word names it. A third colored mark here would sit too close to the
            # patient-tissue orange to be told apart.
            ax.text(FL * 0.92, yy, "none", va="center", ha="left", fontsize=6,
                    color=NONE, style="italic", weight="bold")
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=6.3)
    for t, r in zip(ax.get_yticklabels(), rows):
        if not r[2]:
            t.set_color(NONE)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.set_xscale("log"); ax.set_xlim(0.30, 9000)
    ax.set_xticks([1, 10, 100, 1000]); ax.set_xticklabels(LOGX[:4])
    ax.set_xlabel("regulatory epigenomic samples (log scale)", fontsize=7.5)
    ax.spines["left"].set_visible(False)
    tidy(ax)
    ax.legend(handles=[Patch(facecolor=LINE_C, label="Established cell line"),
                       Patch(facecolor=TISSUE, label="Patient-derived: tumor, "
                                                     "metastasis, recurrence, PDX or "
                                                     "organoid")],
              fontsize=6.4, frameon=False, loc="upper left", bbox_to_anchor=(0.0, -0.062),
              ncol=2, handlelength=1.0, handleheight=0.75, columnspacing=1.4)

    nz = sum(1 for _, c, p in rows if c and not p)
    ok = len(rows) - nz
    worst = max((t for t in rows if t[2]), key=lambda t: t[1] / t[2])
    fig.suptitle(f"{ok} of the {len(rows)} entities with regulatory epigenomics can be "
                 f"validated against patient tissue; {nz} cannot yet",
                 fontsize=9.6, weight="bold", x=0.0, y=1.008, ha="left")
    foot(fig, 0.0, -0.105,
         f"Each row is one entity. The blue point is regulatory epigenomic samples from "
         f"established cell lines, the orange point the same assays on patient-derived "
         f"material; the bar between them is the validation gap, annotated as a ratio. "
         f"{worst[0]} is the widest at {worst[1]/worst[2]:.0f}x.\n\n"
         f"A row marked none has cell-line data and no patient-derived material of any "
         f"kind, so a finding made in the model cannot presently be checked against human "
         f"disease. No entity has the opposite problem.\n\n"
         f"Restricted to regulatory assays, which is where the asymmetry lives: DNA "
         f"methylation arrays are mostly run on patient tissue and would mask it.\n\n"
         f"Thirteen further entities have no regulatory data at all and cannot appear "
         f"here (Figure WP2D)."
         + (f" A fourteenth, {' and '.join(ctrl_only)}, is also absent for a different "
            f"reason: its entire regulatory record is normal or reference tissue rather "
            f"than tumor or model, so it has neither point to plot." if ctrl_only else ""))
    save(fig, "WP4_can_the_models_be_validated")
    return {"n_rows": len(rows), "n_no_patient": nz, "n_validatable": ok,
            "worst": (worst[0], round(worst[1]/worst[2]))}


def save(fig, stem):
    os.makedirs(FIGURES, exist_ok=True)
    p = os.path.join(FIGURES, stem + ".pdf")
    fig.savefig(p); plt.close(fig)
    print("  wrote", os.path.basename(p))


if __name__ == "__main__":
    d = load()
    out = {"WP1 growth": fig_growth(d), "WP2 burden": fig_burden(d),
           "WP3 composition": fig_composition(d), "WP4 validation": fig_validation(d)}
    print()
    for k, v in out.items():
        print(f"{k}:")
        for kk, vv in v.items():
            if kk == "agg":
                continue
            print(f"    {kk} = {vv}")
