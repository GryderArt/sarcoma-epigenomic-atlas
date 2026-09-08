#!/usr/bin/env python3
"""The four figures for the white-paper section "Critical assessment of existing datasets".

Composed from the atlas tables, not from the exploratory figures, so every number on
every panel is computed at render time from T3/T4/T13/T15/T27/T28. Four panels maximum
per figure. Vector PDF, Arial, text left editable (pdf.fonttype 42).

WP1  burden versus data           A burden  B data  C per-case  D per-entity scatter
WP2  what the atlas is made of    A provenance  B zero-coverage by assay  C per-entity
WP3  what it takes to get it      A regulatory by source  B methylation by source  C entities
WP4  the walled resource          the St Jude CSTN panel, model by model
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
LINE, PDER, TUMR = "#c3d5e8", "#2a78d6", "#0d366b"      # sequential: distance from patient
DARK, MUTED, GRID = "#2b3440", "#7c8894", "#e3e2dc"
CLS = {"paediatric": PED, "both": BOTH, "adult": ADULT}
LAB = {"paediatric": "Paediatric-predominant", "both": "Both ages", "adult": "Adult-predominant"}
LOGX = ["1", "10", "100", "1,000", "10,000"]

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
def wp1(d):
    """Burden versus data, by age class and by entity."""
    T13 = d["T13"]
    agg = collections.defaultdict(collections.Counter)
    for r in T13:
        a = agg[r["age_class"]]
        a["cases"] += F(r, "US_cases_per_year_all_ages")
        a["epi"] += I(r, "epigenomic_samples")
        a["reg"] += I(r, "regulatory_epigenomic_samples")
    order = ["paediatric", "both", "adult"]

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
        ax.set_yticklabels(["Paediatric", "Both ages", "Adult"] if ax is axA else [],
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
    axD.axhspan(0.42, 0.72, color=NONE, alpha=0.09, zorder=0)
    for a, pts in zeros.items():
        for x, _ in pts:
            axD.scatter(x, 0.56, marker="v", s=26, color=NONE, zorder=3)
    # Label the four best-covered entities and the four highest-burden ones -- the two
    # ends of the argument. Offsets are placed greedily so labels never sit on each other.
    pts = [(F(r, "US_cases_per_year_all_ages"), I(r, "regulatory_epigenomic_samples"),
            r["display_name"].split(" (")[0]) for r in T13
           if F(r, "US_cases_per_year_all_ages") and I(r, "regulatory_epigenomic_samples")]
    pick = {p[2]: p for p in sorted(pts, key=lambda t: -t[1])[:4]}
    pick.update({p[2]: p for p in sorted(pts, key=lambda t: -t[0])[:4]})
    placed = []
    for x, yv, n in sorted(pick.values(), key=lambda t: -t[1]):
        for dxo, dyo, ha in ((0, 9, "center"), (0, -13, "center"),
                             (7, 3, "left"), (-7, 3, "right"),
                             (7, -8, "left"), (-7, -8, "right")):
            px, py = axD.transData.transform((x, yv))
            box = (px + dxo - 26 * (ha != "left"), py + dyo, 52, 10)
            if not any(abs(box[0] - b[0]) < 52 and abs(box[1] - b[1]) < 11 for b in placed):
                placed.append(box)
                axD.annotate(n, (x, yv), fontsize=6, color=DARK, textcoords="offset points",
                             xytext=(dxo, dyo), ha=ha)
                break
    axD.set_xscale("log"); axD.set_yscale("symlog", linthresh=1)
    axD.set_xlim(8, 9000); axD.set_ylim(0.3, 20000)
    axD.set_xticks([10, 100, 1000]); axD.set_xticklabels(["10", "100", "1,000"])
    axD.set_yticks([1, 10, 100, 1000]); axD.set_yticklabels(LOGX[:4])
    axD.set_xlabel("US cases per year (all ages)", fontsize=7.5)
    axD.set_ylabel("regulatory epigenomic samples\nin public archives", fontsize=7.5)
    nz = sum(len(v) for v in zeros.values())
    axD.text(9.2, 0.86, f"{nz} entities with a published incidence rate and no regulatory "
                        f"epigenomics at all — plotted on the floor",
             fontsize=6.4, color=NONE, va="bottom", weight="bold")
    axD.legend(handles=[Patch(facecolor=CLS[g], label=LAB[g]) for g in order],
               fontsize=6.5, frameon=False, loc="lower right", handlelength=1.0,
               handleheight=0.75)
    tidy(axD, grid="both"); panel_tag(axD, "D", dx=-0.075, dy=1.02)

    fig.suptitle("Sarcoma epigenomics has been generated in inverse proportion to who "
                 "gets the disease", fontsize=10, weight="bold", x=0.0, y=1.005, ha="left")
    foot(fig, 0.0, -0.055,
         "Burden is US cases per year summed over the 45 named entities; topography-coded "
         "registry data excludes sarcomas coded to the organ they arise in, so the adult "
         "total is if anything an undercount.\n\n"
         "Regulatory = ChIP-seq, CUT&RUN, CUT&Tag, ChIP-exo/chip, ATAC, scATAC, DNase, "
         "FAIRE, MNase, Hi-C, HiChIP, Micro-C, Capture-HiC, ChIA-PET and 4C, with input "
         "and IgG controls excluded.\n\n"
         "Four further entities have no population rate published anywhere and are omitted "
         "from panel D only; they are counted everywhere else.")
    save(fig, "WP1_burden_versus_data")
    return {"ratio": ratio, "agg": agg, "n_scatter_zero": nz}


# =====================================================================================
def wp2(d):
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
                        Patch(facecolor=TUMR, label="Patient tumour, metastasis, recurrence")],
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
                   "(tumour, metastasis, recurrence, PDX or organoid)", fontsize=7.5)
    axC.spines["left"].set_visible(False)
    tidy(axC); panel_tag(axC, "C", dx=-0.42, dy=1.008)

    cl, pdx, pt = split(reg)
    fig.suptitle(f"{len(none_)} of {len(T13)} entities have never had a patient's tumour "
                 f"profiled for active chromatin, accessibility or 3D architecture",
                 fontsize=10, weight="bold", x=0.0, y=0.995, ha="left")
    n_prov = cl + pdx + pt
    foot(fig, 0.0, -0.012,
         f"Panel A denominator is the {n_prov:,} de-duplicated sarcoma regulatory samples "
         f"whose material of origin is stated: {100*cl/n_prov:.0f}% are established cell "
         f"lines and {100*pt/n_prov:.0f}% are patient tumour. Samples recorded only as "
         f"normal/reference, mouse model or unspecified are excluded from that denominator "
         f"rather than assigned by guesswork; against the full de-duplicated set the "
         f"cell-line share is {100*cl/len(reg):.0f}%.\n\n"
         f"Panel C separates patient-derived material from cell lines deliberately — a "
         f"synovial sarcoma organoid grown from a patient is not a decades-old line, and "
         f"collapsing the two would call several entities a hard zero that are not.\n\n"
         f"Absence here means no public, entity-labelled deposit. It is not proof that no "
         f"experiment was done: see Figure WP3.")
    save(fig, "WP2_what_the_atlas_is_made_of")
    return {"n_zero_pd": len(none_), "pct_cell_line": 100*cl/(cl+pdx+pt),
            "pct_tumour": 100*pt/(cl+pdx+pt), "zero_by_assay": vals}


# =====================================================================================
def wp3(d):
    """Where the data lives, and what it takes to get it."""
    T13, T15, reg, epi = d["T13"], d["T15"], d["reg"], d["epi"]
    DNA = {"WGBS", "RRBS", "Methyl-array", "MeDIP/hMeDIP", "Bisulfite-PCR"}
    geo_reg = len(reg)
    geo_dna = sum(1 for r in epi if r["assay_class"] in DNA)

    def ctl(src, fam):
        return sum(I(r, "n_samples") for r in T15
                   if r["source"] == src and r["assay_family"] == fam)
    ega_reg, ega_dna = ctl("EGA", "regulatory"), ctl("EGA", "DNA methylation")
    ccdi_dna = ctl("CCDI", "DNA methylation")
    sj = sum(I(r, "n_samples") for r in T15 if r.get("contributor") == "St Jude CSTN")
    noacc = sum(I(r, "n_samples") for r in
                csv.DictReader(topen("T27_access_routes.tsv"), delimiter="\t")
                if r.get("has_accession", "").startswith("no"))
    EBI_REG, EBI_DNA, EBI_HELD = 70, 713, 1974

    S = [("NCBI GEO", "open", geo_reg, geo_dna),
         ("EBI ArrayExpress / ENA\n(studies not in GEO)", "open", EBI_REG, EBI_DNA),
         ("EBI EpiRR / IHEC\nreference epigenomes", "open", 0, 0),
         ("EBI EGA", "wall", ega_reg, ega_dna),
         ("NCI CCDI (dbGaP)", "wall", 0, ccdi_dna),
         ("NCI GDC (TARGET / TCGA)", "wall", 0, 0),
         ("St Jude CSTN\n(no accession)", "wall", 0, noacc)]
    COL = {"open": OPEN, "wall": WALL}

    fig = plt.figure(figsize=(7.4, 5.3))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.7, 0.85], hspace=0.72, wspace=0.06)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, :])

    y = list(range(len(S)))[::-1]
    for ax, idx, ttl, sub, tag in (
            (axA, 2, "Regulatory epigenomics",
             "ChIP-seq · CUT&RUN · CUT&Tag · ATAC · DNase · Hi-C · HiChIP", "A"),
            (axB, 3, "DNA methylation", "arrays · WGBS · RRBS · MeDIP", "B")):
        ax.barh(y, [t[idx] for t in S], color=[COL[t[1]] for t in S], height=0.6,
                edgecolor="none")
        for yy, t in zip(y, S):
            if t[idx]:
                ax.text(t[idx] * 1.10, yy, f"{t[idx]:,}", va="center", fontsize=6.8,
                        color=DARK)
            else:
                ax.scatter(1.15, yy, marker="x", s=18, color=NONE, linewidth=1.1, zorder=3)
                ax.text(1.55, yy, "none", va="center", fontsize=6.5, color=NONE,
                        weight="bold")
        ax.set_xscale("log"); ax.set_xlim(1, 40000)
        ax.set_xticks([1, 10, 100, 1000, 10000]); ax.set_xticklabels(LOGX)
        ax.set_xlabel("samples (log scale)", fontsize=7.5)
        ax.set_title(ttl, fontsize=8.6, weight="bold", loc="left", pad=15)
        ax.text(0, 1.015, sub, transform=ax.transAxes, fontsize=6.2, color=MUTED)
        ax.spines["left"].set_visible(False)
        tidy(ax)
        panel_tag(ax, tag, dx=-0.44 if ax is axA else -0.06, dy=1.13)
    axA.set_yticks(y); axA.set_yticklabels([t[0] for t in S], fontsize=6.8)
    axA.tick_params(axis="y", length=0, pad=3)
    axB.set_yticks(y); axB.set_yticklabels([])

    # ---- C: the wall has layers
    rungs = [("Open", geo_reg + EBI_REG + geo_dna + EBI_DNA, OPEN,
              "download and reanalyse today"),
             ("Formal request,\narchive accession", ega_reg + ega_dna + ccdi_dna, WALL,
              "EGA or dbGaP: a data-access agreement per study, signed by a PI "
              "and an institution"),
             ("Formal request,\nno accession", noacc, WALL,
              "a request route exists, but nothing to cite and no committee of record — "
              "the weakest case within the wall")]
    yy = list(range(len(rungs)))[::-1]
    top = max(n for _, n, _, _ in rungs)
    for i2, (lab, n, col, note) in zip(yy, rungs):
        axC.barh(i2, n, color=col, height=0.55, edgecolor="white", linewidth=0.8,
                 hatch="///" if "no accession" in lab else None)
        axC.text(n + top * 0.02, i2, f"{n:,}", va="center", fontsize=7, weight="bold",
                 color=DARK)
        axC.text(top * 1.30, i2, note, va="center", fontsize=6.1, color=MUTED)
    axC.set_yticks(yy); axC.set_yticklabels([r[0] for r in rungs], fontsize=6.8)
    axC.set_xlim(0, (geo_reg + EBI_REG + geo_dna + EBI_DNA) * 3.15)
    axC.set_xticks([0, 5000, 10000, 15000])
    axC.set_xticklabels(["0", "5,000", "10,000", "15,000"])
    axC.set_xlabel("epigenomic samples, regulatory and methylation combined", fontsize=7.5,
                   loc="left")
    axC.set_title("The wall is not one wall", fontsize=8.6, weight="bold", loc="left", pad=6)
    axC.spines["left"].set_visible(False); axC.tick_params(axis="y", length=0)
    tidy(axC); panel_tag(axC, "C", dx=-0.44, dy=1.08)

    treg_o, treg_w = geo_reg + EBI_REG, ega_reg
    tdna_o, tdna_w = geo_dna + EBI_DNA, ega_dna + ccdi_dna + noacc
    fig.legend(handles=[Patch(facecolor=OPEN, label="Freely available"),
                        Patch(facecolor=WALL, label="Behind a formal request"),
                        Patch(facecolor=NONE, label="None exists")],
               fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(0.09, 0.055),
               ncol=3, handlelength=1.0, handleheight=0.75)
    fig.suptitle(f"Regulatory epigenomics for sarcoma is {100*treg_o/(treg_o+treg_w):.0f}% "
                 f"open. DNA methylation is {100*tdna_o/(tdna_o+tdna_w):.0f}%.",
                 fontsize=10, weight="bold", x=0.0, y=1.05, ha="left")
    foot(fig, 0.0, -0.10,
         f"Unit is samples, except CCDI, which resolves to participants — the level at "
         f"which its metadata resolves. Duplicates removed; ChIP input and IgG controls "
         f"excluded from the regulatory counts.\n\n"
         f"The EGA bar includes the St Jude CSTN deposit ({sj:,} samples), counted there "
         f"rather than separately because EGA is where it lives. E-MTAB-9875 "
         f"({EBI_HELD:,} assays) is held out of the EBI bar: it is the sarcoma methylation "
         f"classifier, the same material as GSE140686 in GEO. Neither CCDI nor the GDC has "
         f"ChIP-seq in its assay vocabulary at all, for any disease.\n\n"
         f"Not one of the "
         f"{sum(1 for r in T13 if I(r,'regulatory_epigenomic_samples')==0)} entities with "
         f"no regulatory epigenomics has anything behind the wall either — controlled "
         f"access explains part of the gap, but none of the hard zeros.")
    save(fig, "WP3_where_the_data_lives")
    return {"reg_open_pct": 100*treg_o/(treg_o+treg_w),
            "dna_open_pct": 100*tdna_o/(tdna_o+tdna_w),
            "noacc": noacc, "sj": sj, "reg_wall": treg_w, "dna_wall": tdna_w,
            "reg_total": treg_o+treg_w, "dna_total": tdna_o+tdna_w}


# =====================================================================================
ACTIVE, REPRESS, METH, DERIV = "#eda100", "#b8791f", "#7a5cc4", "#c9c6bd"
COLS = [("H3K27Ac", "H3K27ac", ACTIVE), ("H3K9-14Ac", "H3K9/14ac", ACTIVE),
        ("H3K4me1", "H3K4me1", ACTIVE), ("H3K4me2", "H3K4me2", ACTIVE),
        ("H3K4me3", "H3K4me3", ACTIVE), ("BRD4", "BRD4", ACTIVE),
        ("RNAPolII", "RNA Pol II", ACTIVE), ("H3K36me3", "H3K36me3", ACTIVE),
        ("H3K27me3", "H3K27me3", REPRESS), ("H3K9me3", "H3K9me3", REPRESS),
        ("CTCF", "CTCF", REPRESS), ("WGBS", "WGBS", METH),
        ("SE*", "super-enhancer\ncalls", DERIV),
        ("ATAC-seq", "ATAC-seq", NONE), ("DNase-seq", "DNase-seq", NONE),
        ("CUT&RUN", "CUT&RUN /\nCUT&Tag", NONE), ("Hi-C", "Hi-C / HiChIP", NONE)]
DERIVED_SET = {"SE", "SE.noK4me3", "SuperEnhancer(SE)"}
ORDER = ["FN-RMS", "FP-RMS", "RMS-MYOD1", "Osteosarcoma", "Ewing", "Liposarcoma-NOS",
         "Sarcoma NOS", "normal/reference"]
PRETTY = {"Sarcoma NOS": "Sarcoma NOS (high grade)", "normal/reference": "normal / reference"}


def wp4(d):
    """The richest walled resource, model by model and mark by mark."""
    rows = list(csv.DictReader(topen("T28_stjude_viz_tracks.tsv"), delimiter="\t"))
    have, meta = collections.defaultdict(set), {}
    for r in rows:
        have[r["model_id"]].add("SE*" if r["assay"] in DERIVED_SET else r["assay"])
        meta.setdefault(r["model_id"], r)
    models = sorted(meta, key=lambda m: (ORDER.index(meta[m]["atlas_entity"]), m))
    groups = collections.OrderedDict()
    for m in models:
        groups.setdefault(meta[m]["atlas_entity"], []).append(m)

    nrow, ncol = len(models), len(COLS)
    fig, ax = plt.subplots(figsize=(7.4, 0.185 * nrow + 1.9))
    y, ylab, bands = 0, [], []
    for ent, ms in groups.items():
        bands.append((y, len(ms), ent))
        for m in ms:
            got, wd = have[m], meta[m]["availability"].startswith("withdrawn")
            for x, (key, _, col) in enumerate(COLS):
                if key in got:
                    ax.add_patch(Rectangle((x + .12, y + .14), .76, .72, facecolor=col,
                                           edgecolor="none", alpha=0.42 if wd else 1.0))
                else:
                    ax.add_patch(Rectangle((x + .12, y + .14), .76, .72, facecolor="none",
                                           edgecolor=GRID, lw=.6))
            f_ = meta[m]["fusion_or_driver"]
            ylab.append((y + .5, m + (f"  {f_}" if f_.endswith("::FOXO1") else ""), wd,
                         "NOT in the administrator" in meta[m]["entity_call_basis"]))
            y += 1
    ax.set_xlim(0, ncol); ax.set_ylim(0, nrow); ax.invert_yaxis()
    ax.set_yticks([p for p, *_ in ylab])
    ax.set_yticklabels([t for _, t, *_ in ylab], fontsize=6.1)
    for lab, (_, _, wd, flag) in zip(ax.get_yticklabels(), ylab):
        if wd: lab.set_color(NONE)
        elif flag: lab.set_color(MUTED); lab.set_style("italic")
    ax.set_xticks([i2 + .5 for i2 in range(ncol)])
    ax.set_xticklabels([c[1] for c in COLS], fontsize=6.3, rotation=55, ha="left",
                       rotation_mode="anchor")
    ax.xaxis.set_ticks_position("top"); ax.tick_params(length=0, pad=2)
    for sp in ax.spines.values():
        sp.set_visible(False)
    tf = mtrans.blended_transform_factory(ax.transAxes, ax.transData)
    for y0, n, ent in bands:
        ax.plot([1.012, 1.012], [y0 + .12, y0 + n - .12], color=DARK, lw=1.6, transform=tf,
                clip_on=False, solid_capstyle="butt")
        ax.text(1.026, y0 + n / 2, PRETTY.get(ent, ent), va="center", ha="left",
                fontsize=6.6, color=DARK, weight="bold", transform=tf, clip_on=False)
    ax.axvline(13, color=DARK, lw=.9)

    universal = set.intersection(*have.values())
    n_uni = len(universal - {"WGBS", "SE*", "INPUT"})
    fig.suptitle("The densest sarcoma epigenome resource that exists is complete, "
                 "uniform — and entirely behind a wall.",
                 fontsize=10, weight="bold", x=0.0, y=1.085, ha="left")
    fig.text(0.0, 1.038,
             f"Every one of the {nrow} St Jude CSTN models carries the same {n_uni}-mark "
             f"panel plus WGBS and super-enhancer calls (CTCF only in the 2018 study). "
             f"Right of the rule: what it does not contain at all.",
             fontsize=7.3, color=DARK, ha="left")
    fig.legend(handles=[Patch(facecolor=ACTIVE, label="active / elongation mark, or factor"),
                        Patch(facecolor=REPRESS, label="repressive mark or insulator"),
                        Patch(facecolor=METH, label="DNA methylation"),
                        Patch(facecolor=DERIV, label="derived calls (super-enhancers)"),
                        Patch(facecolor="none", edgecolor=GRID, label="absent")],
               fontsize=6.6, frameon=False, loc="upper left", bbox_to_anchor=(0.075, 0.055),
               ncol=3, handlelength=1.1, handleheight=0.8, columnspacing=1.4)
    foot(fig, 0.075, -0.02,
         "Access: all of it requires a St Jude / EGA data access agreement. None of it is "
         "added to any atlas total — the same material is already counted once inside the "
         "EGA bar of Figure WP3. What the enumeration adds is this resolution, which EGA "
         "does not publish.\n\n"
         "Entity calls are the St Jude data administrator's; the fusion partner is from "
         "the CSTN portal model table.\n\n"
         "Red label: SJOS010930_X1, reported by the administrator as withdrawn from CSTN "
         "and no longer available — its tracks are still drawn on the browser page and are "
         "shown faded. Grey italic label: SJRHB010463_X16, which is on the browser page "
         "and in the portal as PAX3::FOXO1 but is absent from the administrator's list; "
         "the discrepancy is left standing.", width=138, size=5.8)
    save(fig, "WP4_the_walled_resource")
    absent = [c[1].replace("\n", " ") for c in COLS if not any(c[0] in have[m] for m in models)]
    return {"n_models": nrow, "n_marks": n_uni, "absent": absent,
            "n_tracks": len(rows),
            "n_reg_tracks": sum(1 for r in rows if r["assay_family"] == "regulatory"
                                and r["is_control"] != "Y")}


def save(fig, stem):
    os.makedirs(FIGURES, exist_ok=True)
    p = os.path.join(FIGURES, stem + ".pdf")
    fig.savefig(p); plt.close(fig)
    print("  wrote", os.path.basename(p))


if __name__ == "__main__":
    d = load()
    s1, s2, s3, s4 = wp1(d), wp2(d), wp3(d), wp4(d)
    print("\nWP1:", {k: v for k, v in s1.items() if k != "agg"})
    print("WP2:", {k: v for k, v in s2.items() if k != "zero_by_assay"})
    print("     zero by assay:", s2["zero_by_assay"])
    print("WP3:", s3)
    print("WP4:", s4)
