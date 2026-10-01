#!/usr/bin/env python3
"""F25 -- where sarcoma epigenomic data lives, and what it takes to get it.

Counts every epigenomic sample the atlas knows about, by source archive, split by whether
it can be downloaded today or needs a data-access agreement. A third state emerged during
the harvest and is shown because it is real: data that exists, is described in print or in
a portal, and is deposited nowhere at all.

Unit note, stated on the figure: GEO, EBI and EGA are counted in samples; CCDI is counted
in participants, because that is the level at which its metadata resolves. They are not
perfectly interchangeable and the figure says so.
"""
import csv, gzip, os, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Patch
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
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.04,
})

OPEN, WALL, NONE = "#2a78d6", "#eda100", "#d03b3b"   # NONE now = "none exists"
DARK, MUTED, GRID = "#2b3440", "#898781", "#e1e0d9"
K = FuncFormatter(lambda v, p: f"{int(v):,}")

EPI = {"ChIP-seq","CUT&RUN","CUT&Tag","ChIP-exo","ChIP-chip","ATAC-seq","scATAC-seq",
       "DNase-seq","FAIRE-seq","MNase-seq","Hi-C","HiChIP","Micro-C","Capture-HiC",
       "ChIA-PET","4C-seq","Repli-seq","WGBS","RRBS","Methyl-array","MeDIP/hMeDIP",
       "Bisulfite-PCR"}
REG = {"ChIP-seq","CUT&RUN","CUT&Tag","ChIP-exo","ChIP-chip","ATAC-seq","scATAC-seq",
       "DNase-seq","FAIRE-seq","MNase-seq","Hi-C","HiChIP","Micro-C","Capture-HiC",
       "ChIA-PET","4C-seq"}
DNA = {"WGBS","RRBS","Methyl-array","MeDIP/hMeDIP","Bisulfite-PCR"}

# EBI studies that exist only at EBI, sized from the ENA and BioStudies APIs.
# E-MTAB-9875 (1,974) is held out: it is the sarcoma methylation classifier, essentially
# the same material as GSE140686, and adding it would double-count GEO.
EBI_REG, EBI_DNA = 70, 713
EBI_HELD_OUT = 1974


def main():
    t4 = list(csv.DictReader(topen("T4_samples_atomic.tsv"), delimiter="\t"))
    e = [r for r in t4 if r["assay_class"] in EPI and r["is_duplicate"] != "Y"
         and r["entity_kind"] == "sarcoma"]
    isreg = lambda r: (r["assay_class"] in REG
                       and r["epi_target_norm"] not in ("input/none", "none"))
    geo_reg = sum(1 for r in e if isreg(r))
    geo_dna = sum(1 for r in e if r["assay_class"] in DNA)

    t15 = list(csv.DictReader(topen("T15_controlled_access.tsv"), delimiter="\t"))
    def ctl(src, fam):
        return sum(int(r["n_samples"] or 0) for r in t15
                   if r["source"] == src and r["assay_family"] == fam)
    ega_reg, ega_dna = ctl("EGA", "regulatory"), ctl("EGA", "DNA methylation")
    ccdi_dna = ctl("CCDI", "DNA methylation")
    sj = sum(int(r["n_samples"] or 0) for r in t15 if r.get("contributor") == "St Jude CSTN")

    # Resources with a formal request route but no archive accession. These are behind
    # the same paperwork wall as EGA and dbGaP, not unobtainable -- so they are counted
    # with the wall, and flagged in the footnote as the weakest case within it.
    noacc = 0
    try:
        noacc = sum(int(r["n_samples"] or 0) for r in
                    csv.DictReader(topen("T27_access_routes.tsv"), delimiter="\t")
                    if r.get("has_accession", "").startswith("no"))
    except Exception:
        pass
    # sample x assay manifest recovered from the CSTN ProteinPaint browsers
    viz_reg = viz_samples = 0
    try:
        tk = list(csv.DictReader(topen("T28_stjude_viz_tracks.tsv"), delimiter="\t"))
        live = [r for r in tk if not r.get("availability", "").startswith("withdrawn")]
        viz_reg = sum(1 for r in live if r["assay_family"] == "regulatory"
                      and r["is_control"] != "Y")
        viz_samples = len({r["model_id"] for r in live})
    except Exception:
        pass

    # source, tier, regulatory, methylation, footnote marker
    S = [
        ("NCBI GEO",                       "open", geo_reg, geo_dna, ""),
        ("EBI ArrayExpress / ENA\n(studies not in GEO)", "open", EBI_REG, EBI_DNA, "a"),
        ("EBI EpiRR / IHEC\nreference epigenomes", "open", 0, 0, ""),
        ("EBI EGA",                        "wall", ega_reg, ega_dna, "b"),
        ("NCI CCDI\n(dbGaP)",              "wall", 0, ccdi_dna, "c"),
        ("NCI GDC\n(TARGET / TCGA)",       "wall", 0, 0, "d"),
        ("St Jude CSTN\n(no accession)",   "wall", 0, noacc, "e"),
    ]
    COL = {"open": OPEN, "wall": WALL, "none": NONE}

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.9),
                             gridspec_kw={"width_ratios": [1, 1], "wspace": 0.08})
    y = list(range(len(S)))[::-1]

    for ax, idx, title, sub in (
            (axes[0], 2, "Regulatory epigenomics",
             "ChIP-seq · CUT&RUN · ATAC · DNase · Hi-C · HiChIP"),
            (axes[1], 3, "DNA methylation",
             "arrays · WGBS · RRBS · MBD-seq")):
        vals = [s[idx] for s in S]
        ax.barh(y, vals, color=[COL[s[1]] for s in S], height=0.62, edgecolor="none")
        for yy, s in zip(y, S):
            v = s[idx]
            if v:
                ax.text(v * 1.09, yy, f"{v:,}", va="center", fontsize=7.2, color=DARK)
            else:
                ax.scatter(1.15, yy, marker="x", s=20, color=NONE, linewidth=1.2, zorder=3)
                ax.text(1.6, yy, "none", va="center", fontsize=7, color=NONE, weight="bold")
        ax.set_xscale("log"); ax.set_xlim(1, 40000)
        ax.set_xticks([1, 10, 100, 1000, 10000])
        ax.set_xticklabels(["1", "10", "100", "1,000", "10,000"])
        ax.set_xlabel("samples (log scale)", fontsize=8)
        ax.set_title(title, fontsize=9, weight="bold", loc="left", pad=17)
        ax.text(0, 1.012, sub, transform=ax.transAxes, fontsize=6.6, color=MUTED)
        ax.grid(True, axis="x", color=GRID, linewidth=0.6); ax.set_axisbelow(True)
        ax.spines["left"].set_visible(False)
        ax.tick_params(length=3, pad=2)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([s[0] + (f"  ({s[4]})" if s[4] else "") for s in S], fontsize=7.2)
    axes[0].tick_params(axis="y", length=0, pad=3)
    axes[1].set_yticks(y); axes[1].set_yticklabels([])

    tot_reg_open = geo_reg + EBI_REG
    tot_reg_wall = ega_reg
    tot_dna_open = geo_dna + EBI_DNA
    tot_dna_wall = ega_dna + ccdi_dna + noacc

    fig.legend(handles=[
        Patch(facecolor=OPEN, label="Freely available — download and reanalyze today"),
        Patch(facecolor=WALL, label="Behind a formal request wall — a data-access "
                                    "agreement per study, signed by the PI and an institution"),
        Patch(facecolor=NONE, label="None exists")],
        fontsize=6.9, frameon=False, loc="upper left", bbox_to_anchor=(0.005, -0.005),
        ncol=1, handlelength=1.1, handleheight=0.8)

    fig.suptitle("Regulatory epigenomics for sarcoma is mostly open. DNA methylation is "
                 "mostly not.", fontsize=9.6, weight="bold", x=0.0, y=1.13, ha="left")
    fig.text(0.0, 1.075,
             f"{100*tot_reg_open/(tot_reg_open+tot_reg_wall):.0f}% of the "
             f"{tot_reg_open+tot_reg_wall:,} regulatory samples can be downloaded today, "
             f"against only {100*tot_dna_open/(tot_dna_open+tot_dna_wall):.0f}% of "
             f"the {tot_dna_open+tot_dna_wall:,} methylation samples.",
             fontsize=7.6, color=DARK, ha="left")

    fig.text(0.0, -0.20,
        "Unit: samples, except CCDI, which resolves to participants. Duplicates removed; "
        "ChIP input and IgG controls excluded from the regulatory counts.\n"
        f"(a) The 23 studies that exist only at EBI, sized from the ENA and BioStudies "
        f"APIs. E-MTAB-9875 ({EBI_HELD_OUT:,} assays) is held out — it is the sarcoma "
        f"methylation classifier, the same material as GSE140686 in GEO.\n"
        f"(b) Includes the St Jude CSTN deposit ({sj:,} samples), which is counted here "
        f"rather than separately because EGA is where it lives.  (c) CCDI has no ChIP-seq "
        f"in its assay vocabulary at all, for any disease.\n"
        "(d) The GDC has no ChIP-seq in its vocabulary either, and its 410 ATAC-seq files "
        "belong to 23 TCGA cohorts, none of them SARC. Methylation-array files exist for "
        "TARGET and TCGA sarcoma projects but are\n      not sample-resolved here, so they "
        "are not counted.\n"
        f"(e) St Jude's COMET methylation project and the St Jude Cloud CSTN dataset. A "
        f"request route exists for both, so they sit behind the same wall as EGA and "
        f"dbGaP — but neither has an archive\n      accession, which makes them the "
        f"weakest case within it.\n"
        f"Also recovered this round: the CSTN ProteinPaint browsers resolve to "
        f"{viz_reg} regulatory tracks on {viz_samples} models (T28, F26), with the "
        f"per-model diagnosis supplied by the St Jude data administrator. That material "
        f"is already\n      counted once inside the EGA bar, so it is not added again — "
        f"what it adds is assay-level resolution EGA does not publish.",
        fontsize=5.9, color=MUTED, ha="left", va="top", linespacing=1.55)

    os.makedirs(FIGURES, exist_ok=True)
    p = os.path.join(FIGURES, "F25_where_the_data_lives.pdf")
    fig.savefig(p); plt.close(fig)
    print("wrote", p)
    print(f"\n  regulatory : {tot_reg_open:,} open / {tot_reg_wall:,} behind a wall "
          f"({100*tot_reg_open/(tot_reg_open+tot_reg_wall):.1f}% free)")
    print(f"  methylation: {tot_dna_open:,} open / {tot_dna_wall:,} behind a formal "
          f"request ({100*tot_dna_open/(tot_dna_open+tot_dna_wall):.1f}% free), of which "
          f"{noacc:,} have no archive accession")
    print(f"  CSTN browsers recovered: {viz_reg} regulatory tracks on {viz_samples} "
          f"models (not added -- already inside the EGA bar)")


if __name__ == "__main__":
    main()
