#!/usr/bin/env python3
"""F26 -- what the richest walled sarcoma epigenome resource does and does not contain.

The St Jude CSTN browsers carry a uniform assay panel across 28 models. Enumerating them
(stage 15) makes two things visible that a sample count cannot:

  1. the panel is complete and identical on every model -- eleven regulatory marks plus
     WGBS, with no missing cells. Nothing in the open literature approaches this for
     sarcoma.
  2. the panel stops where the field has moved. There is no chromatin accessibility, no
     3D genome, and no CUT&RUN or CUT&Tag anywhere in it, and CTCF exists on only the
     nine models in the 2018 osteosarcoma/rare-tumour study.

Entity labels are the St Jude data administrator's, supplied by email; the fusion partner
comes from the CSTN portal's own model table. Both are recorded per row in T28.
"""
import csv, os, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Patch, Rectangle
import matplotlib.transforms as mtrans
from _paths import DATA, FIGURES, topen

fm._load_fontmanager(try_read_cache=False)
plt.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.linewidth": 0.7, "figure.dpi": 200,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.04,
})

WALL, NONE = "#eda100", "#d03b3b"
DARK, MUTED, GRID = "#2b3440", "#898781", "#e1e0d9"
ACTIVE, REPRESS, METH, DERIV = "#eda100", "#b8791f", "#7a5cc4", "#c9c6bd"

# columns, in biological order. Assays the resource does NOT contain are listed too --
# an empty column is the point of the figure, not an omission from it.
COLS = [
    ("H3K27Ac",  "H3K27ac",  ACTIVE), ("H3K9-14Ac", "H3K9/14ac", ACTIVE),
    ("H3K4me1",  "H3K4me1",  ACTIVE), ("H3K4me2",   "H3K4me2",   ACTIVE),
    ("H3K4me3",  "H3K4me3",  ACTIVE), ("BRD4",      "BRD4",      ACTIVE),
    ("RNAPolII", "RNA Pol II", ACTIVE),
    ("H3K36me3", "H3K36me3", ACTIVE),
    ("H3K27me3", "H3K27me3", REPRESS), ("H3K9me3",  "H3K9me3",   REPRESS),
    ("CTCF",     "CTCF",     REPRESS),
    ("WGBS",     "WGBS",     METH),
    ("SE*",      "super-enhancer\ncalls", DERIV),
    ("ATAC-seq", "ATAC-seq", NONE), ("DNase-seq", "DNase-seq", NONE),
    ("CUT&RUN",  "CUT&RUN /\nCUT&Tag", NONE),
    ("Hi-C",     "Hi-C / HiChIP", NONE),
]
DERIVED_SET = {"SE", "SE.noK4me3", "SuperEnhancer(SE)"}
ORDER = ["FN-RMS", "FP-RMS", "RMS-MYOD1", "Osteosarcoma", "Ewing", "Liposarcoma-NOS",
         "Sarcoma NOS", "normal/reference"]
PRETTY = {"Sarcoma NOS": "Sarcoma NOS\n(high grade)", "normal/reference": "normal /\nreference"}


def main():
    rows = list(csv.DictReader(topen("T28_stjude_viz_tracks.tsv"), delimiter="\t"))
    have = collections.defaultdict(set)
    meta = {}
    for r in rows:
        m = r["model_id"]
        have[m].add("SE*" if r["assay"] in DERIVED_SET else r["assay"])
        meta.setdefault(m, r)

    models = sorted(meta, key=lambda m: (ORDER.index(meta[m]["atlas_entity"]), m))
    groups = collections.OrderedDict()
    for m in models:
        groups.setdefault(meta[m]["atlas_entity"], []).append(m)

    nrow, ncol = len(models), len(COLS)
    fig, ax = plt.subplots(figsize=(7.6, 0.185 * nrow + 1.9))

    y = 0; ylab = []; bands = []
    for ent, ms in groups.items():
        bands.append((y, len(ms), ent))
        for m in ms:
            got = have[m]
            withdrawn = meta[m]["availability"].startswith("withdrawn")
            for x, (key, _, col) in enumerate(COLS):
                if key in got:
                    ax.add_patch(Rectangle((x + .12, y + .14), .76, .72,
                                           facecolor=col, edgecolor="none",
                                           alpha=0.45 if withdrawn else 1.0))
                else:
                    ax.add_patch(Rectangle((x + .12, y + .14), .76, .72,
                                           facecolor="none", edgecolor=GRID, lw=.6))
            f = meta[m]["fusion_or_driver"]
            tag = f"  {f}" if f.endswith("::FOXO1") else ""
            ylab.append((y + .5, m + tag, withdrawn,
                         "NOT in the administrator's list" in meta[m]["entity_call_basis"]))
            y += 1

    ax.set_xlim(0, ncol); ax.set_ylim(0, nrow); ax.invert_yaxis()
    ax.set_yticks([p for p, *_ in ylab])
    ax.set_yticklabels([t for _, t, *_ in ylab], fontsize=6.1)
    for lab, (_, _, wd, flag) in zip(ax.get_yticklabels(), ylab):
        if wd:   lab.set_color(NONE)
        elif flag: lab.set_color(MUTED); lab.set_style("italic")
    ax.set_xticks([i + .5 for i in range(ncol)])
    ax.set_xticklabels([c[1] for c in COLS], fontsize=6.3, rotation=55,
                       ha="left", rotation_mode="anchor")
    ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
    ax.tick_params(length=0, pad=2)
    for s in ax.spines.values():
        s.set_visible(False)

    # entity bands on the right, where single-model groups still have room to be read
    tf = mtrans.blended_transform_factory(ax.transAxes, ax.transData)
    for y0, n, ent in bands:
        ax.plot([1.012, 1.012], [y0 + .12, y0 + n - .12], color=DARK, lw=1.6,
                transform=tf, clip_on=False, solid_capstyle="butt")
        ax.text(1.026, y0 + n / 2, PRETTY.get(ent, ent).replace("\n", " "),
                va="center", ha="left", fontsize=6.6, color=DARK, weight="bold",
                transform=tf, clip_on=False)
    ax.axvline(13, color=DARK, lw=.9)

    fig.suptitle("The densest sarcoma epigenome resource that exists is complete, "
                 "uniform — and entirely behind a wall.",
                 fontsize=9.4, weight="bold", x=0.0, y=1.085, ha="left")
    universal = set.intersection(*have.values())
    n_uni = len(universal - {"WGBS", "SE*", "INPUT"})
    fig.text(0.0, 1.038,
             f"Every one of the {nrow} St Jude CSTN models carries the same "
             f"{n_uni}-mark panel plus WGBS and super-enhancer calls (CTCF only in the "
             f"2018 study). Right of the rule: what it does not contain at all.",
             fontsize=7.3, color=DARK, ha="left")

    fig.legend(handles=[
        Patch(facecolor=ACTIVE, label="active / elongation mark, or factor"),
        Patch(facecolor=REPRESS, label="repressive mark or insulator"),
        Patch(facecolor=METH, label="DNA methylation"),
        Patch(facecolor=DERIV, label="derived calls (super-enhancers)"),
        Patch(facecolor="none", edgecolor=GRID, label="absent")],
        fontsize=6.6, frameon=False, loc="upper left", bbox_to_anchor=(0.075, 0.055),
        ncol=3, handlelength=1.1, handleheight=0.8, columnspacing=1.4)

    fig.text(0.075, -0.02,
        "Access: all of it requires a St Jude / EGA data access agreement. None of it is "
        "added to any atlas total — the same material is already counted once under EGA "
        "(F25). What the\nenumeration adds is this resolution, which EGA does not publish.\n"
        "Entity calls are the St Jude data administrator's (email); the fusion partner is "
        "from the CSTN portal model table (T17). Both are recorded per track in T28.\n"
        "Red label: SJOS010930_X1, which the administrator reports as withdrawn from CSTN "
        "and no longer available — its tracks are still drawn on the browser page and are "
        "shown faded.\nGrey italic label: SJRHB010463_X16, which is on the browser page "
        "and in the portal as PAX3::FOXO1 but is absent from the administrator's list; "
        "the discrepancy is left standing.",
        fontsize=5.8, color=MUTED, ha="left", va="top", linespacing=1.6)

    os.makedirs(FIGURES, exist_ok=True)
    p = os.path.join(FIGURES, "F26_cstn_assay_matrix.pdf")
    fig.savefig(p); plt.close(fig)
    print("wrote", p)
    absent = [c[1].replace("\n", " ") for c in COLS if not any(c[0] in have[m] for m in models)]
    print(f"  {nrow} models x {ncol} assay classes")
    print(f"  present on every model: "
          f"{sorted(k for k in set.intersection(*have.values()))}")
    print(f"  absent from the whole resource: {absent}")


if __name__ == "__main__":
    main()
