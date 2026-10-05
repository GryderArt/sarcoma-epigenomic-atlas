#!/usr/bin/env python3
"""WP5 -- the burden quilt: area is disease, colour is data.

Every other figure in this set asks the reader to hold two numbers in mind at once and
divide them. This one does the division geometrically. Each patch is one entity, its area
proportional to US cases per year, so the page is a picture of the disease burden; colour
then carries how much epigenomic data exists for it. A large pale patch is, by
construction, a common sarcoma nobody has profiled.

Two panels, identical tiling:

  A  colour = epigenomic samples, absolute. What a literature search would find.
  B  colour = epigenomic samples per new case per year. Coverage against need.

The pair is the point. A patch that darkens from A to B is an entity whose apparent
coverage was real; one that pales is an entity whose sample count only looked adequate
because the disease is rare. GIST is the clearest case: 360 samples reads as respectable
until it is spread over 4,670 new patients a year.

Patches are grouped into three regions by age class, sized by the burden each carries, so
the pediatric/adult asymmetry is visible before a single label is read.

Area and colour are both on log scales -- burden spans 15 to 4,670 cases and per-case
coverage spans 800-fold -- and zero is drawn in the status red rather than at the pale end
of the ramp, because "none" is a different statement from "few".
"""
import csv, os, sys, math, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "wpfig", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "28_whitepaper_figures.py"))
WP = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(WP)
plt, FIGURES, topen = WP.plt, WP.FIGURES, WP.topen
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap, to_rgb
import matplotlib.transforms as mtrans

I, F = WP.I, WP.F

# Both panels run LOW to HIGH with the loud colour at the LOW end. That inverts the usual
# sequential convention, where dark means more -- deliberately. The figure exists to show
# what is missing, so the entities with the least data have to be the ones that reach the
# eye first; a conventional ramp buries them in the palest corner.
#
# Four candidates, each validated with the data-viz skill's checker. Traffic-light
# red-amber-green is absent on merit: its poles separate by only dE 4.1 under
# deuteranopia, so the two ends of the argument would look the same to roughly one man in
# twelve. Every ramp below keeps its endpoints well clear of that.
SCHEMES = {
  "A_red_amber_pale": dict(
      label="Red → amber → pale",
      ramp=["#7f0d0d", "#b51b1b", "#d9442b", "#ef7e35", "#f8b25c", "#fcd9a0", "#f3ece0"],
      note="The alarm reading: deepest red is the least-served entity, and coverage "
           "fades to near-paper as it improves. Loudest of the four, and the one that "
           "needs no key to interpret."),
  "B_magma": dict(
      label="Magma",
      ramp=["#3b0f52", "#7d1e6d", "#b5306b", "#df5454", "#f3834a", "#fbb961", "#fcf0ba"],
      note="Perceptually uniform and colourblind-safe along its whole length, so equal "
           "steps in the data look like equal steps on the page. Reads as a heat map, "
           "which some reviewers will expect and some will find busy."),
  "C_red_grey_blue": dict(
      label="Red → grey → blue",
      ramp=["#8c1d1d", "#c0392b", "#e08b7a", "#ddd9d2", "#7fb0e4", "#2a78d6", "#0d366b"],
      note="Deficit in red, adequacy in the house blue, a grey passage between. The "
           "quietest of the four and the closest to the other four figures in the set."),
  "D_purple_gold": dict(
      label="Purple → gold",
      ramp=["#3d0f6b", "#6a2596", "#9b4a9e", "#c47191", "#dd9a6f",
            "#edc25f", "#f7e9a8"],
      note="Avoids red entirely, which some journals reserve for statistical "
           "significance. Separates well for every form of colour vision."),
}

SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
       "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("atlas_blue", SEQ)
NONE = "#d03b3b"          # status: critical. Zero is not the pale end of a ramp.
INK, MUTED, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
GAP = 0.055               # surface gap between patches, in data units

AGE = [("pediatric", "Pediatric-predominant"),
       ("both", "Both ages"),
       ("adult", "Adult-predominant")]
SHORT = {"Undifferentiated pleomorphic sarcoma": "UPS / MFH",
         "Dermatofibrosarcoma protuberans": "DFSP",
         "Desmoid / aggressive fibromatosis": "Desmoid",
         "Liposarcoma, well-differentiated": "LPS, WD",
         "Liposarcoma, dedifferentiated": "LPS, DD",
         "Liposarcoma, myxoid / round cell": "LPS, myxoid",
         "Liposarcoma, pleomorphic": "LPS, pleo",
         "Liposarcoma, NOS / mixed": "LPS, NOS",
         "GCTB / chondroblastoma": "GCTB",
         "Rhabdoid tumor / ATRT": "Rhabdoid / ATRT",
         "Inflammatory myofibroblastic tumor": "IMT",
         "Extraskeletal myxoid chondrosarcoma": "EMC",
         "Endometrial stromal sarcoma": "ESS",
         "Epithelioid hemangioendothelioma": "EHE",
         "Clear cell sarcoma of soft tissue": "CCS, soft tissue",
         "Clear cell sarcoma of the kidney": "CCSK",
         "Myoepithelial carcinoma of soft tissue": "Myoepith. carcinoma",
         "Atypical fibroxanthoma / pleomorphic dermal sarcoma": "AFX / PDS",
         "Angiomatoid fibrous histiocytoma": "AFH",
         "Ossifying fibromyxoid tumor": "OFMT",
         "FP-RMS (PAX3/7-FOXO1)": "FP-RMS",
         "FN-RMS (fusion-negative)": "FN-RMS",
         "Fibrosarcoma (adult)": "Fibrosarcoma",
         "Solitary fibrous tumor": "SFT",
         "PEComa (malignant)": "PEComa"}
short = lambda n: SHORT.get(n, n.split(" (")[0])


# ------------------------------------------------------------------ squarified treemap
def squarify(values, x, y, w, h):
    """Bruls, Huizing & van Wijk's squarified treemap.

    Returns [(x, y, w, h), ...] in the order given. Laying rows out shortest-side-first
    keeps patches near square, which is what makes their areas comparable by eye -- a
    slice-and-dice treemap of the same data produces slivers no reader can compare.
    """
    vals = list(values)
    total = sum(vals)
    if total <= 0:
        return [(x, y, 0, 0) for _ in vals]
    scale = (w * h) / total
    norm = [v * scale for v in vals]
    out, idx = [None] * len(vals), list(range(len(vals)))
    # largest first, which is what the algorithm assumes
    idx.sort(key=lambda i: -norm[i])
    rects = _squarify(([norm[i] for i in idx]), x, y, w, h)
    for i, r in zip(idx, rects):
        out[i] = r
    return out


def _worst(row, side):
    s = sum(row)
    if s <= 0 or side <= 0:
        return float("inf")
    mx, mn = max(row), min(row)
    return max((side * side * mx) / (s * s), (s * s) / (side * side * mn))


def _layout_row(row, x, y, w, h):
    """Place one row against the rectangle's shorter side and return the remainder.

    Wide rectangle -> the row becomes a column down the left edge; tall rectangle -> a
    strip across the top. Putting it against the longer side instead is what turns a
    squarified treemap back into slice-and-dice, which is how this first rendered.
    """
    s = sum(row)
    out = []
    if w >= h:                          # column down the left, width set by the row
        d = s / h if h else 0
        cy = y
        for v in row:
            ch = v / d if d else 0
            out.append((x, cy, d, ch)); cy += ch
        return out, (x + d, y, w - d, h)
    d = s / w if w else 0               # strip across the top, height set by the row
    cx = x
    for v in row:
        cw = v / d if d else 0
        out.append((cx, y, cw, d)); cx += cw
    return out, (x, y + d, w, h - d)


def _squarify(norm, x, y, w, h):
    out, row, i = [], [], 0
    while i < len(norm):
        side = min(w, h)
        if row and _worst(row, side) < _worst(row + [norm[i]], side):
            placed, (x, y, w, h) = _layout_row(row, x, y, w, h)
            out += placed; row = []
            continue
        row.append(norm[i]); i += 1
    if row:
        placed, _ = _layout_row(row, x, y, w, h)
        out += placed
    return out


# ------------------------------------------------------------------ colour
def shade(v, lo, hi, cmap=None):
    """Log-scaled position on a ramp. Zero takes the ramp's loud end and is additionally
    cross-hatched by the caller: on a scheme whose deficit end is already red, a red
    swatch for "none" says the same thing as "very little", which is the one distinction
    this figure cannot afford to lose."""
    if v <= 0:
        return (cmap if cmap is not None else CMAP)(0.0)
    c = cmap if cmap is not None else CMAP
    t = (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))
    return c(min(max(t, 0.0), 1.0))


def diverging(sch):
    """Build the two-armed map for a scheme, deficit first."""
    warm = list(sch["poles"])[::-1]          # pale -> deep, reversed to run deep -> pale
    return LinearSegmentedColormap.from_list(
        "div", list(sch["poles"]) + [sch["mid"]] + list(sch["cool"]))


def diverge_at(v, lo, hi, mid, cmap):
    """Position v on a diverging ramp whose neutral sits at `mid`, each arm scaled to its
    own span so the midpoint lands at 0.5 even when the arms differ in length."""
    if v <= 0:
        return NONE
    lv, lm = math.log10(v), math.log10(mid)
    if lv <= lm:
        span = lm - math.log10(lo)
        t = 0.5 * (lv - math.log10(lo)) / span if span else 0.5
    else:
        span = math.log10(hi) - lm
        t = 0.5 + 0.5 * (lv - lm) / span if span else 0.5
    return cmap(min(max(t, 0.0), 1.0))


def ink_on(rgb):
    """Label colour by the patch's luminance, so text never sits at 2:1 on its own fill."""
    if isinstance(rgb, str):
        rgb = to_rgb(rgb)
    lum = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return "#ffffff" if lum < 0.55 else INK


# ------------------------------------------------------------------ one panel
def panel(ax, ents, norate, key, lo, hi, title, sub, fmt, col_of=None):
    """Draw the quilt once. `key` picks the value that colours each patch."""
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    pending = []
    STRIP = 7.2                       # footer for entities with no published rate
    body_h = 100 - STRIP - 1.6

    # Outer tiling: three age-class regions, sized by the burden each carries.
    groups = [(lab, [e for e in ents if e["age"] == a]) for a, lab in AGE]
    gcases = [sum(e["cases"] for e in g) for _, g in groups]   # warped: area weight
    outer = squarify(gcases, 0, 0, 100, body_h)

    for (lab, g), (gx, gy, gw, gh) in zip(groups, outer):
        pad = 0.55
        hx, hy, hw, hh = gx + pad, gy + pad + 3.0, gw - 2 * pad, gh - 2 * pad - 3.0
        ax.add_patch(Rectangle((gx, gy), gw, gh, facecolor="none",
                               edgecolor="#c3c2b7", linewidth=0.8, zorder=4))
        ax.text(gx + 1.0, gy + 2.2, " ".join(lab.upper()), fontsize=6.0, weight="bold",
                color=MUTED, va="center", ha="left", zorder=5)
        tot = sum(e.get("cases_true", e["cases"]) for e in g)
        ax.text(gx + gw - 1.0, gy + 2.2, f"{tot:,.0f} cases/yr", fontsize=5.8,
                color=MUTED, va="center", ha="right", zorder=5)

        g = sorted(g, key=lambda e: -e["cases"])
        for e, (px, py, pw, ph) in zip(g, squarify([e["cases"] for e in g],
                                                   hx, hy, hw, hh)):
            v = e[key]
            col = col_of(v) if col_of else shade(v, lo, hi)
            ax.add_patch(Rectangle((px + GAP, py + GAP),
                                   max(pw - 2 * GAP, 0), max(ph - 2 * GAP, 0),
                                   facecolor=col, edgecolor="none", zorder=3,
                                   hatch="xxx" if v <= 0 else None,
                                   linewidth=0))
            pending.append((e, px, py, pw, ph, col, fmt(v)))

    # Entities with no population rate have no area to occupy; they get a strip rather
    # than a fabricated denominator.
    sy = body_h + 1.6
    ax.add_patch(Rectangle((0, sy), 100, STRIP, facecolor="none",
                           edgecolor="#c3c2b7", linewidth=0.8, linestyle=(0, (2, 2)),
                           zorder=4))
    ax.text(0.9, sy + 1.9, "NO PUBLISHED INCIDENCE RATE — SHOWN AT EQUAL WIDTH",
            fontsize=5.4, color=MUTED, va="center", ha="left", zorder=5)
    cw = 100 / max(len(norate), 1)
    for i, e in enumerate(norate):
        # Panel B is a ratio, and these entities have no denominator. Colouring them red
        # there would say "no data" when what is missing is the incidence rate: three of
        # the four do have epigenomic samples. They are drawn unfilled instead.
        if key == "epi":
            col = col_of(e["epi"]) if col_of else shade(e["epi"], lo, hi)
            txt = fmt(e["epi"])
        else:
            col, txt = "#f0efec", "no incidence rate"
        px = i * cw
        ax.add_patch(Rectangle((px + GAP, sy + 3.1), cw - 2 * GAP, STRIP - 3.1 - GAP,
                               facecolor=col, edgecolor="none", zorder=3))
        ic = ink_on(col)
        ax.text(px + cw / 2, sy + 4.6, short(e["name"]),
                fontsize=5.0, color=ic, ha="center", va="center", zorder=5)
        ax.text(px + cw / 2, sy + 6.0, txt, fontsize=5.0, weight="bold",
                color=ic, ha="center", va="center", zorder=5)

    place_labels(ax, pending)
    ax.set_title(title, fontsize=8.8, weight="bold", loc="left", pad=11, color=INK)
    ax.text(0, -0.9, sub, transform=mtrans.blended_transform_factory(
        ax.transAxes, ax.transData), fontsize=6.6, color=MUTED, va="bottom", ha="left")


def place_labels(ax, pending):
    """Name and value inside each patch, measured rather than estimated.

    A character-count estimate put "Leiomyosarcoma" on a patch that could hold it and
    reduced it to "L" -- and to a second "L" elsewhere on the page. So each candidate is
    drawn and measured with the real renderer, exactly as the scatter labeller does, and
    the first that fits inside the patch wins. A patch that fits nothing stays blank; its
    entity is still recoverable from the table and from the interactive map.
    """
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    probe = ax.text(0, 0, "", ha="center", va="center", zorder=7)

    for e, px, py, pw, ph, col, val in pending:
        # patch size in display pixels
        (x0, y0), (x1, y1) = ax.transData.transform([(px, py), (px + pw, py + ph)])
        bw, bh = abs(x1 - x0) - 4, abs(y1 - y0) - 4
        if bw < 14 or bh < 9:
            continue
        full = short(e["name"])
        initials = "".join(w[0] for w in full.replace("/", " ").replace("-", " ").split()
                           if w[:1].isalnum()).upper()
        ic = ink_on(col)
        best = None
        for text in ([full] if len(initials) < 2 else [full, initials]):
            for fs in (9.5, 8.0, 7.0, 6.2, 5.6, 5.0, 4.4, 4.0):
                probe.set_text(text); probe.set_fontsize(fs)
                bb = probe.get_window_extent(rend)
                if bb.width <= bw and bb.height * 2.25 <= bh:
                    best = (text, fs, True); break        # room for a value line too
                if bb.width <= bw and bb.height <= bh:
                    best = best or (text, fs, False)
            if best and best[2]:
                break
        if not best:
            continue
        text, fs, two = best
        cx, cy = px + pw / 2, py + ph / 2
        dy = ph * 0.145 if two else 0
        ax.text(cx, cy - dy, text, fontsize=fs, color=ic, ha="center", va="center",
                zorder=6)
        if two:
            ax.text(cx, cy + dy * 1.45, val, fontsize=fs * 0.82, weight="bold",
                    color=ic, ha="center", va="center", zorder=6, alpha=0.93)
    probe.remove()


# ------------------------------------------------------------------ build
def load():
    T13 = list(csv.DictReader(topen("T13_incidence_vs_data.tsv"), delimiter="\t"))
    ents, norate = [], []
    for r in T13:
        c = F(r, "US_cases_per_year_all_ages")
        e = {"name": r["display_name"], "age": r["age_class"], "cases": c,
             "epi": I(r, "epigenomic_samples"),
             "reg": I(r, "regulatory_epigenomic_samples"),
             "meth": I(r, "DNA_methylation_samples")}
        e["per"] = (e["epi"] / c) if c else 0.0
        (ents if c else norate).append(e)
    return ents, norate


def build(power, stem, note, scheme=None):
    """One figure. `power` warps patch area: 1.0 is true area, below 1 compresses the
    range so rare entities stay legible. The exponent is stated on the figure, because an
    area encoding that is not proportional to its quantity and does not say so is a lie.
    `scheme` selects a diverging palette for panel B; None keeps the one-hue ramp."""
    ents, norate = load()
    for e in ents:
        e["cases"] = e["cases_true"] = e.get("cases_true", e["cases"])
    sized = [dict(e, cases=e["cases_true"] ** power) for e in ents]

    epis = [e["epi"] for e in ents + norate if e["epi"] > 0]
    pers = [e["per"] for e in ents if e["per"] > 0]
    elo, ehi = min(epis), max(epis)
    plo, phi = min(pers), max(pers)
    # The reference: samples per case if the whole atlas were spread evenly over the
    # whole burden. Computed, not chosen.
    mid = sum(e["epi"] for e in ents) / sum(e["cases_true"] for e in ents)

    sch = SCHEMES[scheme] if scheme else None
    # One ramp, both panels, loud end at the low end: least data reaches the eye first.
    cmap = LinearSegmentedColormap.from_list("sch", sch["ramp"]) if sch else CMAP
    colA = lambda v: shade(v, elo, ehi, cmap)
    colB = lambda v: shade(v, plo, phi, cmap)

    fig, axes = plt.subplots(1, 2, figsize=(14.2, 7.4))
    fig.subplots_adjust(left=0.012, right=0.988, top=0.800, bottom=0.120, wspace=0.055)

    panel(axes[0], sized, norate, "epi", elo, ehi,
          "A   Coloured by epigenomic samples",
          f"absolute count, log scale \u00b7 {elo:,} to {ehi:,} \u00b7 loud = least data",
          lambda v: f"{v:,.0f}" if v else "0", colA)
    panel(axes[1], sized, norate, "per", plo, phi,
          "B   Coloured by samples per new case per year",
          f"coverage against need, log scale \u00b7 {plo:.2f} to {phi:.1f} \u00b7 "
          f"loud = least data",
          lambda v: f"{v:.2f}" if v and v < 10 else (f"{v:.0f}" if v else "0"), colB)

    ratio = (max(e["cases_true"] for e in ents) / min(e["cases_true"] for e in ents))
    shown = ratio ** power
    fig.suptitle("The quilt of unmet need: patch area is how common the disease is, "
                 "colour is how much of its chromatin anyone has read",
                 fontsize=11.2, weight="bold", x=0.012, y=0.980, ha="left", color=INK)
    sub = ("Each patch is one of the 45 named entities, grouped into three regions by "
           "age class and sized within them by US incidence. ")
    sub += ("The same tiling is coloured two ways. Colour runs loud-to-quiet as coverage "
            "improves, so the entities with\nthe least data are the ones that reach the "
            "eye first; a patch that quietens from A to B had only looked well-covered "
            "because the disease is rare.")
    fig.text(0.012, 0.952, sub, fontsize=7.6, color=MUTED, ha="left", va="top",
             linespacing=1.55)

    legend(fig, elo, ehi, plo, phi, sch, cmap)
    WP.foot(fig, 0.012, 0.080,
            note + " Colour is log-scaled in both panels and runs loud to quiet as "
            "coverage improves; an entity with no epigenomic data of any kind takes the "
            "loud end and is cross-hatched, "
            "because none is a different statement from few. Epigenomic means every "
            "assay in the atlas, DNA methylation included. Incidence is US cases per "
            f"year over the 45 named entities; {len(norate)} "
            f"{'entity has' if len(norate) == 1 else 'entities have'} no population rate "
            "published anywhere and "
            f"{'is' if len(norate) == 1 else 'are'} shown at equal width in the strip "
            "beneath, excluded from every area and from panel B.", width=168, size=6.2)
    p = os.path.join(FIGURES, stem + ".pdf")
    fig.savefig(p); plt.close(fig)
    print(f"  wrote {os.path.basename(p)}")
    return shown, ratio


def legend(fig, elo, ehi, plo, phi, sch=None, cmap=None):
    """Two ramps and the zero swatch, drawn as discrete steps rather than a colorbar so
    the log stepping is legible, and so the diverging midpoint is a mark a reader can
    point at rather than a position they have to infer."""
    n = 9
    # --- panel A, sequential
    x0, y0 = 0.012, 0.878
    fig.text(x0, y0 + 0.030, "samples (panel A)", fontsize=6.2, color=MUTED,
             ha="left", va="bottom")
    for k in range(n):
        fig.patches.append(plt.Rectangle(
            (x0 + k * 0.0155, y0), 0.0145, 0.016, transform=fig.transFigure,
            facecolor=(cmap or CMAP)(k / (n - 1)), edgecolor="none", zorder=5))
    fig.text(x0, y0 - 0.004, f"{elo:,.0f}", fontsize=5.6, color=MUTED, ha="left", va="top")
    fig.text(x0 + n * 0.0155 - 0.002, y0 - 0.004, f"{ehi:,.0f}", fontsize=5.6,
             color=MUTED, ha="right", va="top")

    # --- panel B
    x1 = 0.012 + 0.215
    lab = "per case (panel B)  \u2014 least data on the left"
    fig.text(x1, y0 + 0.030, lab, fontsize=6.2, color=MUTED, ha="left", va="bottom")
    m = 11
    for k in range(m):
        c = (cmap or CMAP)(k / (m - 1))
        fig.patches.append(plt.Rectangle(
            (x1 + k * 0.0155, y0), 0.0145, 0.016, transform=fig.transFigure,
            facecolor=c, edgecolor="none", zorder=5))
    fig.text(x1, y0 - 0.004, f"{plo:.2f}", fontsize=5.6, color=MUTED, ha="left", va="top")
    fig.text(x1 + m * 0.0155 - 0.002, y0 - 0.004, f"{phi:.1f}", fontsize=5.6,
             color=MUTED, ha="right", va="top")

    x2 = 0.012 + 2 * 0.215
    fig.patches.append(plt.Rectangle((x2, y0), 0.0145, 0.016, transform=fig.transFigure,
                                     facecolor=(cmap or CMAP)(0.0), edgecolor="none",
                                     hatch="xxx", linewidth=0, zorder=5))
    fig.text(x2 + 0.020, y0 + 0.008, "cross-hatched: no epigenomic data of any kind",
             fontsize=6.2, color=MUTED, ha="left", va="center")


BUF = ("Patch area is proportional to (US cases per year) raised to the power 0.55, "
       "which compresses a 311-fold spread in incidence to roughly 24-fold so that "
       "the rarest entities remain legible. Area still ranks burden correctly, but "
       "differences in area understate differences in incidence; the case count is "
       "printed on every patch that has room for it.")

# Chosen for the paper: red for the deficit, the house blue for adequacy, a grey passage
# between. The quietest of the four candidates and the closest to WP1-WP4, so the figure
# set still reads as one system.
CHOSEN = "C_red_grey_blue"

if __name__ == "__main__":
    import sys
    if "--schemes" in sys.argv:
        print("WP5 colour-scheme candidates:")
        for k, v in SCHEMES.items():
            build(0.55, "WP5_scheme_" + k, BUF, scheme=k)
            print(f"      {v['label']}: {v['note']}")
        raise SystemExit
    print("WP5 burden quilt:")
    build(1.00, "WP5_burden_quilt",
          "Patch area is directly proportional to US cases per year.", scheme=CHOSEN)
    build(0.55, "WP5b_burden_quilt_buffered", scheme=CHOSEN, note=
          "Patch area is proportional to (US cases per year) raised to the power 0.55, "
          "which compresses a 311-fold spread in incidence to roughly 24-fold so that "
          "the rarest entities remain legible. Area still ranks burden correctly, but "
          "differences in area understate differences in incidence; the case count is "
          "printed on every patch that has room for it.")
