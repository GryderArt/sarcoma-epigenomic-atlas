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

# Sequential blue, steps 100-700 of the validated ramp (docs/methods.md; the same scale
# the gap map's coverage matrix uses, so the two read as one system).
RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
        "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("atlas_blue", RAMP)
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
def shade(v, lo, hi):
    """Log-scaled position on the ramp; zero is handled by the caller."""
    if v <= 0:
        return NONE
    t = (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))
    return CMAP(min(max(t, 0.0), 1.0))


def ink_on(rgb):
    """Label colour by the patch's luminance, so text never sits at 2:1 on its own fill."""
    if isinstance(rgb, str):
        rgb = to_rgb(rgb)
    lum = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return "#ffffff" if lum < 0.55 else INK


# ------------------------------------------------------------------ one panel
def panel(ax, ents, norate, key, lo, hi, title, sub, fmt):
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
            col = shade(v, lo, hi)
            ax.add_patch(Rectangle((px + GAP, py + GAP),
                                   max(pw - 2 * GAP, 0), max(ph - 2 * GAP, 0),
                                   facecolor=col, edgecolor="none", zorder=3))
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
            col, txt = shade(e["epi"], lo, hi), fmt(e["epi"])
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


def build(power, stem, note):
    """One figure. `power` warps patch area: 1.0 is true area, below 1 compresses the
    range so rare entities stay legible. The exponent is stated on the figure, because an
    area encoding that is not proportional to its quantity and does not say so is a lie."""
    ents, norate = load()
    for e in ents:
        e["cases"] = e["cases_true"] = e.get("cases_true", e["cases"])
    sized = [dict(e, cases=e["cases_true"] ** power) for e in ents]

    epis = [e["epi"] for e in ents + norate if e["epi"] > 0]
    pers = [e["per"] for e in ents if e["per"] > 0]
    elo, ehi = min(epis), max(epis)
    plo, phi = min(pers), max(pers)

    fig, axes = plt.subplots(1, 2, figsize=(14.2, 7.4))
    fig.subplots_adjust(left=0.012, right=0.988, top=0.800, bottom=0.120, wspace=0.055)

    panel(axes[0], sized, norate, "epi", elo, ehi,
          "A   Coloured by epigenomic samples",
          f"absolute count, log scale · {elo:,} to {ehi:,} · red = none",
          lambda v: f"{v:,.0f}" if v else "0")
    panel(axes[1], sized, norate, "per", plo, phi,
          "B   Coloured by samples per new case per year",
          f"coverage against need, log scale · {plo:.2f} to {phi:.1f} · red = none",
          lambda v: f"{v:.2f}" if v and v < 10 else (f"{v:.0f}" if v else "0"))

    ratio = (max(e["cases_true"] for e in ents) / min(e["cases_true"] for e in ents))
    shown = ratio ** power
    fig.suptitle("The quilt of unmet need: patch area is how common the disease is, "
                 "colour is how much of its chromatin anyone has read",
                 fontsize=11.2, weight="bold", x=0.012, y=0.980, ha="left", color=INK)
    fig.text(0.012, 0.952,
             "Each patch is one of the 45 named entities, grouped into three regions by "
             "age class and sized within them by US incidence. The same tiling is "
             "coloured two ways: a patch that pales from A to B\nhad only looked "
             "well-covered because the disease is rare.",
             fontsize=7.6, color=MUTED, ha="left", va="top", linespacing=1.55,
             wrap=False)

    legend(fig, elo, ehi, plo, phi)
    WP.foot(fig, 0.012, 0.080,
            note + " Colour is log-scaled in both panels; an entity with no epigenomic "
            "data of any kind is drawn in red rather than at the pale end of the ramp, "
            "because none is a different statement from few. Epigenomic means every "
            "assay in the atlas, DNA methylation included. Incidence is US cases per "
            "year over the 45 named entities; four entities have no population rate "
            "published anywhere and are shown at equal width in the strip beneath, "
            "excluded from every area and from panel B.", width=168, size=6.2)
    p = os.path.join(FIGURES, stem + ".pdf")
    fig.savefig(p); plt.close(fig)
    print(f"  wrote {os.path.basename(p)}   area exponent {power:g} · "
          f"largest:smallest patch {shown:,.0f}:1 (true burden ratio {ratio:,.0f}:1)")
    return shown, ratio


def legend(fig, elo, ehi, plo, phi):
    """Two ramps and the zero swatch, drawn as real marks rather than a colorbar so the
    log steps are legible as steps."""
    for i, (lab, lo, hi, fmt) in enumerate(
            [("samples (panel A)", elo, ehi, lambda v: f"{v:,.0f}"),
             ("per case (panel B)", plo, phi,
              lambda v: f"{v:.2f}".rstrip("0").rstrip("."))]):
        x0 = 0.012 + i * 0.215
        y0 = 0.878
        fig.text(x0, y0 + 0.030, lab, fontsize=6.2, color=MUTED, ha="left", va="bottom")
        n = 9
        for k in range(n):
            t = k / (n - 1)
            v = 10 ** (math.log10(lo) + t * (math.log10(hi) - math.log10(lo)))
            fig.patches.append(plt.Rectangle(
                (x0 + k * 0.0155, y0), 0.0145, 0.016, transform=fig.transFigure,
                facecolor=CMAP(t), edgecolor="none", zorder=5))
        fig.text(x0, y0 - 0.004, fmt(lo), fontsize=5.6, color=MUTED, ha="left", va="top")
        fig.text(x0 + n * 0.0155 - 0.002, y0 - 0.004, fmt(hi), fontsize=5.6,
                 color=MUTED, ha="right", va="top")
    x0 = 0.012 + 2 * 0.215
    fig.patches.append(plt.Rectangle((x0, 0.878), 0.0145, 0.016,
                                     transform=fig.transFigure, facecolor=NONE,
                                     edgecolor="none", zorder=5))
    fig.text(x0 + 0.020, 0.886, "no epigenomic data of any kind",
             fontsize=6.2, color=MUTED, ha="left", va="center")


if __name__ == "__main__":
    print("WP5 burden quilt:")
    build(1.00, "WP5_burden_quilt",
          "Patch area is directly proportional to US cases per year.")
    build(0.55, "WP5b_burden_quilt_buffered",
          "Patch area is proportional to (US cases per year) raised to the power 0.55, "
          "which compresses a 311-fold spread in incidence to roughly 24-fold so that "
          "the rarest entities remain legible. Area still ranks burden correctly, but "
          "differences in area understate differences in incidence; the case count is "
          "printed on every patch that has room for it.")
