"""Collision-aware point labeling for scatter panels, with leader lines.

Why this exists. Panel D of WP1 originally placed labels by testing a guessed 52x11 pixel
box against previously placed labels. "Rhabdoid tumor / ATRT" (95 US cases per year,
1,052 regulatory samples) and "FP-RMS" (110, 640) sit about 50 display pixels apart; both
took the slot directly above themselves, their guessed boxes cleared each other by a
pixel, and the result read as a single two-line label hovering over the upper point. The
two entities were never pooled -- they are separate rows in every atlas table -- but the
figure could not be used to tell.

Four rules fix that, and all four matter:

  1. Measure, never guess. Each candidate is drawn and measured with the real renderer.
     Note that Annotation.get_window_extent includes the leader line, which would make
     every box tall enough to collide with everything -- so the search uses a bare label
     and the leader is attached only to the winner.
  2. Ownership. A candidate is rejected if the label box sits nearer to somebody else's
     marker than to its own. This is the rule that actually separates a stacked pair.
  3. Direction. Candidates are tried in order of how well they point away from the
     point's nearest neighbour, so crowded pairs open outwards instead of both going up.
  4. A leader line on every label, so the reading never depends on proximity alone.

Placement that cannot be made clean is reported to stdout rather than silently drawn.
"""
CAND = [(0, 7, "center", "bottom"), (0, -7, "center", "top"),
        (7, 0, "left", "center"), (-7, 0, "right", "center"),
        (6, 6, "left", "bottom"), (-6, 6, "right", "bottom"),
        (6, -6, "left", "top"), (-6, -6, "right", "top"),
        (0, 15, "center", "bottom"), (0, -15, "center", "top"),
        (15, 7, "left", "bottom"), (-15, 7, "right", "bottom"),
        (15, -7, "left", "top"), (-15, -7, "right", "top"),
        (0, 24, "center", "bottom"), (0, -24, "center", "top"),
        (26, 14, "left", "bottom"), (-26, 14, "right", "bottom"),
        (26, -14, "left", "top"), (-26, -14, "right", "top"),
        (30, 0, "left", "center"), (-30, 0, "right", "center"),
        (0, 34, "center", "bottom"), (0, -34, "center", "top")]


def _near(bb, qx, qy):
    """Distance in pixels from a point to the nearest edge of a label box."""
    dx = max(bb.x0 - qx, 0, qx - bb.x1)
    dy = max(bb.y0 - qy, 0, qy - bb.y1)
    return (dx * dx + dy * dy) ** 0.5


def place_labels(fig, ax, items, points, avoid=(), fontsize=6, color="#2b3440",
                 leader="#7c8894", pad=10.0, own_margin=1.35, where=""):
    """Label `items` = [(x, y, text), ...] on `ax`, avoiding every marker in
    `points` = [(x, y), ...] and every artist bbox in `avoid`.

    Returns the list of label bounding boxes actually used.
    """
    arrow = dict(arrowstyle="-", color=leader, linewidth=0.5, shrinkA=0.5, shrinkB=2.5)
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    marks = [ax.transData.transform(p) for p in points]
    placed = [a.get_window_extent(rend) for a in avoid]
    frame = ax.get_window_extent(rend)

    def clashes(bb, ox, oy):
        if (bb.x0 < frame.x0 or bb.x1 > frame.x1
                or bb.y0 < frame.y0 or bb.y1 > frame.y1):
            return True                               # a label that runs off the panel
        for b in placed:
            if (bb.x0 < b.x1 + pad and bb.x1 > b.x0 - pad
                    and bb.y0 < b.y1 + pad and bb.y1 > b.y0 - pad):
                return True
        mine = _near(bb, ox, oy)
        for mx, my in marks:
            if abs(mx - ox) < 0.5 and abs(my - oy) < 0.5:
                continue                              # our own marker
            if bb.x0 - 3 < mx < bb.x1 + 3 and bb.y0 - 3 < my < bb.y1 + 3:
                return True                           # never cover a foreign point
            if _near(bb, mx, my) < mine * own_margin + 2:
                return True                           # nearer a stranger than to us
        return False

    def draw(n, x, yv, dxo, dyo, ha, va, arrowed=False):
        return ax.annotate(n, (x, yv), fontsize=fontsize, color=color,
                           textcoords="offset points", xytext=(dxo, dyo), ha=ha, va=va,
                           zorder=5, **({"arrowprops": arrow} if arrowed else {}))

    for x, yv, n in sorted(items, key=lambda t: -t[1]):
        px, py = ax.transData.transform((x, yv))
        others = [m for m in marks if abs(m[0] - px) > 0.5 or abs(m[1] - py) > 0.5]
        if others:
            nx_, ny_ = min(others, key=lambda m: (m[0] - px) ** 2 + (m[1] - py) ** 2)
            ax_, ay_ = px - nx_, py - ny_
            mag = (ax_ * ax_ + ay_ * ay_) ** 0.5 or 1.0
            ax_, ay_ = ax_ / mag, ay_ / mag
            cands = sorted(CAND, key=lambda c: -((c[0] * ax_ + c[1] * ay_)
                                                 / ((c[0] ** 2 + c[1] ** 2) ** 0.5 or 1.0)))
        else:
            cands = list(CAND)

        chosen = None
        for dxo, dyo, ha, va in cands:
            t = draw(n, x, yv, dxo, dyo, ha, va)
            bb = t.get_window_extent(rend)
            t.remove()
            if not clashes(bb, px, py):
                chosen = (bb, dxo, dyo, ha, va)
                break
        if chosen is None:
            # Take the position that overlaps least, rather than a fixed fallback.
            best, bestpen = None, None
            for dxo, dyo, ha, va in cands:
                t = draw(n, x, yv, dxo, dyo, ha, va)
                bb = t.get_window_extent(rend)
                t.remove()
                pen = sum(max(0, min(bb.x1, b.x1) - max(bb.x0, b.x0))
                          * max(0, min(bb.y1, b.y1) - max(bb.y0, b.y0)) for b in placed)
                pen += 600 * (bb.x0 < frame.x0 or bb.x1 > frame.x1
                              or bb.y0 < frame.y0 or bb.y1 > frame.y1)
                pen += sum(400 for mx, my in marks
                           if (abs(mx - px) > 0.5 or abs(my - py) > 0.5)
                           and bb.x0 - 3 < mx < bb.x1 + 3 and bb.y0 - 3 < my < bb.y1 + 3)
                if bestpen is None or pen < bestpen:
                    best, bestpen = (bb, dxo, dyo, ha, va), pen
            chosen = best
            print(f"  !! {where or 'scatter'}: no clean position for {n!r} "
                  f"-- labels may overlap")
        bb, dxo, dyo, ha, va = chosen
        draw(n, x, yv, dxo, dyo, ha, va, arrowed=True)
        placed.append(bb)
    return placed
