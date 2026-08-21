#!/usr/bin/env python3
"""Bind the figure set into one PDF for circulation to the working group.

Kept as a stage rather than a one-off command so the circulated file cannot drift from
the figures in the repository. Pages stay vector and Arial stays embedded; pdfunite only
concatenates.
"""
import os, subprocess, glob
from _paths import FIGURES

OUT = "SASS_EPIGENETICS_FIGURES_combined.pdf"
LEGACY = "ADULT_BURDEN_FIGURES_combined.pdf"   # the name this file used to have


def main():
    pages = sorted(glob.glob(os.path.join(FIGURES, "F[0-9][0-9]_*.pdf")),
                   key=lambda p: int(os.path.basename(p)[1:3]))
    if not pages:
        raise SystemExit("no figures to combine")
    out = os.path.join(FIGURES, OUT)
    subprocess.run(["pdfunite", *pages, out], check=True)
    for p in pages:
        print("  +", os.path.basename(p))
    print(f"wrote {out}  ({os.path.getsize(out)/1e6:.2f} MB, {len(pages)} pages)")
    old = os.path.join(FIGURES, LEGACY)
    if os.path.exists(old):
        os.remove(old)
        print(f"  removed the superseded {LEGACY}")


if __name__ == "__main__":
    main()
