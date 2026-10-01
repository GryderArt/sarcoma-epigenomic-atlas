#!/usr/bin/env python3
"""One-shot pass converting the atlas's own English to US spelling.

Rewrites entity labels and our own prose in the data tables, the pipeline source and
the documentation. Leaves untouched: the verbatim contents of third-party records, and
every regular expression that matches them -- see _enus.py for why both matter.

Idempotent: running it twice changes nothing the second time.
"""
import csv, gzip, os, sys, glob, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _enus import convert, convert_python, VERBATIM_COLS, VERBATIM_FILES
from _paths import ROOT, DATA
csv.field_size_limit(10**7)


def do_tsv(path):
    base = os.path.basename(path).replace(".gz", "")
    if base in VERBATIM_FILES:
        return 0, "skipped (third-party inventory)"
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as f:
        rd = csv.reader(f, delimiter="\t")
        rows = list(rd)
    if not rows:
        return 0, "empty"
    head = rows[0]
    safe = [i for i, c in enumerate(head) if c not in VERBATIM_COLS]
    n = 0
    for r in rows[1:]:
        for i in safe:
            if i < len(r) and r[i]:
                v = convert(r[i])
                if v != r[i]:
                    r[i] = v
                    n += 1
    if not n:
        return 0, "unchanged"
    tmp = path + ".tmp"
    with (gzip.open if path.endswith(".gz") else open)(tmp, "wt", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerows(rows)
    shutil.move(tmp, path)
    held = [c for c in head if c in VERBATIM_COLS]
    return n, f"{n} cells" + (f"; held {', '.join(held)}" if held else "")


def do_text(path, python=False):
    src = open(path, encoding="utf-8").read()
    out = convert_python(src) if python else convert(src)
    if out == src:
        return 0
    open(path, "w", encoding="utf-8").write(out)
    return sum(1 for a, b in zip(src.split("\n"), out.split("\n")) if a != b)


def main():
    print("data tables")
    tot = 0
    for p in sorted(glob.glob(os.path.join(DATA, "**", "*.tsv*"), recursive=True)):
        n, note = do_tsv(p)
        tot += n
        if n or "held" in note or "skipped" in note:
            print(f"  {os.path.basename(p):44s} {note}")
    print(f"  -> {tot:,} cells rewritten\n")

    print("pipeline source (raw-string literals masked)")
    for p in sorted(glob.glob(os.path.join(ROOT, "pipeline", "*.py"))):
        if os.path.basename(p) in ("_enus.py", "29_normalize_spelling.py"):
            continue
        n = do_text(p, python=True)
        if n:
            print(f"  {os.path.basename(p):44s} {n} lines")

    print("\ndocumentation")
    for p in (sorted(glob.glob(os.path.join(ROOT, "docs", "*.md")))
              + [os.path.join(ROOT, "README.md"),
                 os.path.join(ROOT, "data", "README.md"),
                 os.path.join(ROOT, "pipeline", "gapmap_template.html")]):
        if os.path.exists(p):
            n = do_text(p)
            if n:
                print(f"  {os.path.relpath(p, ROOT):44s} {n} lines")


if __name__ == "__main__":
    main()
