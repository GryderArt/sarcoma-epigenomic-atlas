"""Shared paths and table IO for the sarcoma epigenomic atlas pipeline.

Every stage reads and writes through here, so the pipeline runs from a checkout without
editing paths. Override the root with SASS_ROOT if you keep data elsewhere.

Large tables are stored gzipped in the repository; `topen` and `twrite` handle that
transparently, so the stage scripts never need to know which form is on disk.
"""
import os, gzip

ROOT = os.environ.get(
    "SASS_ROOT",
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA      = os.environ.get("SASS_DATA",      os.path.join(ROOT, "data"))
SAMPLES   = os.path.join(DATA, "samples")
INCIDENCE = os.path.join(DATA, "incidence")
FIGURES   = os.environ.get("SASS_FIGURES",   os.path.join(ROOT, "figures"))
WORKBOOKS = os.environ.get("SASS_WORKBOOKS", os.path.join(ROOT, "workbooks"))
DOCS      = os.environ.get("SASS_DOCS",      os.path.join(ROOT, "docs"))
# scratch for multi-gigabyte harvest intermediates; not tracked in git
WORK      = os.environ.get("SASS_WORK",      os.path.join(ROOT, ".work"))

for _d in (DATA, SAMPLES, FIGURES, WORKBOOKS, DOCS, WORK):
    os.makedirs(_d, exist_ok=True)


def _candidates(name):
    """Where a table might live, in preference order."""
    if os.path.isabs(name):
        yield name; yield name + ".gz"; return
    base = os.path.basename(name)
    for d in (DATA, SAMPLES, WORK):
        yield os.path.join(d, base)
        yield os.path.join(d, base + ".gz")


def resolve(name):
    for p in _candidates(name):
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        f"{name} not found under {DATA} or {WORK}. "
        f"Large tables are gzipped in the repository; run the earlier pipeline stage "
        f"if this is an intermediate.")


def topen(name):
    """Open a data table for reading, transparently handling gzip."""
    p = resolve(name)
    return gzip.open(p, "rt", newline="") if p.endswith(".gz") \
        else open(p, "r", newline="")


def twrite(name, gz=None, work=False):
    """Open a data table for writing. Tables over a few MB are gzipped by default so the
    repository stays reviewable; pass gz=False to force plain text."""
    d = WORK if work else DATA
    p = name if os.path.isabs(name) else os.path.join(d, os.path.basename(name))
    if gz is None:
        gz = p.endswith(".gz")
    if gz and not p.endswith(".gz"):
        p += ".gz"
    os.makedirs(os.path.dirname(p), exist_ok=True)
    return gzip.open(p, "wt", newline="") if gz else open(p, "w", newline="")


def dpath(name, work=False):
    """Path a table should be written to (for callers that need the name, not a handle)."""
    return os.path.join(WORK if work else DATA, os.path.basename(name))
