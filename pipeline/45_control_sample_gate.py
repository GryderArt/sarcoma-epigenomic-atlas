#!/usr/bin/env python3
"""Stop counting ChIP input and IgG controls as if they were experiments.

Correction #24. Correction #8 excluded input controls from the regulatory counts, but
could only find them where a word boundary saw the word. Underscore is a word character,
so `_input_` -- the commonest field separator in a GEO title -- was invisible:

  GSM2967370  "160831_CHOR_input (MGH_4)"      counted as a chordoma H3K27ac profile
  GSM3100242  "GIST-T1_DMSO_input"             counted as a GIST H3K27ac profile
  GSM3326011  "TTC1240_Input_dBRD9_ChIP-Seq"   counted with target "n/a"

Each says what it is. The first two carry a mark anyway -- the depositor left the
series-wide antibody field filled in on the input, and the supplementary file is named
for the mark it controls -- which is why the title has to be read as a separator-split
string rather than as prose.

`04_classify.control_sample()` holds the rule and this stage applies it to the shipped
table; every change is logged in T38 with the word that decided it. Two words are
decisive even against a named antibody, because neither is ever a target: `input` and
`WCE`. Two are not, and are handled narrowly:

  mock  "ARID1A_MOCK_G401" is an ARID1A ChIP in mock-treated cells, not a mock IP, so
        the word is not read from a title at all.
  IgG   "Stanford_ChipSeq_IMR90_CTCF_(SC-15914)_IgG-rab" names the antibody's class on a
        real CTCF ChIP, so IgG decides only where the record names no target at all.

Run after 44, before 08. Idempotent.
"""
import csv, gzip, os, shutil, sys, collections, importlib.util
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import resolve, topen, twrite
csv.field_size_limit(10**7)


def _mod(name, fn):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(os.path.dirname(os.path.abspath(__file__)), fn))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


C4 = _mod("c04", "04_classify.py")
LOG = "T38_control_sample_corrections.tsv"
COLS = ["gsm", "gse", "decided_by", "assay_class", "target_from", "target_to",
        "antibody_raw", "disease", "entity_kind", "title", "gse_title"]
CHIPFAM = {"ChIP-seq", "CUT&RUN", "CUT&Tag", "ChIP-exo", "ChIP-chip", "HiChIP", "ChIA-PET"}
CTRL = ("input/none", "none", "None", "")


def regulatory(rows):
    """Regulatory samples per entity as every downstream stage counts them: a regulatory
    assay whose target is not a control."""
    c = collections.Counter()
    for r in rows:
        if r["entity_kind"] == "sarcoma" and r["relevance"] != "off_target" \
                and r["disease"] and r["assay_class"] in C4.REG_ASSAYS \
                and r["epi_target_norm"] not in ("input/none", "none"):
            c[r["disease"]] += 1
    return c


def main():
    path = resolve("T4_samples_atomic.tsv")
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        hdr, rows = list(rd.fieldnames), list(rd)
    before = regulatory(rows)

    log = []
    for r in rows:
        if r["assay_class"] not in CHIPFAM or r["epi_target_norm"] in CTRL:
            continue
        ab = (r.get("antibody_raw") or "").strip()
        sblob = " ;; ".join([r.get("title", ""), r.get("source_name", ""),
                             r.get("characteristics", "")])
        mk = C4.first(C4.MARK_PATS, ab + " ;; " + sblob)
        tf = C4.first(C4.TF_PATS, ab + " ;; " + sblob)
        named = bool(mk or tf or (ab and ab.lower().strip(".") not in C4._NO_ANTIBODY))
        if not C4.control_sample(r, named):
            continue
        flat = C4._SEP.sub(" ", " ;; ".join([r.get("title", ""), r.get("source_name", "")]))
        hard = C4._CONTROL_HARD.search(flat)
        a = C4.AB_RE.search(r.get("characteristics", "") or "")
        why = (hard.group(0).lower() if hard else
               "antibody: " + a.group(1).strip() if a and
               a.group(1).strip().lower().strip(".") in C4._NO_ANTIBODY else "IgG")
        frm = r["epi_target_norm"]
        r["epi_target_norm"] = "input/none"
        r["antibody_raw"] = ""
        log.append({"gsm": r["gsm"], "gse": r["gse"], "decided_by": why,
                    "assay_class": r["assay_class"], "target_from": frm,
                    "target_to": "input/none", "antibody_raw": ab,
                    "disease": r["disease"], "entity_kind": r["entity_kind"],
                    "title": r["title"][:90], "gse_title": r["gse_title"][:90]})

    done = {x["gsm"] for x in log}
    state = {r["gsm"]: r["epi_target_norm"] for r in rows}
    try:
        for x in csv.DictReader(topen(LOG), delimiter="\t"):
            if x["gsm"] not in done and state.get(x["gsm"]) == "input/none":
                log.append(x)
    except FileNotFoundError:
        pass
    log.sort(key=lambda x: (x["decided_by"], x["gse"], x["gsm"]))

    tmp = path + ".tmp"
    with op(tmp, "wt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader(); [w.writerow(r) for r in rows]
    shutil.move(tmp, path)
    with twrite(LOG, gz=False) as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader(); [w.writerow(x) for x in log]

    after = regulatory(rows)
    print(f"re-typed as controls this run: {len(done)}; {len(log)} logged -> {LOG}")
    print("  decided by:", dict(collections.Counter(
        x["decided_by"].split(":")[0] for x in log).most_common()))
    print("  had been credited to:", dict(collections.Counter(
        x["target_from"] for x in log).most_common(6)))
    print(f"\nregulatory samples: {sum(before.values()):,} -> {sum(after.values()):,}")
    for d in sorted(set(before) | set(after),
                    key=lambda d: (after[d] - before[d]) / max(before[d], 1)):
        if before[d] != after[d]:
            print(f"   {d:32s} {before[d]:5} -> {after[d]:5}   "
                  f"{100*(before[d]-after[d])/before[d]:5.1f}% control")


if __name__ == "__main__":
    main()
