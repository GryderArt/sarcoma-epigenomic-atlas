#!/usr/bin/env python3
"""Render README.md from the atlas tables.

The README is the first thing anyone sees, and it was the last thing anyone re-ran. By
1 September 2026 it carried a sample total, a study count, a per-case ratio and a
GSE140686 entity count from three different runs -- every one of them wrong, and none of
them wrong in a way a reader could detect. A document that states numbers it cannot
re-derive is a liability in a repository whose entire claim is reproducibility.

So the prose lives in README.tmpl.md with {placeholders}, and this stage fills them from
the tables. Any figure the tables stop supporting fails here, loudly, instead of sitting
in the README for three runs.

Run after 33.
"""
import csv, json, os, re, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import ROOT, DOCS, topen
csv.field_size_limit(10**7)

REG = {"ChIP-seq", "CUT&RUN", "CUT&Tag", "ChIP-exo", "ChIP-chip", "ATAC-seq", "scATAC-seq",
       "DNase-seq", "FAIRE-seq", "MNase-seq", "Hi-C", "HiChIP", "Micro-C", "Capture-HiC",
       "ChIA-PET", "4C-seq"}
PDERIVED = {"primary_tumor", "metastasis", "recurrence", "PDX", "organoid"}
I = lambda r, k: int(float(r.get(k) or 0))
F = lambda r, k: float(r.get(k) or 0)
N = lambda v: f"{int(round(v)):,}"


def rewrap(md, width=92):
    """Re-flow prose paragraphs after substitution.

    The template is hand-wrapped, and a substituted value is almost never the width of
    its placeholder, so without this the rendered file has 40-character lines next to
    300-character ones. Tables, fenced code and link-only lines are left alone.
    """
    import textwrap
    # Markdown links and file paths are single tokens; breaking inside one silently
    # corrupts the link, which is exactly what the first run of this did.
    OPT = dict(break_long_words=False, break_on_hyphens=False)
    out, buf, fence = [], [], False

    def flush():
        if not buf:
            return
        para = " ".join(" ".join(buf).split())
        first = buf[0]
        if first.lstrip().startswith(("- ", "* ")) or first.lstrip()[:2].rstrip(".").isdigit():
            ind = " " * (len(first) - len(first.lstrip()) + 2)
            out.extend(textwrap.wrap(para, width, subsequent_indent=ind, **OPT))
        else:
            out.extend(textwrap.wrap(para, width, **OPT))
        buf.clear()

    for line in md.split("\n"):
        if line.startswith("```"):
            flush(); fence = not fence; out.append(line); continue
        if fence or line.startswith(("|", "#", "---")) or not line.strip():
            flush(); out.append(line); continue
        if line.lstrip().startswith(("- ", "* ")) or (buf and not line.startswith(" ") and
                                                      buf[0].lstrip().startswith(("- ", "* "))):
            if line.lstrip().startswith(("- ", "* ")):
                flush()
            buf.append(line); continue
        buf.append(line)
    flush()
    return "\n".join(out)


def refresh_data_readme():
    """Re-measure the Size and Rows columns of data/README.md against the files on disk.

    The prose in that file is hand-written and stays hand-written; only the measurements
    are regenerated, plus any row count the prose repeats from its own row. After the
    delta harvest T4 was described as 85,698 rows in a file sitting beside a 92,138-row
    table.
    """
    import gzip, re as _re
    p = os.path.join(os.path.dirname(DOCS), "data", "README.md")
    lines = open(p, encoding="utf-8").read().split("\n")
    out, touched = [], 0
    for ln in lines:
        m = _re.match(r"\| `([^`]+)` \| ([^|]+) \| ([\d,]+) rows \| (.*)\|\s*$", ln)
        if not m:
            out.append(ln); continue
        rel, _sz, old, desc = m.groups()
        fp = os.path.join(os.path.dirname(DOCS), "data", rel)
        if not os.path.exists(fp):
            out.append(ln); continue
        mb = os.path.getsize(fp) / 1e6
        size = f"{mb:.1f} MB" if mb >= 1 else f"{mb*1000:.0f} KB"
        op = gzip.open if fp.endswith(".gz") else open
        with op(fp, "rt", errors="replace") as fh:
            rows = sum(1 for _ in fh) - 1
        # the description often repeats its own row count; keep the two in step
        desc = desc.replace(old, f"{rows:,}")
        out.append(f"| `{rel}` | {size} | {rows:,} rows | {desc}|")
        touched += 1
    open(p, "w", encoding="utf-8").write("\n".join(out))
    print(f"re-measured {touched} rows of {p}")


def refresh_citation(f, v):
    """Keep CITATION.cff's abstract in step with the atlas it describes."""
    import re as _re
    p = os.path.join(os.path.dirname(DOCS), "CITATION.cff")
    s = open(p, encoding="utf-8").read()
    s = _re.sub(r"spanning [\d,]+ de-duplicated samples across [\d,]+ studies",
                f"spanning {f['samples']:,} de-duplicated samples across "
                f"{f['studies']:,} studies", s)
    s = _re.sub(r"and \d+ named diagnostic entities",
                f"and {f['entities']} named diagnostic entities", s)
    open(p, "w", encoding="utf-8").write(s)
    print(f"refreshed {p}")


def main():
    f = json.load(open(os.path.join(DOCS, "whitepaper_facts.json")))
    v = dict(f)
    t4 = list(csv.DictReader(topen("T4_samples_atomic.tsv"), delimiter="\t"))
    T13 = list(csv.DictReader(topen("T13_incidence_vs_data.tsv"), delimiter="\t"))
    T14 = list(csv.DictReader(topen("T14_ega_sarcoma_datasets.tsv"), delimiter="\t"))
    T20 = list(csv.DictReader(topen("T20_ccdi_studies.tsv"), delimiter="\t"))
    T22b = list(csv.DictReader(topen("T22b_ccdi_file_census.tsv"), delimiter="\t"))
    T24 = list(csv.DictReader(topen("T24_ccdi_entity_counts.tsv"), delimiter="\t"))
    T28 = list(csv.DictReader(topen("T28_stjude_viz_tracks.tsv"), delimiter="\t"))
    epi = [r for r in t4 if r["is_epigenomic"] == "Y" and r["entity_kind"] == "sarcoma"
           and r["is_duplicate"] != "Y"]
    reg = [r for r in epi if r["assay_class"] in REG
           and r["epi_target_norm"] not in ("input/none", "none")]

    v["t4_rows"] = N(len(t4))
    v["samples_n"] = N(f["samples"])
    v["entities_n"] = f["entities"]
    v["models_n"] = N(f["n_models"])
    v["problematic_n"] = f["n_problematic"]
    v["studies_all"] = N(f["studies"] + f["controlled_datasets"])
    v["studies_geo"] = N(f["studies"])
    v["controlled_datasets_n"] = f["controlled_datasets"]
    # Read the record count out of the gap map's own payload rather than recomputing it.
    # Stage 23 adds anchored zero-entity rows and the controlled-access datasets to the
    # GEO grain; a recomputation here was 3,280 against the map's 4,573, and the number a
    # reader can check is the one the map shows them.
    gm = open(os.path.join(DOCS, "index.html"), encoding="utf-8").read()
    tot = json.loads(re.search(r'"totals":(\{[^}]*\})', gm).group(1))
    v["grain"] = N(tot["records"])
    # Residual bins are not a constant: the delta harvest added two. The README said
    # "four" through three runs in which there were five, then six.
    v["bins"] = tot["entities"] - f["entities"] - 1
    v["bins_word"] = {4: "four", 5: "five", 6: "six", 7: "seven"}.get(v["bins"],
                                                                      str(v["bins"]))
    for k, lit in (("samples", f["samples"]), ("series", f["studies"])):
        if tot[k] != lit:
            sys.exit(f"gap map {k}={tot[k]:,} but the tables say {lit:,} -- re-run 23")

    # ---- burden
    ped = sum(F(r, "US_cases_per_year_all_ages") for r in T13 if r["age_class"] == "pediatric")
    v["pct_adult_cases"] = f["pct_cases_adult"]
    v["cases_total_n"] = N(f["cases_total_round"])
    v["pct_reg_pediatric"] = round(100 * f["reg_ped"] / f["reg_rated_total"])
    v["gap_n"] = f["gap"]
    # f["zero_reg_named"] is already prose-cased and already ordered by burden, so the
    # rates pair with it positionally. Re-lowercasing here is what turned BCOR-CCNB3 into
    # "bcor-ccnb3" in the first draft of this file.
    zero = sorted((r for r in T13 if I(r, "regulatory_epigenomic_samples") == 0),
                  key=lambda r: -F(r, "US_cases_per_year_all_ages"))
    costed = list(zip(f["zero_reg_named"],
                      [F(r, "US_cases_per_year_all_ages") for r in zero]))
    v["zero_reg_costed"] = ", ".join(f"{nm} ({N(c)}/yr)" for nm, c in costed[:4] if c)
    v["zero_reg_rest"] = ", ".join(nm for nm, c in costed[4:] if c)
    v["zero_reg_cases_n"] = N(f["zero_reg_cases"])

    # ---- the sharpest single-entity case for model predominance
    per = collections.defaultdict(collections.Counter)
    for r in reg:
        per[r["disease"]]["all"] += 1
        if r["sample_type"] in PDERIVED:
            per[r["disease"]]["pd"] += 1
    name = {r["atlas_disease"]: r["display_name"] for r in T13}
    sharp = max(((c["all"], name.get(d, d)) for d, c in per.items()
                 if d in name and not c["pd"]), default=(0, ""))
    v["sharp_entity"] = sharp[1]
    v["sharp_total"] = N(sharp[0])

    # ---- GSE140686
    g = [r for r in t4 if r["gse"] == "GSE140686"]
    v["g686_n"] = N(f["gse140686_n"])
    v["g686_labels"] = f["gse140686_labels"]
    v["g686_entities"] = len({r["disease"] for r in g if r["disease"] in name})
    v["g686_only"] = f["gse140686_only"]

    # ---- CCDI and GDC
    v["ccdi_studies"] = len(T20)
    v["ccdi_participants"] = N(sum(I(r, "num_of_participants") for r in T20))
    v["ccdi_files"] = N(sum(I(r, "num_of_files") for r in T20))
    cen = {r["library_strategy"]: 0 for r in T22b}
    for r in T22b:
        cen[r["library_strategy"]] += I(r, "files")
    v["ccdi_wgs"] = N(cen.get("WGS", 0)); v["ccdi_wxs"] = N(cen.get("WXS", 0))
    v["ccdi_sarc_participants"] = N(sum(I(r, "ccdi_participants") for r in T24))
    v["ccdi_sarc_samples"] = N(sum(I(r, "ccdi_samples") for r in T24))
    v["ccdi_sarc_methyl"] = N(sum(I(r, "ccdi_participants_with_methylation") for r in T24))
    v["ccdi_sarc_reg"] = N(sum(I(r, "ccdi_regulatory_epigenomic_samples") for r in T24))
    v["ccdi_sarc_entities"] = len(T24)

    # ---- controlled access
    v["ega_datasets"] = N(len(T14))
    v["controlled_total_n"] = N(f["controlled_total"])
    v["controlled_reg_n"] = N(f["controlled_reg"])
    # Every CSTN track is the same material as the EGA deposit, so none is added to any
    # total; the inventory exists to say what the walled resource actually contains.
    v["cstn_tracks"] = N(len(T28))
    v["cstn_reg"] = N(sum(1 for r in T28 if r["assay_family"] == "regulatory"))
    v["cstn_models"] = len({r["model_id"] for r in T28 if r["model_id"]})
    v["cstn_marks"] = len({r["assay"] for r in T28 if r["assay"]})

    # Capitalized variants, so the template never has to do case work inline.
    for k in [k for k in v if k.endswith("_word")]:
        v[k + "_cap"] = str(v[k])[:1].upper() + str(v[k])[1:]
    v["n_corrections_cap"] = f["n_corrections_word"].capitalize()

    tmpl = open(os.path.join(ROOT, "pipeline", "README.tmpl.md"), encoding="utf-8").read()
    try:
        out = tmpl.format(**v)
    except KeyError as e:
        sys.exit(f"README template wants {e} and the tables do not provide it")
    p = os.path.join(ROOT, "README.md")
    open(p, "w", encoding="utf-8").write(rewrap(out))
    print(f"wrote {p} ({len(out.split())} words, {len(v)} values substituted)")
    refresh_data_readme()
    refresh_citation(f, v)


if __name__ == "__main__":
    main()
