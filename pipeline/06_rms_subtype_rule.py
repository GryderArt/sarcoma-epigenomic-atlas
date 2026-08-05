#!/usr/bin/env python3
"""Re-apply the RMS subtype rule to the v3 atlas.

The user's rule: an RMS sample is FUSION-NEGATIVE when the record establishes that it
carries neither a PAX3/PAX7 fusion nor a mutant MYOD1 -- some FN-RMS carry no RAS
mutation at all, so a RAS call must not be the criterion. GEO almost never writes
"fusion-negative"; it writes the driver ("RMS primary tumor, NRAS"), so within a series
that systematically genotypes its cases, the absence of a named PAX fusion is the
usable signal.

Guard rails:
  * only the RMS-NOS bin is re-derived, so calls already made are never demoted;
  * for MODELS the curated model table decides -- driver-genotype inference is never
    allowed to overrule a known cell line (RH4/RH30/RH41/CW9019 are PAX-fusion lines
    that carry TP53/other mutations and would otherwise be swept into FN-RMS);
  * driver-absence inference applies only to patient-derived material.
"""
import csv, re, collections
import os
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
csv.field_size_limit(10**7)
O = DATA

_STOP = r"(?!foxo1|fkhr|ncoa|maml|ino80|foxo4)"
PAX_FUS = re.compile(
    r"pax\s*[37]\s*(?:::|--|[-/:])\s*" + _STOP + r"[A-Za-z][A-Za-z0-9]{2,8}|"
    r"pax\s*[37]\s*[-/:_ ]*\s*(foxo1|fkhr|ncoa[123]|ino80d|maml3|foxo4)|"
    r"\bp3f\b|\bp7f\b|pax3\s*foxo1|pax7\s*foxo1|"
    r"fusion[- ]positive|fp[- ]rms|\bfprms\b|t\(2;13\)|t\(1;13\)", re.I)
OTHER_FUS = re.compile(r"vgll2|srf\s*[-:]{1,2}\s*ncoa|tead1\s*[-:]{1,2}\s*ncoa2|"
                       r"ncoa2\s*[-:]{1,2}\s*(vgll2|srf|tead1)", re.I)
MYOD1_MUT = re.compile(r"myod1\s*[- ]?(l122r|p\.?l122r|mutant|mutation|mut\b)|"
                       r"mutant\s+myod1|myod1[- ]mutant", re.I)
MYOD1_WT  = re.compile(r"myod1\s*(wild[- ]?type|wt\b)", re.I)
FN_WORDS  = re.compile(r"fusion[- ]negative|\bfn[- ]rms\b|\bfnrms\b|\berms\b|embryonal", re.I)
RMSY = re.compile(r"rhabdomyosarc|\brms\b|\barms\b|\berms\b|myogenic|myoblast", re.I)
# a NAMED driver gene -- generic words like "genotype" or "mutation" are not enough
DRIVER = re.compile(r"\b[nkh]ras\b|\bpik3ca\b|\btp53\b|\bfgfr4\b|\bfbxw7\b|\bctnnb1\b|"
                    r"\bbcor\b|\bmyod1\b|\bmycn\b|\bcdkn2a\b|\bnf1\b|\bpten\b|\bdicer1\b|"
                    r"\btert\b|\bigf1r\b|\barid1a\b|\batrx\b|\bmet\b", re.I)
MODEL_ST = {"cell_line","organoid","xenograft(CDX)","mouse_model"}
PATIENT_ST = {"primary_tumor","metastasis","recurrence","PDX"}

def model_subtype(s):
    if not s: return None
    t = s.strip().lower()
    if t.startswith(("ambiguous","not rms")) or "unspecified" in t or "not specified" in t \
       or "but fusion-posit" in t or "questionable" in t:
        return None
    if "myod1" in t and "wild" not in t: return "RMS-MYOD1"
    if t.startswith("fp-rms") or t.startswith("fp rms"): return "FP-RMS"
    if t.startswith("fn-rms") or t.startswith("assigned fn-rms"): return "FN-RMS"
    if t.startswith("sc-rms-ncoa2"): return "FN-RMS"
    return None

def blob(r):
    """Sample-level text only -- never the series title, which leaks the study's
    headline subtype onto every sample in the series."""
    return " ".join((r["title"], r["source_name"], r["characteristics"]))

def main():
    MSUB = {}
    for m in csv.DictReader(topen("T2_models.tsv"), delimiter="\t"):
        s = model_subtype(m.get("subtype_or_fusion",""))
        if s:
            for nm in [m["model_name"]] + [a.strip() for a in (m.get("aliases") or "").split(";")]:
                if nm: MSUB[nm.strip().lower()] = s

    rows = list(csv.DictReader(topen("T4_samples_atomic.tsv"), delimiter="\t"))
    hdr = list(rows[0].keys())
    if "rms_call_basis" not in hdr: hdr.append("rms_call_basis")

    # which series systematically genotype their RMS patient material?
    gsys = collections.defaultdict(lambda: [0, 0])
    for r in rows:
        if r["sample_type"] in PATIENT_ST and RMSY.search(blob(r)):
            gsys[r["gse"]][1] += 1
            if DRIVER.search(blob(r)) or PAX_FUS.search(blob(r)): gsys[r["gse"]][0] += 1
    SYS = {g for g, (a, b) in gsys.items() if b >= 4 and a / b >= 0.6}

    changed = collections.Counter(); log = []
    for r in rows:
        r.setdefault("rms_call_basis", "")
        if r["disease"] != "RMS-NOS": continue
        b = blob(r)
        if not RMSY.search(b) and not RMSY.search(r["gse_title"]): continue
        ms = MSUB.get((r["model_matched"] or "").strip().lower()) or \
             model_subtype(r.get("model_subtype",""))
        new = basis = None
        if MYOD1_MUT.search(b) and not MYOD1_WT.search(b):
            new, basis = "RMS-MYOD1", "mutant MYOD1 named in sample record"
        elif PAX_FUS.search(b):
            new, basis = "FP-RMS", "PAX3/7 fusion named in sample record"
        elif r["sample_type"] in MODEL_ST and ms:
            new, basis = ms, f"curated subtype of model {r['model_matched']}"
        elif OTHER_FUS.search(b):
            new, basis = "FN-RMS", "non-PAX RMS fusion named (no PAX, no mutant MYOD1)"
        elif FN_WORDS.search(b):
            new, basis = "FN-RMS", "explicit fusion-negative / embryonal wording"
        elif ms:
            new, basis = ms, f"curated subtype of model {r['model_matched']}"
        elif r["sample_type"] in PATIENT_ST and DRIVER.search(b):
            new, basis = "FN-RMS", "driver genotype annotated, no PAX fusion and no mutant MYOD1"
        elif r["sample_type"] in PATIENT_ST and r["gse"] in SYS:
            new, basis = "FN-RMS", ("series systematically reports fusion status; "
                                    "no PAX fusion and no mutant MYOD1 reported for this case")
        if new and new != "RMS-NOS":
            log.append((r["gsm"], r["gse"], "RMS-NOS", new, basis, r["sample_type"],
                        r["title"][:48], r["source_name"][:48]))
            r["disease"] = new; r["rms_call_basis"] = basis
            changed[new] += 1

    # ---- model authority: a curated cell line identity outranks series context.
    # v2 assigned some series-level subtypes to samples that are unambiguously an
    # FN-RMS or FP-RMS line (RD and SMS-CTR sitting in FP-RMS, for instance).
    RMSBINS = {"FP-RMS","FN-RMS","RMS-MYOD1","RMS-NOS"}
    over = collections.Counter()
    for r in rows:
        if r["disease"] not in RMSBINS: continue
        if r["sample_type"] not in MODEL_ST | {"PDX"}: continue
        ms = MSUB.get((r["model_matched"] or "").strip().lower())
        if ms and ms != r["disease"]:
            over[(r["disease"], ms)] += 1
            log.append((r["gsm"], r["gse"], r["disease"], ms,
                        f"model authority: {r['model_matched']} is a curated {ms} line",
                        r["sample_type"], r["title"][:48], r["source_name"][:48]))
            r["disease"] = ms
            r["rms_call_basis"] = f"model authority: curated subtype of {r['model_matched']}"
    print("model-authority overrides:", {f"{a} -> {b}": n for (a, b), n in over.most_common()})

    print("RMS-NOS resolved ->", dict(changed))
    print(f"{len(rows)} rows; {sum(changed.values())} reassigned; "
          f"{len(SYS)} systematically-genotyped series")
    with twrite("T9b_rms_subtype_log.tsv") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gsm","gse","from","to","basis","sample_type","title","source_name"])
        [w.writerow(x) for x in log]
    with twrite("T4_samples_atomic.tsv", gz=True) as f:
        w = csv.DictWriter(f, fieldnames=hdr, delimiter="\t", extrasaction="ignore")
        w.writeheader(); [w.writerow(r) for r in rows]

    print("\ncheck 1 -- Gryder/Yohe primary tumours (GSE83728), user's ground truth:")
    for r in rows:
        if r["gse"] == "GSE83728" and r["sample_type"] == "primary_tumor" \
           and r["epi_target_norm"] == "H3K27ac":
            print(f"   {r['gsm']}  {r['disease']:10s}  {r['source_name'][:44]}")
    print("\ncheck 2 -- no PAX-fusion cell line may sit in FN-RMS:")
    FP = {"rh4","rh30","rh41","rh28","rh3","rh5","rh10","rh65","kfr","cw9019","tc212",
          "arms1","rmzrc2","scmcrm2","garvin1","humems","nrs1"}
    bad = collections.Counter(r["model_matched"] for r in rows
                              if r["disease"] == "FN-RMS"
                              and re.sub(r"[^a-z0-9]","", (r["model_matched"] or "").lower()) in FP)
    print("   violations:", dict(bad) if bad else "none")
    print("\ncheck 3 -- no FN cell line may sit in FP-RMS:")
    FN = {"rd","smsctr","rh36","rmsym","rms559","jr1","ruch2","ruch3","ttc442","ttc516"}
    bad2 = collections.Counter(r["model_matched"] for r in rows
                               if r["disease"] == "FP-RMS"
                               and re.sub(r"[^a-z0-9]","", (r["model_matched"] or "").lower()) in FN)
    print("   violations:", dict(bad2) if bad2 else "none")

if __name__ == "__main__":
    main()
