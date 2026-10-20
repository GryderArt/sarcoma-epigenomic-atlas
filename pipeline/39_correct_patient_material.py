#!/usr/bin/env python3
"""Keep cell lines, models, other cancers and non-human material out of the patient counts.

The gap map showed 143B -- an osteosarcoma cell line -- with a PATIENT badge. It was one of
several hundred. `sample_type` was decided by the first phrase that matched, and tumour
words were tested before cell-line words, so "highly metastatic osteosarcoma cell line"
read as a metastasis; the catalogued model type could only overrule a record that had
matched nothing at all. Three further routes put non-patient material in the patient
columns, and all three surfaced while auditing the first:

  a model-name collision   "BT12" (breast tumour 12) matched the ATRT line BT-12 and
                           carried breast cancers into ATRT; "Giant cell tumor" is an alias of
                           the GCT line and claimed bone tumour biopsies; reference-RNA pools
                           listing U2OS claimed 78 osteosarcoma biopsies for U-2 OS.
  non-human material       a mouse, zebrafish or canine "primary tumour" is not a patient's;
                           neither is mouse liver, zebrafish embryo or Arabidopsis.
  comparison samples       Wilms tumours inside a CCSK series, breast tumours inside an Ewing
                           series: their own records name a different cancer.

The rules live in `04_classify.patient_material()` and are applied here to the shipped
table, so a fresh classification and this correction cannot drift:

  20a  the record presents the sample as its catalogued model  -> the model's own type
       (cell line; CDX when engrafted; PDX; organoid)
  20b  the model name is a patient label or an incidental match -> model link removed, the
       sample stays patient tissue (an RMS subtype that came only from that model reverts to
       RMS-NOS). Kept as a donor link when the same series profiles the derived model.
  20c  ... and the model was the only link to a sarcoma entity    -> EXCLUDED - not a sarcoma
  20d  non-human host                                           -> mouse_model (rodents) or
                                                                   non-human
  20e  the record names a different cancer and no sarcoma         -> EXCLUDED - not a sarcoma
  20f  the record's own title or source calls it a cell line      -> cell_line

Each change is logged with its basis in T33. Run after 37, before 08; re-run after 32, 06,
07, 30, 36 and 37 whenever the table is rebuilt. Idempotent: a second run changes nothing
and keeps the log.
"""
import csv, gzip, os, sys, shutil, collections, importlib.util
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import resolve, topen, twrite
csv.field_size_limit(10**7)


def _mod(name, fn):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(os.path.dirname(os.path.abspath(__file__)), fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


C4 = _mod("c04", "04_classify.py")
LOG = "T33_patient_material_corrections.tsv"
COLS = ["gsm", "gse", "rule", "sample_type_from", "sample_type_to", "disease_from",
        "disease_to", "model_from", "model_to", "organism", "assay_class", "is_epigenomic",
        "basis", "title"]
EXCLUDED = "EXCLUDED - not a sarcoma"
RMS = {"FP-RMS", "FN-RMS", "RMS-MYOD1"}
PATIENT = C4.PATIENT_ST
MODEL_ST = {"cell_line", "PDX", "organoid", "xenograft(CDX)", "mouse_model"}
REG = {"ChIP-seq", "CUT&RUN", "CUT&Tag", "ChIP-exo", "ChIP-chip", "ATAC-seq", "scATAC-seq",
       "DNase-seq", "FAIRE-seq", "MNase-seq", "Hi-C", "HiChIP", "Micro-C", "Capture-HiC",
       "ChIA-PET", "4C-seq"}


def patient_reg(rows):
    """Regulatory patient-tumour profiles per entity, as the white-paper figures count them."""
    c = collections.Counter()
    for r in rows:
        if r["entity_kind"] == "sarcoma" and r["is_epigenomic"] == "Y" \
                and r["is_duplicate"] != "Y" and r["assay_class"] in REG \
                and r["epi_target_norm"] not in ("input/none", "none") \
                and r["sample_type"] in PATIENT:
            c[r["disease"]] += 1
    return c


def main():
    path = resolve("T4_samples_atomic.tsv")
    rows = list(csv.DictReader(topen("T4_samples_atomic.tsv"), delimiter="\t"))
    hdr = list(rows[0].keys())
    before = patient_reg(rows)

    # a model profiled as a model somewhere in the same series: a patient tumour labelled with
    # its name is the tumour the model was derived from (GSE268193: ARMS1_T beside ARMS1_X)
    derived = collections.defaultdict(set)
    for r in rows:
        if r["model_matched"] and r["sample_type"] in MODEL_ST \
                and C4.judge_model(r)[0] == "model":
            derived[r["gse"]].add(r["model_matched"])

    log = []
    for r in rows:
        if r["entity_kind"] not in ("sarcoma", "benign"):
            continue
        acts = C4.patient_material(r)
        if not acts:
            continue
        frm = (r["sample_type"], r["disease"], r["model_matched"])
        rules, why = [], []
        for rule, ch, basis in acts:
            if rule == "20b" and r["model_matched"] in derived[r["gse"]]:
                continue
            if ch.get("clear_model"):
                m = r["model_matched"]
                if r["donor_key"] in (m, "donor:" + m):
                    r["donor_key"] = ""
                for k in ("model_matched", "model_rrid", "model_subtype", "model_problematic"):
                    r[k] = ""
                if not ch.get("exclude") and r["disease"] in RMS \
                        and m and m in (r.get("rms_call_basis") or ""):
                    r["disease"] = "RMS-NOS"
                    r["subtype_evidence"] = ""
                    r["rms_call_basis"] = (f"reverted to RMS-NOS: the subtype came from {m}, "
                                           f"a model-name collision (correction #20)")
                    basis += f"; the {frm[1]} call rested on {m} and reverts to RMS-NOS"
                if r["relevance"] == "model_confirmed":
                    r["relevance"] = {"sample": "sample_explicit",
                                      "series": "series_explicit"}.get(r["disease_source"],
                                                                       "sarcoma_nos")
            if ch.get("exclude"):
                r["disease"] = EXCLUDED
                r["entity_kind"] = "nonsarcoma"
                r["disease_source"] = f"correction #20: {basis}"
                r["relevance"] = "off_target"
            if ch.get("sample_type"):
                r["sample_type"] = ch["sample_type"]
            rules.append(rule); why.append(basis)
        if not rules:
            continue
        log.append({"gsm": r["gsm"], "gse": r["gse"], "rule": "+".join(rules),
                    "sample_type_from": frm[0], "sample_type_to": r["sample_type"],
                    "disease_from": frm[1], "disease_to": r["disease"],
                    "model_from": frm[2], "model_to": r["model_matched"],
                    "organism": r["organism"], "assay_class": r["assay_class"],
                    "is_epigenomic": r["is_epigenomic"], "basis": " | ".join(why),
                    "title": r["title"][:90]})

    # keep earlier entries whose correction still stands, so a re-run does not erase the log
    state = {r["gsm"]: (r["sample_type"], r["disease"]) for r in rows}
    done = {x["gsm"] for x in log}
    try:
        for x in csv.DictReader(topen(LOG), delimiter="\t"):
            if x["gsm"] not in done and state.get(x["gsm"]) == (x["sample_type_to"], x["disease_to"]):
                log.append(x)
    except FileNotFoundError:
        pass
    fresh = len(done)
    log.sort(key=lambda x: (x["rule"], x["gse"], x["gsm"]))

    tmp = path + ".tmp"
    op = gzip.open if path.endswith(".gz") else open
    with op(tmp, "wt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader(); [w.writerow(r) for r in rows]
    shutil.move(tmp, path)
    with twrite(LOG, gz=False) as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader(); [w.writerow(x) for x in log]

    after = patient_reg(rows)
    by = collections.Counter(x["rule"] for x in log)
    epi = collections.Counter(x["rule"] for x in log if x["is_epigenomic"] == "Y")
    print(f"corrected this run: {fresh} samples; logged in total: {len(log)} "
          f"({sum(epi.values())} epigenomic) -> {LOG}")
    for k in sorted(by):
        print(f"   {k:10s} {by[k]:6d}   epigenomic {epi.get(k, 0):5d}")
    print(f"\nregulatory profiles on patient tumour tissue: {sum(before.values()):,} -> "
          f"{sum(after.values()):,}")
    for d in sorted(set(before) | set(after), key=lambda d: after[d] - before[d]):
        if before[d] != after[d]:
            print(f"   {d:36s} {before[d]:5d} -> {after[d]:5d}")


if __name__ == "__main__":
    main()
