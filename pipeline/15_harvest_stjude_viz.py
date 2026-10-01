#!/usr/bin/env python3
"""Enumerate the St Jude CSTN epigenetic browsers at sample x assay resolution.

The CSTN "epigenetic landscape" pages on viz.stjude.cloud are ProteinPaint embeds. The
page ships its whole track manifest inline -- every sample, every assay, every underlying
file -- so the browser contents can be enumerated exactly, which the CSTN portal's own
Shiny interface does not allow.

What this is for: EGA describes these deposits as, for example, "ChIP-Seq files for RMS,
242 samples". That is opaque. The manifest resolves the same material to the atlas's own
grain -- H3K27ac for SJRHB013758_X1 -- so the gap map can say which marks exist on which
model rather than only how many files there are.

ENTITY CALLS ARE CURATED, NOT INFERRED. An earlier version of this stage read the subtype
off the label suffix -- (ERMS), (ARMS), (SCLEROS). The St Jude data administrator then
supplied the definitive per-sample list by email, which is transcribed verbatim in
ADMIN_EMAIL below and is the primary source for every call. `T17_stjude_opdx_models.tsv`
(the CSTN portal's own model table, harvested independently) is read as a second,
corroborating source and supplies the fusion partner, which the email does not give.

Where the two sources disagree, or where a sample appears in one and not the other, the
discrepancy is reported and left standing -- see the report at the end of main(). Three
are open at the time of writing and are listed in the module docstring of the report.

ACCESS: the tracks are viewable in the browser; the underlying bigWig and bed files sit
behind Cloudflare bot protection and, for real access, the St Jude / EGA data access
agreement. They are recorded as "formal request" and are NOT added to any sample total,
because the same material is already counted once under EGA.
"""
import urllib.request, html, re, json, csv, os, collections

from _paths import DATA, twrite, topen

UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                    "(KHTML, like Gecko) Chrome/120 Safari/537.36"}

PAGES = {
    "RHB": ("Epigenetic landscape of rhabdomyosarcoma subtypes",
            "https://viz.stjude.cloud/st-jude-childrens-research-hospital/visualization/"
            "epigenetic-landscape-of-rhabdomyosarcoma-subtypes~67"),
    "SJOS2018": ("Epigenetic landscape of xenografts from osteosarcoma and rare tumors",
                 "https://viz.stjude.cloud/st-jude-childrens-research-hospital/visualization/"
                 "epigenetic-landscape-of-xenografts-from-osteosarcoma-and-rare-tumors~69"),
}

# assay -> atlas assay class / target
REG_ASSAYS = {"BRD4", "H3K27Ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3",
              "H3K9-14Ac", "H3K9me3", "RNAPolII", "CTCF", "INPUT"}
DERIVED = {"SE", "SuperEnhancer(SE)", "SE.noK4me3"}      # called peaks, not raw assays
DNAME_ASSAYS = {"WGBS"}

# --------------------------------------------------------------------------------------
# GROUND TRUTH -- transcribed verbatim from the St Jude data administrator's email.
# Heading -> the samples the administrator placed under it, plus any note they attached.
# Nothing here is inferred. The atlas entity each heading maps to is in ADMIN_TO_ATLAS.
# --------------------------------------------------------------------------------------
ADMIN_EMAIL = {
    "Rhabdomyosarcoma / Fusion negative": [
        "SJRHB013758_X1", "SJRHB013758_X2", "SJRHB000026_X1", "SJRHB000026_X2",
        "SJRHB010927_X1", "SJRHB012405_X1", "SJRHB010928_X1",
        "SJRHB011_X", "SJRHB011_Y", "SJRHB012_Y", "SJRHB012_Z"],
    "Rhabdomyosarcoma / Fusion positive": [
        "SJRHB010468_X1", "SJRHB013757_X1", "SJRHB013759_X1"],
    "Sclerosal RMS": ["SJRHB015720_X1"],
    "Osteosarcoma": ["SJOS016015_X1", "SJOS010930_X1", "SJOS001105_X1",
                     "SJOS001108_X1", "SJOS001107_X3"],
    "Ewings sarcoma": ["SJEWS001321_X1"],
    "Liposarcoma": ["SJLPS014753_X1"],
    "High grade sarcoma": ["SJHGS015726_X2"],
}
# the one availability statement in the email that is not derivable from the page
ADMIN_NOTES = {"SJOS010930_X1": "removed from CSTN and not available"}

# administrator's wording -> (atlas entity, evidence grade for the subtype call)
#
# "Sclerosal RMS" is spindle cell/sclerosing RMS. The atlas bin for it is RMS-MYOD1,
# because that is the genotype the WHO entity is defined by -- but neither the email nor
# the CSTN portal states MYOD1 status for this model, so the evidence is graded as
# histology, not genotype. The portal reports PIK3CA, which co-occurs in roughly a third
# of MYOD1 L122R cases (PMID 24793135) and is suggestive, not confirmatory.
#
# "High grade sarcoma" is deliberately NOT mapped to UPS/MFH. Undifferentiated pleomorphic
# sarcoma is a specific diagnosis; "high grade sarcoma" is the absence of one. It goes to
# the Sarcoma NOS residual bin, which is excluded from every "of 45 entities" denominator.
ADMIN_TO_ATLAS = {
    "Rhabdomyosarcoma / Fusion negative": ("FN-RMS",         "institutional_diagnosis"),
    "Rhabdomyosarcoma / Fusion positive": ("FP-RMS",         "institutional_diagnosis"),
    "Sclerosal RMS":                      ("RMS-MYOD1",      "histology_spindle_sclerosing"),
    "Osteosarcoma":                       ("Osteosarcoma",   "institutional_diagnosis"),
    "Ewings sarcoma":                     ("Ewing",          "institutional_diagnosis"),
    "Liposarcoma":                        ("Liposarcoma-NOS", "institutional_diagnosis"),
    "High grade sarcoma":                 ("Sarcoma NOS",    "institutional_diagnosis"),
}

# Normal and reference material. The administrator's list covers tumor models only;
# these four are read off the browser, where they are unambiguous.
NORMAL = re.compile(r"^SJNORM", re.I)

MODEL_ID = re.compile(r"^(SJ[A-Z]+\d+_[XYZ]\d*)")
LABEL_SUBTYPE = re.compile(r"\((ERMS|ARMS|SCLEROS)\)", re.I)
LABEL_TO_ATLAS = {"ERMS": "FN-RMS", "ARMS": "FP-RMS", "SCLEROS": "RMS-MYOD1"}


def model_id(label):
    """SJRHB013758_X1_Prim_(ERMS) -> SJRHB013758_X1. Normals keep their whole label."""
    m = MODEL_ID.match(label)
    return m.group(1) if m else label


def build_curation():
    """model id -> curated record, from the administrator's email."""
    cur = {}
    for heading, ids in ADMIN_EMAIL.items():
        ent, ev = ADMIN_TO_ATLAS[heading]
        for i in ids:
            cur[i] = {
                "atlas_entity": ent,
                "admin_label": heading,
                "subtype_evidence": ev,
                "entity_call_basis": "St Jude data administrator (email)",
                "availability": ("withdrawn - " + ADMIN_NOTES[i] if i in ADMIN_NOTES
                                 else "in CSTN"),
            }
    return cur


def read_portal():
    """The CSTN portal's own model table, harvested independently in stage 17.

    Used only to corroborate and to supply the fusion partner. It never overrides the
    administrator; a conflict is reported instead.
    """
    out = {}
    try:
        for r in csv.DictReader(topen("T17_stjude_opdx_models.tsv"), delimiter="\t"):
            out[r["model_name"]] = r
    except Exception as e:
        print(f"  !! could not read T17 for corroboration: {e}")
    return out


PATIENT = re.compile(r"^(SJ[A-Z]+\d+)_")


def portal_record(mid, portal):
    """The portal row for this model, or the same patient at another passage.

    SJRHB013757_X1 is on the browser page; the portal carries only SJRHB013757_X2. A
    fusion is a property of the patient's tumor, not of the passage, so the sibling row
    is allowed to supply it -- flagged, so no reader mistakes it for an exact match.
    """
    if mid in portal:
        return portal[mid], ""
    m = PATIENT.match(mid)
    if m:
        sib = sorted(k for k in portal if k.startswith(m.group(1) + "_"))
        if sib:
            return portal[sib[0]], sib[0]
    return None, ""


def portal_fusion(rec):
    """Pull the fusion partner out of the portal's free-text subtype field."""
    s = (rec.get("subtype") or "") if rec else ""
    m = re.search(r"(PAX[37])[- ]?FOXO1", s, re.I)
    if m:
        return f"{m.group(1).upper()}::FOXO1"
    if re.search(r"fusion negative", s, re.I):
        return "none (fusion negative)"
    return ""


def fetch(url):
    r = urllib.request.Request(url, headers=UA)
    return urllib.request.urlopen(r, timeout=120).read().decode("utf-8", "replace")


def extract_config(body):
    i = body.find("runproteinpaint(")
    if i < 0:
        raise RuntimeError("no runproteinpaint( call on the page")
    seg = html.unescape(body[i + len("runproteinpaint("):])
    start = seg.find("{")
    depth, end = 0, None
    for j, ch in enumerate(seg[start:], start):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                end = j + 1
                break
    raw = seg[start:end]
    # the holder is a live DOM handle, not JSON
    raw = re.sub(r'"holder"\s*:\s*[A-Za-z0-9_\[\]\(\)\.]+\s*,', '"holder": null,', raw)
    return json.loads(raw)


def main():
    curation = build_curation()
    portal = read_portal()
    rows, uncurated = [], collections.OrderedDict()

    for study, (title, url) in PAGES.items():
        cfg = extract_config(fetch(url))
        sv = cfg["studyview"]
        for assay in sv.get("assays", []):
            blk = sv.get(assay)
            if not isinstance(blk, dict):
                continue
            for sample, inner in blk.items():
                if sample == "config" or not isinstance(inner, dict):
                    continue
                mid = model_id(sample)
                c = curation.get(mid)
                if c is None and NORMAL.match(sample):
                    c = {"atlas_entity": "normal/reference",
                         "admin_label": "",
                         "subtype_evidence": "",
                         "entity_call_basis": "browser label (normal/reference material)",
                         "availability": "in CSTN"}
                if c is None:
                    # not in the administrator's list -- fall back to the browser label
                    sfx = LABEL_SUBTYPE.search(sample)
                    c = {"atlas_entity": (LABEL_TO_ATLAS[sfx.group(1).upper()] if sfx
                                          else "Sarcoma NOS"),
                         "admin_label": "",
                         "subtype_evidence": "histology" if sfx else "",
                         "entity_call_basis": "browser label - NOT in the administrator's "
                                              "list (see discrepancy report)",
                         "availability": "in CSTN"}
                    uncurated[mid] = sample

                prec, via = portal_record(mid, portal)
                for t in inner.values():
                    if not isinstance(t, dict):
                        continue
                    fam = ("regulatory" if assay in REG_ASSAYS else
                           "DNA methylation" if assay in DNAME_ASSAYS else
                           "derived (called peaks)")
                    rows.append({
                        "study": study, "viz_title": title,
                        "sample_label": sample,
                        "model_id": mid,
                        "sample_id": t.get("sampleID", ""),
                        "atlas_entity": c["atlas_entity"],
                        "admin_label": c["admin_label"],
                        "entity_call_basis": c["entity_call_basis"],
                        "subtype_evidence": c["subtype_evidence"],
                        "fusion_or_driver": portal_fusion(prec),
                        "fusion_source": ("CSTN portal, exact model" if prec and not via
                                          else f"CSTN portal, sibling passage {via} "
                                               f"(patient-level)" if via else ""),
                        "portal_diagnosis": (prec or {}).get("diagnosis", ""),
                        "portal_subtype": (prec or {}).get("subtype", ""),
                        "availability": c["availability"],
                        "assay": assay, "assay_family": fam,
                        "is_control": "Y" if assay == "INPUT" else "",
                        "genome": sv.get("genome", ""),
                        "track_type": t.get("type", ""),
                        "file": t.get("file", "") or t.get("url", ""),
                        "access": "FORMAL REQUEST (St Jude CSTN / EGA data access agreement)",
                        "counted_in_totals": "no - same material as the EGA deposit",
                        "viz_url": url,
                    })

    with twrite("T28_stjude_viz_tracks.tsv", gz=False) as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); [w.writerow(r) for r in rows]

    # ---------------------------------------------------------------- summary
    print(f"{len(rows)} tracks across {len({r['sample_label'] for r in rows})} samples")
    for study in PAGES:
        g = [r for r in rows if r["study"] == study]
        print(f"  {study:9s} {len(g):4d} tracks · {len({r['sample_label'] for r in g})} samples "
              f"· {len({r['assay'] for r in g})} assays")
    print("\nassay family:", dict(collections.Counter(r["assay_family"] for r in rows)))
    print("\nby atlas entity (samples), with the basis of the call:")
    ent = collections.defaultdict(set); basis = collections.defaultdict(set)
    for r in rows:
        ent[r["atlas_entity"]].add(r["model_id"])
        basis[r["atlas_entity"]].add(r["entity_call_basis"].split(" - ")[0])
    for e, s in sorted(ent.items(), key=lambda kv: -len(kv[1])):
        print(f"   {len(s):3d}  {e:16s} {'; '.join(sorted(basis[e]))}")
    k27 = {r["model_id"] for r in rows if r["assay"] == "H3K27Ac"}
    print(f"\nmodels with H3K27ac: {len(k27)}")

    # ---------------------------------------------------------------- discrepancies
    print("\n" + "=" * 78)
    print("DISCREPANCY REPORT -- reported, not resolved")
    print("=" * 78)
    seen = {r["model_id"] for r in rows}
    n = 0

    for mid, label in uncurated.items():
        n += 1
        p = portal.get(mid)
        print(f"\n[{n}] {mid} is in the browser but NOT in the administrator's email.")
        print(f"    browser label : {label}")
        if p:
            print(f"    CSTN portal   : {p.get('diagnosis','')} / {p.get('subtype','')}")
            if "ChIP-seq" in (p.get("assays_available") or ""):
                acc = re.search(r"ChIP-seq\[([^\]]+)\]", p["assays_available"])
                print(f"    portal lists ChIP-seq for it{' (' + acc.group(1) + ')' if acc else ''}")
        print("    -> kept, entity taken from the browser label and flagged as such.")

    for mid, note in ADMIN_NOTES.items():
        n += 1
        print(f"\n[{n}] {mid}: administrator says \"{note}\".")
        print(f"    Its {sum(1 for r in rows if r['model_id']==mid)} tracks are still on "
              f"the browser page; they are kept and marked withdrawn.")
        p = portal.get(mid)
        if p and "not returned by live CSTN portal API" in (p.get("notes") or ""):
            print("    Corroborated: the live CSTN portal API does not return this model "
                  "either (T17).")

    for mid in sorted(seen):
        c, (p, via) = curation.get(mid), portal_record(mid, portal)
        if via:
            p = None
        if not c:
            continue
        if not p:
            alt = sorted(k for k in portal if k.split("_")[0] == mid.split("_")[0])
            if alt:
                n += 1
                print(f"\n[{n}] {mid} is not in the CSTN portal table; the same patient "
                      f"appears at a different passage: {', '.join(alt)}.")
                print(f"    -> entity kept from the administrator; fusion taken from "
                      f"{alt[0]}, which is patient-level.")
            continue
        pf = portal_fusion(p)
        neg = pf.startswith("none") or not pf
        if c["atlas_entity"] == "FP-RMS" and neg:
            n += 1; print(f"\n[{n}] CONFLICT {mid}: administrator says fusion positive, "
                          f"portal says '{p.get('subtype','')}'.")
        if c["atlas_entity"] == "FN-RMS" and pf.endswith("FOXO1"):
            n += 1; print(f"\n[{n}] CONFLICT {mid}: administrator says fusion negative, "
                          f"portal says '{p.get('subtype','')}'.")
    if not n:
        print("\nnone.")

    # fusion partners the portal adds on top of the administrator's binary call
    fus = {r["model_id"]: r["fusion_or_driver"] for r in rows
           if r["fusion_or_driver"].endswith("::FOXO1")}
    if fus:
        print("\nFusion partner recovered from the CSTN portal (the email gives only "
              "positive/negative):")
        for m, f_ in sorted(fus.items()):
            print(f"    {m:20s} {f_}")


if __name__ == "__main__":
    main()
