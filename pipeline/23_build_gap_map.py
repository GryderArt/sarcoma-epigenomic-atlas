#!/usr/bin/env python3
"""Build the interactive all-ages Sarcoma Epigenomic Gap Map (self-contained HTML).

Payload is dictionary-encoded columnar arrays so the whole atlas -- every
model x assay x target x series combination -- travels inside one file with no
network dependency."""
import csv, json, collections, os, re
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
csv.field_size_limit(10**7)
O = DATA
OUT = os.path.join(DOCS, "index.html")

EPI = {"ChIP-seq","CUT&RUN","CUT&Tag","ChIP-exo","ChIP-chip","ATAC-seq","scATAC-seq",
       "DNase-seq","FAIRE-seq","MNase-seq","Hi-C","HiChIP","Micro-C","Capture-HiC",
       "ChIA-PET","4C-seq","Repli-seq","WGBS","RRBS","Methyl-array","MeDIP/hMeDIP",
       "Bisulfite-PCR"}
ACC  = {"ATAC-seq","scATAC-seq","DNase-seq","FAIRE-seq","MNase-seq"}
D3   = {"Hi-C","HiChIP","Micro-C","Capture-HiC","ChIA-PET","4C-seq"}
DNA  = {"WGBS","RRBS","Methyl-array","MeDIP/hMeDIP","Bisulfite-PCR"}
CHIP = {"ChIP-seq","CUT&RUN","CUT&Tag","ChIP-exo","ChIP-chip"}
CTRL_TARGET = {"input/none","none"}
ACTIVE = {"H3K27ac","H3K4me1","H3K4me2","H3K9ac","H3K9-14ac","H4ac","H3K18ac","H2BK20ac"}
REPRESS = {"H3K27me3","H3K9me3","H3K9me2","H2AK119ub","H4K20me3"}
PROMOTER = {"H3K4me3","H3K36me3","H3K79me2"}
PATIENT = {"primary_tumor","metastasis","recurrence"}
PDER = PATIENT | {"PDX","organoid"}
MODELS = {"cell_line","organoid","xenograft(CDX)","mouse_model","PDX"}

def lane(assay, target):
    if assay in DNA: return "DNA methylation"
    if assay in ACC: return "Accessibility"
    if assay in D3:  return "3D genome"
    if assay in CHIP:
        if target in CTRL_TARGET: return "Input / control"
        if target == "H3K27ac":   return "H3K27ac"
        if target in ACTIVE:      return "Other active marks"
        if target in REPRESS:     return "Repressive marks"
        if target in PROMOTER:    return "Promoter / body marks"
        return "TF / cofactor"
    return "Other"

LANES = ["H3K27ac","Other active marks","Promoter / body marks","Repressive marks",
         "TF / cofactor","Accessibility","3D genome","DNA methylation","Input / control"]

# Not every bin in the atlas is a disease. Three kinds:
#   named    -- a diagnostic entity in the WHO-2020 sense (RMS and liposarcoma split to
#               molecular subtype, since that is what the paper is about). Denominators use
#               ONLY these.
#   residual -- a bin for samples whose entity could not be resolved from the deposit
#               ("Sarcoma NOS", "RMS-NOS") or one-off cases with nowhere else to sit.
#               Counting these as entities would inflate every denominator.
#   control  -- non-neoplastic reference tissue that is not a disease at all.
RESIDUAL = {"Sarcoma NOS", "RMS-NOS", "Other rare bone tumor", "Myoepithelial tumor",
            "Sarcoma, unspecified / mixed cohort", "RMS (any subtype)",
            "Liposarcoma (any subtype)"}
def kindof(name, anchored):
    n = name.lower()
    if "control" in n or "non-malignant" in n or "normal" in n: return "control"
    if name in RESIDUAL or n.endswith(" nos") or "unspecified" in n: return "residual"
    return "named" if anchored else "residual"
REG_LANES = set(LANES[:7])

def main():
    rows = list(csv.DictReader(topen("T4_samples_atomic.tsv"), delimiter="\t"))
    epi = [r for r in rows if r["assay_class"] in EPI and r["is_duplicate"] != "Y"
           and r["entity_kind"] == "sarcoma"]
    T13 = {r["atlas_disease"]: r for r in
           csv.DictReader(topen("T13_incidence_vs_data.tsv"), delimiter="\t")}
    G3 = {r["disease"]: r for r in
          csv.DictReader(topen("T3_entity_counts.tsv"), delimiter="\t")}
    CTRL = list(csv.DictReader(topen("T15_controlled_access.tsv"),
                               delimiter="\t"))
    K686 = {r["gsm"]: r for r in csv.DictReader(topen("GSE140686_recovered_diagnoses.tsv"),
                                                delimiter="\t")}

    # ---- aggregate to the "H3K27ac ChIP-seq for RH4" grain
    agg = collections.OrderedDict()
    for r in epi:
        if r["gse"] == "GSE140686":
            k = K686.get(r["gsm"], {})
            name = "DKFZ classifier — " + (k.get("institutional_diagnosis") or r["disease"])
        elif r["model_matched"]:
            name = r["model_matched"]
        elif r["sample_type"] in PATIENT:
            name = (r["source_name"][:44] or r["title"][:44] or "patient tumor")
        else:
            name = (r["source_name"][:44] or r["title"][:44] or "unnamed")
        key = (r["disease"], r["sample_type"], name, r["assay_class"],
               r["epi_target_norm"] or "unspecified", r["gse"])
        a = agg.setdefault(key, {"n": 0, "r": r})
        a["n"] += 1

    D = {"disease": [], "stype": [], "model": [], "assay": [], "target": [], "gse": [],
         "lane": [], "title": []}
    idx = {k: {} for k in D}
    def enc(field, v):
        m = idx[field]
        if v not in m:
            m[v] = len(D[field]); D[field].append(v)
        return m[v]

    recs = []
    for (dis, st, name, assay, tgt, gse), a in agg.items():
        r = a["r"]
        recs.append([enc("disease", dis), enc("stype", st), enc("model", name),
                     enc("assay", assay), enc("target", tgt), enc("gse", gse),
                     a["n"], enc("lane", lane(assay, tgt)),
                     int(re.split(r"[;,\s]+", (r["gse_pubmed"] or "0").strip())[0] or 0),
                     int(re.sub(r"\D", "", (r["gse_date"] or "")[:4]) or 0),
                     enc("title", r["gse_title"][:120])])

    # ---- per-entity summary
    byd = collections.defaultdict(list)
    for i, rec in enumerate(recs): byd[D["disease"][rec[0]]].append(i)
    creg, cdna, cds = collections.Counter(), collections.Counter(), collections.defaultdict(list)
    CMAP = {"UPS / undifferentiated sarcoma": "UPS/MFH"}
    UNRES = {"RMS (any subtype)", "Liposarcoma (any subtype)",
             "Sarcoma, unspecified / mixed cohort"}
    mixed = collections.defaultdict(list)
    # Stage 12 rewrites T15 from scratch and stage 36 adds these columns afterwards, so
    # running them out of order silently drops every cohort attribution. Say so loudly.
    if CTRL and "contains_entities" not in CTRL[0]:
        print("  !! T15 has no cohort attribution -- run 36_attribute_mixed_cohorts.py "
              "after 12 and before this stage, or 430 controlled regulatory samples "
              "will reach no entity page")
    for c in CTRL:
        n = int(c["n_samples"] or 0)
        e = CMAP.get(c["entity"], c["entity"])
        cds[e].append({"acc": c["accession"], "n": n, "fam": c["assay_family"],
                       "title": c["title"][:130], "dac": c["governing_body"][:70],
                       "src": c["source"], "url": c["url"]})
        if e in UNRES:
            # A cohort deposited above the subtype still bears on the subtypes inside it.
            # Stage 36 records which, and on what evidence; it reaches every entity it
            # names but is added to none of their totals, because the same samples cannot
            # be counted once under each. Without this the St Jude RMS ChIP-seq -- 400
            # samples -- appeared on no entity page at all.
            models = dict(kv.split(":") for kv in (c.get("contains_models") or "").split(";")
                          if ":" in kv)
            for t in (c.get("contains_entities") or "").split(";"):
                if not t:
                    continue
                mixed[t].append({"acc": c["accession"], "n": n, "fam": c["assay_family"],
                                 "title": c["title"][:130], "dac": c["governing_body"][:70],
                                 "src": c["source"], "url": c["url"], "bin": e,
                                 "basis": c.get("contains_basis", ""),
                                 "models": int(models.get(t, 0))})
            continue
        (creg if c["assay_family"] == "regulatory" else cdna)[e] += n

    ents = []
    for dis, ids in sorted(byd.items(), key=lambda kv: -sum(recs[i][6] for i in kv[1])):
        t = T13.get(dis, {}); g = G3.get(dis, {})
        lanes = collections.Counter()
        for i in ids: lanes[D["lane"][recs[i][7]]] += recs[i][6]
        reg = sum(v for k, v in lanes.items() if k in REG_LANES)
        pd_reg = sum(recs[i][6] for i in ids
                     if D["stype"][recs[i][1]] in PDER and D["lane"][recs[i][7]] in REG_LANES)
        pt_reg = sum(recs[i][6] for i in ids
                     if D["stype"][recs[i][1]] in PATIENT and D["lane"][recs[i][7]] in REG_LANES)
        ents.append({
            "d": dis,
            "name": t.get("display_name") or dis,
            "age": t.get("age_class") or g.get("age_class") or "",
            "cases": int(t.get("US_cases_per_year_all_ages") or 0),
            "lo": int(t.get("US_low") or 0), "hi": int(t.get("US_high") or 0),
            "basis": t.get("incidence_basis", ""),
            "caveat": t.get("verification_caveat", ""),
            "n": sum(recs[i][6] for i in ids),
            "reg": reg, "pd_reg": pd_reg, "pt_reg": pt_reg,
            "lanes": dict(lanes),
            "creg": creg.get(dis, 0), "cdna": cdna.get(dis, 0),
            "cds": cds.get(dis, []), "mix": mixed.get(dis, []),
            "anch": 1 if dis in T13 else 0,
            "kind": kindof(t.get("display_name") or dis, dis in T13),
            "models": int(g.get("model_units") or 0), "pdx": int(g.get("pdx_units") or 0),
            "series": len({D["gse"][recs[i][5]] for i in ids}),
            "nrec": len(ids)})
    # anchored entities with no records at all -- the most important gaps, so they must
    # appear as explicit zero rows rather than vanish
    have = {e["d"] for e in ents}
    for dis, t13 in T13.items():
        if dis in have: continue
        ents.append({"d": dis, "name": t13.get("display_name") or dis,
                     "age": t13.get("age_class",""), "anch": 1,
                     "cases": int(t13.get("US_cases_per_year_all_ages") or 0),
                     "lo": int(t13.get("US_low") or 0), "hi": int(t13.get("US_high") or 0),
                     "basis": t13.get("incidence_basis",""),
                     "caveat": t13.get("verification_caveat",""),
                     "kind": kindof(t13.get("display_name") or dis, True),
                     "n": 0, "reg": 0, "pd_reg": 0, "pt_reg": 0, "lanes": {},
                     "creg": creg.get(dis,0), "cdna": cdna.get(dis,0), "cds": cds.get(dis,[]),
                     "mix": mixed.get(dis,[]),
                     "models": 0, "pdx": 0, "series": 0, "nrec": 0})

    # entities that exist only as controlled-access holdings
    for e, ds in cds.items():
        if e in UNRES or any(x["d"] == e for x in ents): continue
        ents.append({"d": e, "name": e, "age": "", "anch": 0, "kind": kindof(e, False),
                     "cases": 0, "lo": 0, "hi": 0,
                     "basis": "", "caveat": "", "n": 0, "reg": 0, "pd_reg": 0, "pt_reg": 0,
                     "lanes": {}, "creg": creg.get(e, 0), "cdna": cdna.get(e, 0),
                     "cds": ds, "mix": mixed.get(e,[]), "models": 0, "pdx": 0, "series": 0, "nrec": 0})
    unres = [{"entity": c["entity"], "acc": c["accession"], "n": int(c["n_samples"] or 0),
              "fam": c["assay_family"], "title": c["title"][:130],
              "dac": c["governing_body"][:70]}
             for c in CTRL if c["entity"] in UNRES]

    payload = {"dict": D, "recs": recs, "ents": ents, "lanes": LANES,
               "reglanes": sorted(REG_LANES), "unres": unres,
               "built": "atlas v3 (all ages)",
               "kinds": {"named": sum(1 for e in ents if e["kind"]=="named"),
                         "residual": sum(1 for e in ents if e["kind"]=="residual"),
                         "control": sum(1 for e in ents if e["kind"]=="control")},
               "totals": {"samples": sum(r[6] for r in recs), "records": len(recs),
                          "series": len(D["gse"]), "entities": len(ents),
                          "ctrl_reg": sum(int(c["n_samples"] or 0) for c in CTRL
                                          if c["assay_family"] == "regulatory"),
                          "ctrl_dna": sum(int(c["n_samples"] or 0) for c in CTRL
                                          if c["assay_family"] == "DNA methylation")}}
    js = json.dumps(payload, separators=(",", ":"))
    print(f"payload {len(js)/1e6:.2f} MB · {len(recs):,} records · {len(ents)} entities")
    tmpl = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gapmap_template.html")).read()
    open(OUT, "w").write(tmpl.replace("/*__PAYLOAD__*/null", js))
    print(f"wrote {OUT} ({os.path.getsize(OUT)/1e6:.2f} MB)")

if __name__ == "__main__":
    main()
