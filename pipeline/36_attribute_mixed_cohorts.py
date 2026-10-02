#!/usr/bin/env python3
"""Attribute controlled-access cohorts that are deposited above the subtype.

A PhD student looked up FN-RMS in the gap map and could not find the St Jude PDX
ChIP-seq. It was in the atlas the whole time -- EGAD00001004312 and EGAD00001006398,
400 samples -- filed under the entity "RMS (any subtype)", which is an unresolved bin
rather than a named entity, so it attached to no entity page at all. It is counted in the
atlas-wide controlled-access total and invisible everywhere a reader would look for it.

That is not one stranded dataset. 430 of the 905 controlled regulatory samples -- 48% --
sit in cohort-level bins for the same reason: the depositor stated "rhabdomyosarcoma" or
"sarcoma" and the atlas, correctly, refused to guess a subtype.

Refusing to guess was right. Refusing to *say* was not. This stage adds what is known
about each such cohort, at two clearly separated levels of evidence:

  verified   the St Jude data administrator's per-model list (T28, correction #11) names
             the exact models in the deposit, so the named entities it contains are
             established fact, and so is the model count per entity.

  parent     the deposit states only the parent diagnosis. "Rhabdomyosarcoma" is by
             definition one or more of the three RMS entities, so the cohort necessarily
             bears on them, but which samples belong to which is not recoverable.

Neither tier is added to any entity's totals: the same 400 samples cannot be counted once
under FN-RMS and again under FP-RMS. They are listed, with their size and their basis, so
a reader can see what exists and request it. The headline numbers do not move.

Run after 12, before 23.
"""
import csv, os, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import DATA, topen, twrite, resolve
csv.field_size_limit(10**7)

# Unresolved bins in T15, and the named entities each one's parent diagnosis covers.
PARENT = {
    "RMS (any subtype)": ["FP-RMS", "FN-RMS", "RMS-MYOD1"],
    "Liposarcoma (any subtype)": ["Liposarcoma-WD", "Liposarcoma-dediff", "Liposarcoma-NOS",
                                  "Liposarcoma-myxoid", "Liposarcoma-pleomorphic"],
}
# "Sarcoma, unspecified / mixed cohort" is deliberately absent: it spans the whole
# taxonomy, so attaching it to all 45 entities would be noise rather than information.

NEW = ["contains_entities", "contains_basis", "contains_models"]


def main():
    T15 = list(csv.DictReader(topen("T15_controlled_access.tsv"), delimiter="\t"))
    T28 = list(csv.DictReader(topen("T28_stjude_viz_tracks.tsv"), delimiter="\t"))

    # ---- tier 1: the administrator's per-model list, by entity
    admin = collections.defaultdict(set)
    for r in T28:
        e, m = r.get("atlas_entity", ""), r.get("model_id", "")
        if e and m and e not in ("normal/reference", "Sarcoma NOS") \
                and r.get("entity_call_basis", "").startswith("St Jude data administrator"):
            admin[e].add(m)
    print("administrator-verified models by entity:")
    for e, ms in sorted(admin.items(), key=lambda kv: -len(kv[1])):
        print(f"  {e:22s} {len(ms):3d} models")

    # The deposits that material sits in. Keyed on accession so the link is explicit
    # rather than inferred from a title match.
    CSTN_CHIP = {"EGAD00001004312", "EGAD00001006398", "EGAD00001003432"}

    hdr = list(T15[0].keys())
    for c in NEW:
        if c not in hdr:
            hdr.append(c)

    n_ver = n_par = 0
    for r in T15:
        for c in NEW:
            r.setdefault(c, "")
        bin_ = r["entity"]
        if bin_ not in PARENT:
            continue
        covered = PARENT[bin_]
        if r["accession"] in CSTN_CHIP:
            # restrict to entities the administrator actually placed in this material
            ents = [e for e in covered if e in admin]
            if ents:
                r["contains_entities"] = ";".join(ents)
                r["contains_basis"] = "verified"
                r["contains_models"] = ";".join(f"{e}:{len(admin[e])}" for e in ents)
                n_ver += 1
                continue
        r["contains_entities"] = ";".join(covered)
        r["contains_basis"] = "parent"
        n_par += 1

    out = resolve("T15_controlled_access.tsv")
    with twrite(out, gz=False) as f:
        w = csv.DictWriter(f, fieldnames=hdr, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader(); [w.writerow(r) for r in T15]
    print(f"\n{n_ver} cohorts attributed from the administrator's model list, "
          f"{n_par} from the parent diagnosis alone")

    # ---- report: what this makes reachable, and what is still stranded
    rows, still = [], collections.Counter()
    for r in T15:
        n = int(float(r["n_samples"] or 0))
        if not n:
            continue
        if r["contains_entities"]:
            rows.append({"accession": r["accession"], "bin": r["entity"],
                         "assay_family": r["assay_family"], "n_samples": n,
                         "basis": r["contains_basis"],
                         "entities": r["contains_entities"],
                         "models": r["contains_models"], "title": r["title"][:120]})
        elif r["entity"] in ("Sarcoma, unspecified / mixed cohort",):
            still[r["assay_family"]] += n
    with twrite("T31_mixed_cohort_attribution.tsv", gz=False) as f:
        w = csv.DictWriter(f, fieldnames=["accession", "bin", "assay_family", "n_samples",
                                          "basis", "entities", "models", "title"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); [w.writerow(x) for x in sorted(rows, key=lambda x: -x["n_samples"])]

    reg = sum(x["n_samples"] for x in rows if x["assay_family"] == "regulatory")
    dna = sum(x["n_samples"] for x in rows if x["assay_family"] == "DNA methylation")
    print(f"\nnow reachable from an entity page: {reg:,} regulatory, {dna:,} methylation")
    print(f"still unattributable (whole-taxonomy cohorts): {dict(still)}")
    print("  wrote T31_mixed_cohort_attribution.tsv")


if __name__ == "__main__":
    main()
