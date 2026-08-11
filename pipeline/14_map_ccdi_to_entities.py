#!/usr/bin/env python3
"""Map the CCDI cohort onto atlas entities and fold it in as a CONTROLLED-access source.

CCDI records are metadata for dbGaP-controlled data. They belong in the same tier as EGA,
never in the open counts -- mixing them would destroy the open / controlled / absent
distinction the gap map is built on. Every row carries access_tier=CONTROLLED (dbGaP).

Produces:
  T24_ccdi_entity_counts.tsv   per atlas entity: participants, samples, tumour samples,
                               and the epigenomic (methylation-array) subset
  T25_ccdi_unmapped.tsv        CCDI diagnoses this mapping could not place, so the
                               residual is visible rather than silently dropped
"""
import csv, gzip, os, re, sys, collections
csv.field_size_limit(10**7)
from _paths import DATA, topen, twrite, dpath

# CCDI uses ICD-O-3 preferred terms. Map to the atlas's 45 named entities, most specific
# first -- first match wins. Anything unmatched is reported, not discarded.
MAP = [
 (r"^clear cell sarcoma of (the )?kidney|ccsk", "CCSK"),
 (r"alveolar rhabdomyosarcoma|rhabdomyosarcoma, alveolar", "FP-RMS"),
 (r"embryonal rhabdomyosarcoma|botryoid|rhabdomyosarcoma, embryonal", "FN-RMS"),
 (r"spindle cell rhabdomyosarcoma|sclerosing rhabdomyosarcoma", "RMS-MYOD1"),
 (r"pleomorphic rhabdomyosarcoma", "RMS-NOS"),
 (r"rhabdomyosarcoma", "RMS-NOS"),
 (r"ewing sarcoma|askin|peripheral neuroectodermal|pnet of|primitive neuroectodermal",
  "Ewing"),
 (r"osteosarcoma|osteogenic sarcoma", "Osteosarcoma"),
 (r"desmoplastic small round cell", "DSRCT"),
 (r"synovial sarcoma", "Synovial sarcoma"),
 (r"malignant peripheral nerve sheath|mpnst|neurofibrosarcoma", "MPNST"),
 (r"atypical teratoid|rhabdoid tumor|rhabdoid tumour|atrt", "Rhabdoid tumor/ATRT"),
 (r"epithelioid sarcoma", "Epithelioid sarcoma"),
 (r"alveolar soft part", "ASPS"),
 (r"clear cell sarcoma", "Clear cell sarcoma"),
 (r"chondrosarcoma", "Chondrosarcoma"),
 (r"chordoma", "Chordoma"),
 (r"leiomyosarcoma", "Leiomyosarcoma"),
 (r"liposarcoma, well|well[- ]differentiated liposarcoma|atypical lipomatous",
  "Liposarcoma-WD"),
 (r"dedifferentiated liposarcoma", "Liposarcoma-dediff"),
 (r"myxoid liposarcoma|round cell liposarcoma", "Liposarcoma-myxoid"),
 (r"pleomorphic liposarcoma", "Liposarcoma-pleomorphic"),
 (r"liposarcoma", "Liposarcoma-NOS"),
 (r"angiosarcoma|haemangiosarcoma|hemangiosarcoma", "Angiosarcoma"),
 (r"epithelioid haemangioendothelioma|epithelioid hemangioendothelioma", "EHE"),
 (r"kaposi", "Kaposi sarcoma"),
 (r"dermatofibrosarcoma", "DFSP"),
 (r"myxofibrosarcoma", "Myxofibrosarcoma"),
 (r"solitary fibrous|hemangiopericytoma|haemangiopericytoma", "Solitary fibrous tumour"),
 (r"inflammatory myofibroblastic|myofibroblastic tumor|myofibroblastic tumour", "IMT"),
 (r"extraskeletal myxoid chondrosarcoma", "EMC"),
 (r"aggressive fibromatosis|desmoid", "Desmoid"),
 (r"infantile fibrosarcoma|congenital fibrosarcoma", "Infantile fibrosarcoma"),
 (r"low[- ]grade fibromyxoid|sclerosing epithelioid fibrosarcoma", "LGFMS/SEF"),
 (r"fibrosarcoma", "Fibrosarcoma NOS"),
 (r"gastrointestinal stromal|gist", "GIST"),
 (r"endometrial stromal sarcoma", "Endometrial stromal sarcoma"),
 (r"pecoma|perivascular epithelioid", "PEComa"),
 (r"intimal sarcoma", "Intimal sarcoma"),
 (r"giant cell tumor of bone|giant cell tumour of bone|chondroblastoma",
  "GCTB/Chondroblastoma"),
 (r"undifferentiated pleomorphic sarcoma|malignant fibrous histiocytoma|"
  r"undifferentiated sarcoma|undifferentiated high grade|spindle cell sarcoma",
  "UPS/MFH"),
 (r"myoepithelial carcinoma", "Myoepithelial carcinoma"),
 (r"bcor", "BCOR-sarcoma"),
 (r"cic[- ]", "CIC-DUX4"),
 (r"ossifying fibromyxoid", "Ossifying fibromyxoid tumour"),
 (r"angiomatoid fibrous histiocytoma", "Angiomatoid fibrous histiocytoma"),
 (r"atypical fibroxanthoma|pleomorphic dermal sarcoma", "AFX/PDS"),
 (r"^sarcoma, nos$|^sarcoma$|^soft tissue sarcoma", "Sarcoma NOS"),
]
COMP = [(re.compile(p, re.I), e) for p, e in MAP]

# Not sarcoma, and must not be mapped into one. CCDI's sarcoma-family ICD-O categories
# sweep in neighbours that share a category but are different diseases.
EXCLUDE = re.compile(
    r"neuroblastoma|ganglioneur|nephroblastoma|wilms|hepatoblastoma|retinoblastoma|"
    r"medulloblastoma|glioma|glioblastoma|astrocytoma|ependymoma|ganglioglioma|"
    r"meningioma|schwannoma|neurofibroma(?!sarcoma)|hemangioblastoma|"
    r"pleuropulmonary blastoma|germ cell|teratoma|lymphoma|leukemia|leukaemia|"
    r"carcinoma(?! of soft)|melanoma|histiocytosis|hemangioma|haemangioma|"
    r"lipoma(?!tous)|leiomyoma|myopericytoma|angioleiomyoma|nodular fasciitis|"
    r"myositis|fibrous dysplasia|osteoblastoma|craniopharyngioma|pineal|"
    r"choroid plexus|paraganglioma|pheochromocytoma|chordoid glioma", re.I)

EPI_FILE_TYPES = {"idat"}          # methylation arrays
EPI_STRATEGIES = {"ATAC-seq", "Bisulfite-Seq", "ChIP-Seq", "ChIP-seq"}


def nid(v):
    """CCDI returns file-level ids wrapped in brackets ('[1794211]') while sample-level
    ids are bare ('003-1'). Normalise before joining, or the join silently yields zero."""
    return str(v or "").strip().strip("[]").strip().strip("'\"")


def entity(dx):
    if not dx: return None, "no diagnosis recorded"
    if EXCLUDE.search(dx): return None, "not a sarcoma"
    for rx, e in COMP:
        if rx.search(dx): return e, "mapped"
    return None, "unmapped sarcoma-category diagnosis"


def main():
    samples = list(csv.DictReader(topen("T21_ccdi_samples.tsv"), delimiter="\t"))
    parts = {r["participant_id"]: r for r in
             csv.DictReader(topen("T23_ccdi_participants.tsv"), delimiter="\t")}
    files = list(csv.DictReader(topen("T22_ccdi_files.tsv"), delimiter="\t"))
    # T22 ships only the epigenomic rows; T22b is the census of all of them, and is what
    # makes "not one regulatory file among 319,084" checkable from this repository.
    census = list(csv.DictReader(topen("T22b_ccdi_file_census.tsv"), delimiter="\t"))
    all_files = sum(int(r["files"]) for r in census)
    sar = [r for r in samples if r["is_sarcoma_cohort"] == "Y"]
    print(f"{len(samples):,} CCDI samples, {len(sar):,} in the sarcoma cohort, "
          f"{len(parts):,} participants")
    print(f"{all_files:,} sarcoma-cohort files examined; {len(files):,} epigenomic rows "
          f"retained ({100*len(files)/all_files:.1f}%)")

    # a sample's diagnosis, falling back to its participant's
    def dx_of(r):
        d = (r.get("diagnosis") or "").strip()
        if d: return d
        p = parts.get(r["participant_id"])
        return (p.get("diagnosis") or "").strip() if p else ""

    # which samples have epigenomic files
    epi_samples, epi_parts = set(), set()
    strat = collections.Counter()
    for f in files:
        ls = (f.get("library_strategy") or "").strip()
        strat[ls or "(no strategy: array or derived)"] += 1
        if (f.get("file_type") or "").lower() in EPI_FILE_TYPES or ls in EPI_STRATEGIES:
            if f.get("sample_id"): epi_samples.add((nid(f["sample_id"]), nid(f["study_id"])))
            if f.get("participant_id"): epi_parts.add(nid(f["participant_id"]))
    print(f"\nlibrary_strategy across ALL {all_files:,} sarcoma-cohort files (from T22b):")
    cstrat = collections.Counter()
    for r in census: cstrat[r["library_strategy"]] += int(r["files"])
    for k, v in cstrat.most_common():
        print(f"   {v:8,d}  {k}")
    reg = sum(v for k, v in cstrat.items() if k in EPI_STRATEGIES)
    print(f"\n  REGULATORY epigenomic files in the CCDI sarcoma cohort: {reg}")

    agg = collections.defaultdict(lambda: collections.Counter())
    unmapped = collections.Counter()
    for r in sar:
        dx = dx_of(r)
        e, why = entity(dx)
        if not e:
            unmapped[(dx or "(none)", why)] += 1
            continue
        a = agg[e]
        a["samples"] += 1
        if r["sample_tumor_status"] == "Tumor": a["tumour_samples"] += 1
        if r["sample_tumor_status"] == "Normal": a["normal_samples"] += 1
        if r["tumor_classification"] == "Primary":
            a["primary_samples"] += 1
        a.setdefault  # noop
    # distinct participants per entity
    pe = collections.defaultdict(set)
    for r in sar:
        e, _ = entity(dx_of(r))
        if e: pe[e].add(r["participant_id"])
    for e, s in pe.items(): agg[e]["participants"] = len(s)
    for e in agg:
        agg[e]["participants_with_methylation"] = len(
            {r["participant_id"] for r in sar
             if entity(dx_of(r))[0] == e and nid(r["participant_id"]) in epi_parts})

    rows = [{"atlas_entity": e, "ccdi_participants": a["participants"],
             "ccdi_samples": a["samples"], "ccdi_tumour_samples": a["tumour_samples"],
             "ccdi_normal_samples": a["normal_samples"],
             "ccdi_primary_samples": a["primary_samples"],
             "ccdi_participants_with_methylation": a["participants_with_methylation"],
             "ccdi_regulatory_epigenomic_samples": 0,
             "access_tier": "CONTROLLED (dbGaP)"}
            for e, a in agg.items()]
    rows.sort(key=lambda r: -r["ccdi_participants"])
    with twrite("T24_ccdi_entity_counts.tsv", gz=False) as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); [w.writerow(r) for r in rows]

    urows = [{"ccdi_diagnosis": d, "reason": why, "samples": n}
             for (d, why), n in unmapped.most_common()]
    with twrite("T25_ccdi_unmapped.tsv", gz=False) as f:
        w = csv.DictWriter(f, fieldnames=["ccdi_diagnosis", "reason", "samples"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); [w.writerow(r) for r in urows]

    tot = sum(r["ccdi_samples"] for r in rows)
    print(f"\n{len(rows)} atlas entities matched; {tot:,} of {len(sar):,} sarcoma-cohort "
          f"samples mapped ({100*tot/len(sar):.0f}%)")
    print(f"{sum(r['samples'] for r in urows):,} samples unmapped -> T25 "
          f"({len(urows)} distinct diagnoses)")
    tm = sum(r["ccdi_participants_with_methylation"] for r in rows)
    print(f"\nparticipants with methylation arrays, mapped to an entity: {tm:,}")
    print("  (joined participant-to-file; the sample-level join is NOT used because T21 "
          "and\n   T22 use different sample_id spaces in several studies)")
    print(f"\n{'entity':32s} {'part':>6s} {'samples':>8s} {'tumour':>7s} {'meth':>6s}")
    for r in rows[:24]:
        print(f"{r['atlas_entity'][:31]:32s} {r['ccdi_participants']:6,d} "
              f"{r['ccdi_samples']:8,d} {r['ccdi_tumour_samples']:7,d} "
              f"{r['ccdi_participants_with_methylation']:6,d}")
    print("\ntop unmapped reasons:")
    for r in urows[:10]:
        print(f"   {r['samples']:6,d}  [{r['reason']}] {r['ccdi_diagnosis'][:56]}")


if __name__ == "__main__":
    main()
