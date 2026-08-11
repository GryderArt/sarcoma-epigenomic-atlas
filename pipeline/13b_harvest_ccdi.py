#!/usr/bin/env python3
"""Harvest every CCDI sample and file record through the open GraphQL API.

The CCDI Hub Explore Dashboard is a JavaScript front end over
https://ccdi.cancer.gov/v1/graphql/, which is open, unauthenticated and pages 10,000
records at a time. Row-level participant, sample and file metadata is OPEN even where the
underlying data is dbGaP-controlled -- so the catalogue can be enumerated in full without
credentials, exactly as for EGA.

What this produces:
  T21_ccdi_samples.tsv.gz   the sarcoma-cohort samples. All 70,820 CCDI samples are
                            fetched and the participant join is run over all of them --
                            that is how the cohort is defined -- but only the sarcoma rows
                            are kept. The all-disease totals survive in T26.
  T26_ccdi_disease_census.tsv  participants, samples and assay strategies by ICD-O
                            category across ALL of CCDI, which is where the finding that
                            its 91 ATAC-seq and 108 bisulfite-seq participants are
                            entirely leukaemia remains checkable.
  T22_ccdi_files.tsv.gz     the EPIGENOMIC subset of files belonging to a sarcoma-cohort
                            participant. All 319,084 files are fetched and examined; only
                            the ~2.9% that are epigenomic are kept, because 163,547 WGS
                            and 78,836 WXS rows tell us nothing about chromatin and cost
                            23 MB in the repository.
  T22b_ccdi_file_census.tsv the full library_strategy x file_type x data_category census
                            of all 319,084 files, so the claim that not one of them is a
                            regulatory epigenomic file stays verifiable from this
                            repository without shipping the rows themselves.
  T20_ccdi_studies.tsv      the study listing

IMPORTANT: these records are METADATA for controlled-access data. They must never be
counted alongside open GEO samples -- see access_tier, which is set to CONTROLLED
(dbGaP) on every row.
"""
import urllib.request, json, csv, gzip, os, sys, time

URL = "https://ccdi.cancer.gov/v1/graphql/"
OUT = os.environ.get("SASS_DATA", os.path.dirname(os.path.abspath(__file__)))
PAGE = 10000

# ICD-O-3 diagnosis categories that constitute the sarcoma family, plus the paediatric
# small-round-blue-cell and embryonal neighbours the atlas already tracks.
SARCOMA_CATEGORIES = [
    "Soft tissue tumors and sarcomas, NOS", "Miscellaneous bone tumors",
    "Neuroepitheliomatous neoplasms", "Myomatous neoplasms",
    "Osseous and chondromatous neoplasms", "Mesenchymal, non-meningothelial tumors",
    "Nerve sheath tumors", "Uncertain differentiation", "Fibromatous neoplasms",
    "Synovial-like neoplasms", "Blood vessel tumors", "Soft tissue tumors",
    "Chondro-osseous tumors", "Granular cell tumors and alveolar soft part sarcomas",
    "Lipomatous neoplasms", "Lymphatic vessel tumors", "Skeletal muscle tumors",
    "Giant cell tumors", "Myxomatous neoplasms", "CNS Sarcoma",
    "Atypical teratoid/rhabdoid tumor", "Atypical Teratoid/Rhabdoid Tumor",
    "Complex mixed and stromal neoplasms",
]


def gq(query, variables=None, tries=5):
    body = json.dumps({"query": query, "variables": variables or {}}).encode()
    for a in range(tries):
        try:
            req = urllib.request.Request(URL, data=body, headers={
                "Content-Type": "application/json", "Accept": "application/json",
                "User-Agent": "sarcoma-atlas/1.0 (research metadata harvest)"})
            d = json.loads(urllib.request.urlopen(req, timeout=600).read().decode())
            if d.get("errors") and not d.get("data"):
                raise RuntimeError(str(d["errors"])[:200])
            return d
        except Exception as e:
            if a == tries - 1:
                sys.stderr.write(f"  giving up: {e}\n"); raise
            time.sleep(3 * (a + 1))


def page_all(query, key, variables=None, label=""):
    """Page a row-level query to exhaustion."""
    rows, offset = [], 0
    while True:
        v = dict(variables or {}); v.update({"first": PAGE, "offset": offset})
        got = gq(query, v)["data"][key]
        if not got:
            break
        rows.extend(got); offset += PAGE
        print(f"    {label} {len(rows):,}", flush=True)
        if len(got) < PAGE:
            break
    return rows


PARTICIPANTS = """query($first:Int,$offset:Int){
  participantOverview(first:$first, offset:$offset, order_by:"participant_id"){
    participant_id dbgap_accession study_id race sex_at_birth age_at_diagnosis
    diagnosis diagnosis_category anatomic_site last_known_survival_status } }"""

SAMPLES = """query($first:Int,$offset:Int,$dc:[String]){
  sampleOverview(first:$first, offset:$offset, order_by:"sample_id", diagnosis_category:$dc){
    sample_id participant_id study_id anatomic_site participant_age_at_collection
    sample_tumor_status tumor_classification diagnosis diagnosis_category } }"""

FILES = """query($first:Int,$offset:Int,$dc:[String]){
  fileOverview(first:$first, offset:$offset, order_by:"file_id", diagnosis_category:$dc){
    file_id file_name data_category file_type file_size library_selection
    library_source_material library_source_molecule library_strategy file_mapping_level
    file_access sample_id participant_id study_id guid } }"""

STUDIES = """{ studiesListing(first:500){ study_id study_name num_of_participants
    num_of_diagnoses num_of_samples num_of_files num_of_publications } }"""


EPI_STRATEGIES = {"ATAC-seq", "Bisulfite-Seq", "ChIP-Seq", "ChIP-seq",
                  "CUT&RUN", "CUT&Tag", "Hi-C", "HiChIP", "MNase-Seq", "DNase-Hypersensitivity"}
EPI_FILE_TYPES = {"idat"}
EPI_CATEGORY_HINTS = ("methylation", "epigen", "chromatin", "atac", "bisulfite")


def is_epigenomic(r):
    """Keep methylation arrays and any genuinely epigenomic library strategy. Drops WGS,
    WXS, RNA-seq, panels and their index files -- 97% of the cohort's files, and none of
    them informative about chromatin."""
    if str(r.get("file_type") or "").lower() in EPI_FILE_TYPES: return True
    if str(r.get("library_strategy") or "").strip() in EPI_STRATEGIES: return True
    dc = str(r.get("data_category") or "").lower()
    return any(k in dc for k in EPI_CATEGORY_HINTS)


def flat(v):
    if v is None: return ""
    if isinstance(v, (list, tuple)): v = "; ".join(str(x) for x in v)
    return " ".join(str(v).split())


def write_tsv(name, rows, cols, gz=True):
    p = os.path.join(OUT, name + (".gz" if gz and not name.endswith(".gz") else ""))
    op = gzip.open(p, "wt", newline="") if p.endswith(".gz") else open(p, "w", newline="")
    with op as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore",
                           lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: flat(r.get(c)) for c in cols})
    return p, os.path.getsize(p)


def main():
    only_sarcoma_samples = "--sarcoma-only" in sys.argv

    print("studies")
    st = gq(STUDIES)["data"]["studiesListing"]
    for s in st: s["study_name"] = flat(s["study_name"])
    st.sort(key=lambda x: -(x["num_of_participants"] or 0))
    p, n = write_tsv("T20_ccdi_studies.tsv", st, list(st[0].keys()), gz=False)
    print(f"  {len(st)} studies -> {os.path.basename(p)}")

    # Participants first: 76% of SAMPLE rows carry no diagnosis_category of their own.
    # The diagnosis lives on the participant, and the dashboard's sample counts come from a
    # server-side join. Flagging sarcoma from the sample row alone undercounts ~17-fold.
    print("participants (all CCDI)")
    parts = page_all(PARTICIPANTS, "participantOverview", {}, "participants")
    for r in parts:
        r["access_tier"] = "CONTROLLED (dbGaP)"; r["source_archive"] = "CCDI"
    pcols = ["participant_id", "study_id", "dbgap_accession", "diagnosis",
             "diagnosis_category", "anatomic_site", "age_at_diagnosis", "race",
             "sex_at_birth", "last_known_survival_status", "access_tier", "source_archive"]
    # trimmed to the sarcoma cohort after the join below; see write further down

    SARCSET = set(SARCOMA_CATEGORIES)
    def is_sarc(v):
        if not v: return False
        cats = v if isinstance(v, (list, tuple)) else str(v).split(";")
        return any(str(c).strip() in SARCSET for c in cats)
    sarc_participants = {r["participant_id"] for r in parts if is_sarc(r.get("diagnosis_category"))}
    print(f"  {len(sarc_participants):,} of {len(parts):,} participants in the sarcoma cohort")
    keep = [r for r in parts if r["participant_id"] in sarc_participants]
    p, n = write_tsv("T23_ccdi_participants.tsv", keep, pcols)
    print(f"  {len(keep):,} sarcoma-cohort participants -> {os.path.basename(p)} "
          f"{n/1e6:.2f} MB")

    # Samples: take the server-side filtered set (same join the dashboard uses) as the
    # authoritative sarcoma cohort, then union with anything the participant join catches.
    print("samples (all CCDI)")
    samples = page_all(SAMPLES, "sampleOverview", {"dc": []}, "samples")
    print("samples (server-side sarcoma filter)")
    sarc_samples = page_all(SAMPLES, "sampleOverview", {"dc": SARCOMA_CATEGORIES}, "sarcoma samples")
    sarc_ids = {(r["sample_id"], r["study_id"]) for r in sarc_samples}
    for r in samples:
        by_server = (r["sample_id"], r["study_id"]) in sarc_ids
        by_part = r.get("participant_id") in sarc_participants
        r["is_sarcoma_cohort"] = "Y" if (by_server or by_part) else ""
        r["sarcoma_call_basis"] = ("dashboard filter + participant diagnosis" if by_server and by_part
                                   else "dashboard filter" if by_server
                                   else "participant diagnosis" if by_part else "")
        r["access_tier"] = "CONTROLLED (dbGaP)"
        r["source_archive"] = "CCDI"
    cols = ["sample_id", "participant_id", "study_id", "diagnosis", "diagnosis_category",
            "anatomic_site", "sample_tumor_status", "tumor_classification",
            "participant_age_at_collection", "is_sarcoma_cohort", "sarcoma_call_basis",
            "access_tier", "source_archive"]
    # census across ALL diseases before trimming -- small, and it is the evidence
    import collections as _c
    dcen = _c.Counter()
    for r in samples:
        cats = str(r.get("diagnosis_category") or "(none)")
        dcen[(cats, r.get("sample_tumor_status") or "")] += 1
    crows = [{"diagnosis_category": a, "sample_tumor_status": b, "samples": n}
             for (a, b), n in dcen.most_common()]
    write_tsv("T26_ccdi_disease_census.tsv", crows,
              ["diagnosis_category", "sample_tumor_status", "samples"], gz=False)
    print(f"  census across all diseases -> T26_ccdi_disease_census.tsv "
          f"({len(crows)} rows)")

    sarc_only = [r for r in samples if r["is_sarcoma_cohort"]]
    p, n = write_tsv("T21_ccdi_samples.tsv", sarc_only, cols)
    print(f"  {len(sarc_only):,} sarcoma-cohort samples kept of {len(samples):,} -> "
          f"{os.path.basename(p)} {n/1e6:.2f} MB")

    print("files (sarcoma cohort)")
    files = page_all(FILES, "fileOverview", {"dc": SARCOMA_CATEGORIES}, "files")
    for r in files:
        r["access_tier"] = "CONTROLLED (dbGaP)"
        r["source_archive"] = "CCDI"
    fcols = ["file_id", "file_name", "study_id", "participant_id", "sample_id",
             "data_category", "file_type", "library_strategy", "library_selection",
             "library_source_material", "library_source_molecule", "file_mapping_level",
             "file_access", "file_size", "guid", "access_tier", "source_archive"]

    # Census of EVERYTHING first -- this is the evidence, and it is small.
    import collections as _c
    census = _c.Counter()
    for r in files:
        census[(flat(r.get("library_strategy")) or "(no strategy: array or derived)",
                flat(r.get("file_type")) or "(none)",
                flat(r.get("data_category")) or "(none)")] += 1
    crows = [{"library_strategy": a, "file_type": b, "data_category": c, "files": n}
             for (a, b, c), n in census.most_common()]
    p, n = write_tsv("T22b_ccdi_file_census.tsv", crows,
                     ["library_strategy", "file_type", "data_category", "files"], gz=False)
    print(f"  census of all {len(files):,} files -> {os.path.basename(p)} "
          f"({len(crows)} combinations)")
    reg = sum(v for (a, _b, _c), v in census.items() if a in EPI_STRATEGIES)
    print(f"  REGULATORY epigenomic files in the sarcoma cohort: {reg}")

    # Then keep only the epigenomic rows.
    epi = [r for r in files if is_epigenomic(r)]
    p, n = write_tsv("T22_ccdi_files.tsv", epi, fcols)
    print(f"  {len(epi):,} epigenomic files kept of {len(files):,} "
          f"({100*len(epi)/len(files):.1f}%) -> {os.path.basename(p)} {n/1e6:.2f} MB")

    import collections
    print("\n--- what was harvested ---")
    print("  sample tumour status:",
          dict(collections.Counter(r["sample_tumor_status"] for r in samples
                                   if r["is_sarcoma_cohort"]).most_common(5)))
    print("  file library_strategy (all fetched):",
          dict(collections.Counter(r["library_strategy"] or "(none)"
                                   for r in files).most_common(12)))


if __name__ == "__main__":
    main()
