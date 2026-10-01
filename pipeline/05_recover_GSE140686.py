#!/usr/bin/env python3
"""Repair GSE140686 (Koelsche et al., Nat Commun 2021, PMID 33479225).

The single largest sarcoma epigenomic resource in existence -- 1,505 DNA-methylation
arrays -- is deposited in GEO with the entity label stripped out. Every sample reads
"tissue: sarcoma". The per-case diagnosis lives only in the paper's Supplementary Data,
keyed by a REFERENCE_SAMPLE/VALIDATION_SAMPLE id that GEO buries in !Sample_description
(a field most harvesters, including our first pass, never read).

This restores the join and writes a patched atomic table.
"""
import csv, re, collections, os
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
csv.field_size_limit(10**7)
O = DATA

# ---------------------------------------------------------------- entity mapping
# (regex on institutional diagnosis) -> (atlas disease, kind)
# kind: sarcoma | benign  (benign = benign neoplasm, mimic, or non-neoplastic control)
#       nonsarcoma = malignant but not a sarcoma (classifier's negative reference set)
MAP = [
 (r"^control \(", "Non-neoplastic control", "benign"),
 (r"alveolar soft part", "ASPS", "sarcoma"),
 (r"^angiosarcoma", "Angiosarcoma", "sarcoma"),
 (r"^angiomatoid fibrous", "Angiomatoid fibrous histiocytoma", "sarcoma"),
 (r"^atypical fibroxanthoma|pleomorphic dermal sarcoma", "AFX/PDS", "sarcoma"),
 (r"atypical lipomatous|well[- ]differentiated liposarcoma|"
  r"liposarcoma \(well[- ]differentiated\)", "Liposarcoma-WD", "sarcoma"),
 (r"^chondroblastoma|giant cell tumour of bone", "GCTB/Chondroblastoma", "sarcoma"),
 (r"^chondromyxoid fibroma", "Chondromyxoid fibroma", "benign"),
 (r"^chondrosarcoma|mesenchymal chondrosarcoma", "Chondrosarcoma", "sarcoma"),
 (r"^chordoma", "Chordoma", "sarcoma"),
 (r"clear cell sarcoma \(kidney\)|clear cell sarcoma of the kidney", "CCSK", "sarcoma"),
 (r"clear cell sarcoma", "Clear cell sarcoma", "sarcoma"),
 (r"dedifferentiated liposarcoma|liposarcoma \(dedifferentiated\)", "Liposarcoma-dediff", "sarcoma"),
 (r"dermatofibrosarcoma", "DFSP", "sarcoma"),
 (r"desmoid-type", "Desmoid", "sarcoma"),
 (r"desmoplastic small round cell", "DSRCT", "sarcoma"),
 (r"endometrial stromal sarcoma", "Endometrial stromal sarcoma", "sarcoma"),
 (r"epithelioid haemangioendothelioma|^haemangioendothelioma", "EHE", "sarcoma"),
 (r"epithelioid sarcoma", "Epithelioid sarcoma", "sarcoma"),
 (r"ewing", "Ewing", "sarcoma"),
 (r"primitive neuroectodermal", "Ewing", "sarcoma"),
 (r"extraskeletal myxoid chondrosarcoma", "EMC", "sarcoma"),
 (r"fibrocartilaginous mesenchymoma", "Other rare bone tumor", "sarcoma"),
 (r"fibrous dysplasia", "Fibrous dysplasia", "benign"),
 (r"gastroin\w*stinal stromal", "GIST", "sarcoma"),
 (r"^glioblastoma", "Glioblastoma", "nonsarcoma"),
 (r"infantil\w* myofibromatosis", "Infantile myofibromatosis", "benign"),
 (r"infantile fibrosarcoma", "Infantile fibrosarcoma", "sarcoma"),
 (r"inflammatory myofibroblastic", "IMT", "sarcoma"),
 (r"^intimal sarcoma", "Intimal sarcoma", "sarcoma"),
 (r"kaposi", "Kaposi sarcoma", "sarcoma"),
 (r"langerhans cell histiocytosis", "Langerhans cell histiocytosis", "nonsarcoma"),
 (r"^leiomyoma", "Leiomyoma", "benign"),
 (r"^leiomyosarcoma", "Leiomyosarcoma", "sarcoma"),
 (r"lipoblastomatosis|^lipoma", "Lipoma/lipoblastomatosis", "benign"),
 (r"myxoid liposarcoma|liposarcoma \(myxoid\)", "Liposarcoma-myxoid", "sarcoma"),
 (r"pleomorphic liposarcoma|liposarcoma \(pleomorphic\)", "Liposarcoma-pleomorphic", "sarcoma"),
 (r"^liposarcoma$", "Liposarcoma-NOS", "sarcoma"),
 (r"low-grade fibromyxoid|sclerosing epithelioid", "LGFMS/SEF", "sarcoma"),
 (r"malignant mixed mesodermal", "Carcinosarcoma", "nonsarcoma"),
 (r"malignant peripheral nerve sheath", "MPNST", "sarcoma"),
 (r"rhabdoid tumour", "Rhabdoid tumor/ATRT", "sarcoma"),
 (r"^melanoma", "Melanoma", "nonsarcoma"),
 (r"^myoepithelioma", "Myoepithelial tumor", "sarcoma"),
 (r"^myopericytoma|^angioleiomyoma", "Angioleiomyoma/myopericytoma", "benign"),
 (r"^myositis", "Myositis ossificans/proliferans", "benign"),
 (r"^myxofibrosarcoma", "Myxofibrosarcoma", "sarcoma"),
 (r"^neurofibroma", "Neurofibroma", "benign"),
 (r"nodular fasciitis", "Nodular fasciitis", "benign"),
 (r"ossifying fibromyxoid", "Ossifying fibromyxoid tumor", "sarcoma"),
 (r"^osteoblastoma", "Osteoblastoma", "benign"),
 (r"^osteosarcoma", "Osteosarcoma", "sarcoma"),
 (r"^pecoma", "PEComa", "sarcoma"),
 (r"rhabdomyosarcoma", "__RMS__", "sarcoma"),
 (r"small (blue )?round cell tumour with bcor", "BCOR-sarcoma", "sarcoma"),
 (r"small blue round cell tumour with cic", "CIC-DUX4", "sarcoma"),
 (r"solitary fibrous", "Solitary fibrous tumor", "sarcoma"),
 (r"^schwannoma", "Schwannoma", "benign"),
 (r"spindle cell carcinoma", "Spindle cell carcinoma", "nonsarcoma"),
 (r"spindle cell hemangioma", "Spindle cell hemangioma", "benign"),
 (r"squamous cell carcinoma", "Squamous cell carcinoma", "nonsarcoma"),
 (r"^synovial sarcoma", "Synovial sarcoma", "sarcoma"),
 (r"undifferentiated pleomorphic", "UPS/MFH", "sarcoma"),
 (r"undifferentiated epithelioid sarcoma", "Epithelioid sarcoma", "sarcoma"),
 (r"undifferentiated sarcoma|^sarcoma,? ?(not otherwise specified|nos)|"
  r"malignant tumour not otherwise", "Sarcoma NOS", "sarcoma"),
]

# RMS sub-binning follows the user's rule: MYOD1 mutation trumps, then fusion status.
# Histology alone does not establish fusion status, so the evidence level is recorded.
def rms_bin(dx, mc):
    m = mc.lower()
    if "myod1" in m:                        return "RMS-MYOD1", "methylation_class_MYOD1"
    if "alv" in m:                          return "FP-RMS",    "methylation_class_alveolar"
    if "emb" in m:                          return "FN-RMS",    "methylation_class_embryonal"
    d = dx.lower()
    if "alveolar" in d:                     return "FP-RMS",    "histology_alveolar"
    if "embryonal" in d:                    return "FN-RMS",    "histology_embryonal"
    if "spindle" in d or "scleros" in d:    return "RMS-NOS",   "histology_spindle_sclerosing"
    return "RMS-NOS", "histology_NOS"

COMP = [(re.compile(p, re.I), d, k) for p, d, k in MAP]
def classify(dx, mc):
    for rx, d, k in COMP:
        if rx.search(dx):
            if d == "__RMS__":
                b, ev = rms_bin(dx, mc); return b, "sarcoma", ev
            return d, k, "institutional_diagnosis"
    # fall back on the methylation class when the institutional diagnosis is unhelpful
    return "Sarcoma NOS", "sarcoma", "unmapped"

MANIF = {"Primary": "primary_tumor", "Metastasis": "metastasis",
         "Recurrence": "recurrence", "Recurrence/Metastasis": "metastasis", "": "primary_tumor"}

def main():
    lab = {r["gsm"]: r for r in csv.DictReader(topen("GSE140686_sample_keys.tsv"), delimiter="\t")}
    out = []
    for g, r in lab.items():
        d, kind, ev = classify(r["institutional_diagnosis"], r["methylation_class"])
        r["atlas_disease"] = d; r["entity_kind"] = kind; r["subtype_evidence"] = ev
        r["sample_type"] = MANIF.get(r["manifestation"], "primary_tumor")
        out.append(r)
    cols = ["gsm","sample_key","set","institutional_diagnosis","methylation_class",
            "atlas_disease","entity_kind","subtype_evidence","sample_type","site","manifestation"]
    with twrite("GSE140686_recovered_diagnoses.tsv") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); [w.writerow(r) for r in out]

    unm = [r for r in out if r["subtype_evidence"] == "unmapped"]
    print(f"{len(out)} samples relabeled; {len(unm)} unmapped")
    for r in unm[:10]: print("   UNMAPPED:", r["institutional_diagnosis"], "|", r["methylation_class"])
    k = collections.Counter(r["entity_kind"] for r in out)
    print("kind:", dict(k))
    print("\nsarcoma entities recovered (n>=1):")
    c = collections.Counter(r["atlas_disease"] for r in out if r["entity_kind"] == "sarcoma")
    for d, n in c.most_common(): print(f"  {n:5d}  {d}")

    # ---- patch the atomic table
    rows = list(csv.DictReader(topen("T4_samples_atomic_v2.tsv"), delimiter="\t"))
    hdr = list(rows[0].keys())
    for c2 in ("subtype_evidence","entity_kind"):
        if c2 not in hdr: hdr.append(c2)
    n = 0
    for r in rows:
        r.setdefault("subtype_evidence", ""); r.setdefault("entity_kind", "sarcoma")
        p = lab.get(r["gsm"])
        if not p: continue
        r["disease"] = p["atlas_disease"]; r["sample_type"] = p["sample_type"]
        r["subtype_evidence"] = p["subtype_evidence"]; r["entity_kind"] = p["entity_kind"]
        r["disease_source"] = "GSE140686 supplement (Koelsche 2021) recovered via Sample_description"
        n += 1
    with twrite("T4_samples_atomic.tsv", gz=True) as f:
        w = csv.DictWriter(f, fieldnames=hdr, delimiter="\t", extrasaction="ignore")
        w.writeheader(); [w.writerow(r) for r in rows]
    print(f"\npatched {n} rows -> T4_samples_atomic.tsv ({len(rows)} total)")

if __name__ == "__main__":
    main()
