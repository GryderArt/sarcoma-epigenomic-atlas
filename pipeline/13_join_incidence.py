#!/usr/bin/env python3
"""All-ages incidence joined to the v3 atlas (v2 + GSE140686 entity recovery), so adult
and paediatric entities sit on one axis. Anchors curated from data/incidence/*.tsv with
the basis recorded; entities with no population rate anywhere are flagged, not guessed."""
import csv, collections, os
from _paths import (DATA, SAMPLES, INCIDENCE, FIGURES, WORKBOOKS, DOCS,
                    WORK, topen, twrite, dpath)
csv.field_size_limit(10**7)
O = DATA

# atlas_disease, display, US cases/yr ALL AGES, low, high, age_class, basis
A = [
 ("Osteosarcoma","Osteosarcoma",1060,900,1200,"paediatric",
  "SEER all-ages 3.1/million; ~45% of cases are 0-19"),
 ("Ewing","Ewing sarcoma",400,350,450,"paediatric",
  "~400/yr US all ages; ~60% is 0-19"),
 ("FP-RMS","FP-RMS (PAX3/7-FOXO1)",110,80,140,"paediatric",
  "~20% of all-ages RMS (~550-600/yr)"),
 ("FN-RMS","FN-RMS (fusion-negative)",400,330,460,"paediatric",
  "~75-80% of all-ages RMS"),
 ("RMS-MYOD1","MYOD1-mutant RMS",16,10,22,"paediatric",
  "2.7% of all RMS (Shern 2021, PMID 34166060)"),
 ("Rhabdoid tumor/ATRT","Rhabdoid tumour / ATRT",95,70,120,"paediatric",
  "ATRT 0.7/million 0-19 + extrarenal + renal; essentially a disease of infancy"),
 ("CCSK","Clear cell sarcoma of the kidney",20,15,25,"paediatric",
  "~3% of paediatric renal tumours; BCOR-ITD driven"),
 ("DSRCT","DSRCT",67,50,90,"paediatric",
  "NETSARC 0.197/million all ages (PMID 33630918)"),
 ("Infantile fibrosarcoma","Infantile fibrosarcoma",13,10,18,"paediatric",
  "NETSARC 0.038/million; by definition an infant tumour"),
 ("CIC-DUX4","CIC-DUX4 sarcoma",50,20,90,"paediatric",
  "no population rate; CTOS ultra-rare, median age ~40. BOUNDED ESTIMATE"),
 ("BCOR-sarcoma","BCOR-CCNB3 sarcoma",15,5,30,"paediatric",
  "no population rate; CTOS ultra-rare. BOUNDED ESTIMATE"),
 ("GCTB/Chondroblastoma","GCTB / chondroblastoma",570,400,700,"paediatric",
  "GCTB 1.7/million all ages; median age 35"),
 ("Synovial sarcoma","Synovial sarcoma",530,450,620,"both",
  "SEER 530 observed/yr; NETSARC 1.5-1.8/million"),
 ("MPNST","MPNST",390,340,450,"both",
  "NETSARC 1.00/million; Japan NCR 1.3/million"),
 ("ASPS","ASPS",42,35,50,"both","NETSARC 0.117/million"),
 ("Clear cell sarcoma","Clear cell sarcoma of soft tissue",65,41,92,"both",
  "NETSARC 0.269/million"),
 ("Epithelioid sarcoma","Epithelioid sarcoma",156,120,190,"both","NETSARC 0.455/million"),
 ("Chordoma","Chordoma",337,280,400,"both","0.97/million US"),
 ("Chondrosarcoma","Chondrosarcoma",1280,1160,1400,"both","3.4-4.1/million all ages"),
 ("GIST","GIST",4670,4000,5100,"both",
  "13.65/million adults 20+ for GI sites; registry-detected only"),
 ("EMC","Extraskeletal myxoid chondrosarcoma",75,60,90,"both","NETSARC 0.220/million"),
 ("IMT","Inflammatory myofibroblastic tumour",175,150,200,"both","NCI 150-200/yr US"),
 ("DFSP","Dermatofibrosarcoma protuberans",1400,1350,1440,"both",
  "SEER 4.1-4.2/million; two independent series agree"),
 ("Desmoid","Desmoid / aggressive fibromatosis",1740,1090,1830,"both",
  "NETSARC 5.07, Netherlands 5.36, Denmark 3.2/million"),
 ("EHE","Epithelioid haemangioendothelioma",110,79,140,"both",
  "SEER 0.230 (floor) to NETSARC 0.41/million"),
 ("Angiosarcoma","Angiosarcoma",1312,1200,1400,"adult",
  "1,312 observed cases US 2019 vs 657 in 2001 (Wagner, JAMA Netw Open 2024, PMID 38607625). "
  "The CASE COUNT doubled; the age-adjusted RATE rose ~1.6%/yr to 3.3/million, i.e. ~33% - "
  "the remainder is population growth and ageing"),
 ("Kaposi sarcoma","Kaposi sarcoma",1268,1100,1400,"adult",
  "US ~1,268/yr; 73% of the GLOBAL KS burden is in sub-Saharan Africa"),
 ("Liposarcoma-WD","Liposarcoma, well-differentiated",1100,730,1640,"adult",
  "US SEER ~730 (ALT under-reported) vs NETSARC central review ~1,640"),
 ("Liposarcoma-dediff","Liposarcoma, dedifferentiated",900,470,1740,"adult",
  "US 20% of liposarcoma (~470) vs NETSARC 41.4% (~1,740) - central review doubles it"),
 ("Liposarcoma-myxoid","Liposarcoma, myxoid / round cell",480,445,530,"adult",
  "FUS-DDIT3; US ~445, NETSARC ~530"),
 ("Liposarcoma-pleomorphic","Liposarcoma, pleomorphic",185,180,190,"adult",
  "US ~190, NETSARC ~180; below the ultra-rare threshold"),
 ("Liposarcoma-NOS","Liposarcoma, NOS / mixed",540,115,540,"adult",
  "US residual bucket 21-23%; central review collapses it ~8-fold"),
 ("UPS/MFH","Undifferentiated pleomorphic sarcoma",1500,930,2020,"adult",
  "US registry-coded ~930 vs NETSARC central review ~2,020 - coding drift from MFH"),
 ("Myxofibrosarcoma","Myxofibrosarcoma",820,700,900,"adult","NETSARC 2.386/million"),
 ("Solitary fibrous tumour","Solitary fibrous tumour",1200,1000,1410,"adult",
  "NETSARC 3.504/million extracranial + 0.62/million meningeal"),
 ("Leiomyosarcoma","Leiomyosarcoma",3310,2800,3800,"adult",
  "NETSARC 9.679/million all sites incl. uterine"),
 ("Fibrosarcoma NOS","Fibrosarcoma (adult)",36,25,50,"adult",
  "NETSARC 0.106/million - 28 cases in 4 years nationwide in France"),
 ("PEComa","PEComa (malignant)",25,20,112,"adult",
  "NETSARC malignant 0.072/million; 0.326 incl. all PEComa"),
 ("Endometrial stromal sarcoma","Endometrial stromal sarcoma",328,275,380,"adult",
  "low-grade ~275 + high-grade ~53"),
 ("LGFMS/SEF","LGFMS / SEF",229,180,280,"adult","LGFMS 0.515 + SEF 0.155/million"),
 ("Intimal sarcoma","Intimal sarcoma",60,45,75,"adult","NETSARC 0.174/million"),
 # entities now visible in the atlas for which no population rate exists anywhere
 ("Myoepithelial carcinoma","Myoepithelial carcinoma of soft tissue",0,0,0,"adult",
  "NO population rate published; NETSARC coded 2 cases in 4 years"),
 ("Ossifying fibromyxoid tumour","Ossifying fibromyxoid tumour",0,0,0,"adult",
  "NO population rate published"),
 ("Angiomatoid fibrous histiocytoma","Angiomatoid fibrous histiocytoma",0,0,0,"both",
  "NO population rate published"),
 ("AFX/PDS","Atypical fibroxanthoma / pleomorphic dermal sarcoma",0,0,0,"adult",
  "NO sarcoma-registry rate; captured as cutaneous, not sarcoma"),
]

# Adversarial verification (independent GEO/SRA/ENCODE/ChIP-Atlas/ArrayExpress/GDC sweep +
# PubMed full text) of every entity called zero-regulatory. Caveats that survived:
CAVEAT = {
 "Angiosarcoma":
   "Verified zero for HUMAN. Canine haemangiosarcoma, the accepted spontaneous model, does "
   "have regulatory data: GSE304509 (CUT&Tag, H3K27ac/H3K4me3/Kla/RNAPII-Ser5, 16 samples, "
   "PMID 41767677) and GSE150705 (ChRO-seq, 21 samples). Scope the claim to human.",
 "DFSP":
   "Verified zero for regulatory assays. One DFSP tumour has NANOPORE methylation "
   "(GSE320108 / GSM9534876), so say 'DNA methylation profiling, predominantly array-based'.",
 "Myxofibrosarcoma":
   "WEAKEST OF THE ZERO CALLS. PMID 39789291 (O'Donnell, npj Precis Oncol 2025) reports "
   "ATAC-seq on myxofibrosarcoma PDX 918122-036-R (NCI PDMR); no accession is stated and no "
   "matching series exists in GEO or SRA. Phrase as 'no retrievable data' and cite the paper. "
   "FANTOM5 CAGE also exists for MFS lines (CNhs11821 incl. NMFH-1).",
 "Liposarcoma-pleomorphic":
   "Verified zero for bona fide PLS lines (LiSa-2, LS2, NCC-PLPS1-C1, NCC-PLPS2-C1). "
   "Pre-empt SW872, which has 32 ATAC-seq + H3K27ac ChIP-seq runs and is called pleomorphic "
   "in some papers; Cellosaurus/NCIt annotate it as liposarcoma NOS and this atlas bins it "
   "as Liposarcoma-NOS.",
 "IMT":
   "Strongest zero call. Six SRA runs exist in total (RNA-seq / WXS) and Cellosaurus lists "
   "no IMT cell line at all.",
 "EMC":
   "Verified zero (H-EMC-SS, USZ20-EMC1, USZ22-EMC2). FANTOM5 CAGE exists for H-EMC-SS "
   "(CNhs10728); note also that Cellosaurus flags H-EMC-SS as lacking an EWSR1 fusion "
   "(PMID 34413129).",
 "BCOR-sarcoma":
   "Strong zero call, and citable: PMID 40841360 (Nat Commun 2025) built BCOR::CCNB3 "
   "tumoroids, profiled only WGS/WXS/RNA-seq, and states 'For BCOR-rearranged fusion "
   "proteins, no fusion-specific dataset was available'.",
 "CCSK":
   "Strong zero call. TARGET-CCSK holds RNA-seq, genotyping array, methylation array and WGS "
   "only. ChIP-Atlas's 'clear cell sarcoma' rows are all clear cell sarcoma of SOFT TISSUE "
   "(SU-CCS-1, EWSR1-ATF1), not renal. cfDNA RRBS exists at E-MTAB-8770.",
 "Intimal sarcoma":
   "Strong zero call. Thirty SRA runs in total (WXS / RNA-seq); the only cell line, PIS-1, "
   "has nothing.",
}

EPI = {"ChIP-seq","CUT&RUN","CUT&Tag","ChIP-exo","ChIP-chip","ATAC-seq","scATAC-seq",
       "DNase-seq","FAIRE-seq","MNase-seq","Hi-C","HiChIP","Micro-C","Capture-HiC",
       "ChIA-PET","4C-seq","Repli-seq","WGBS","RRBS","Methyl-array","MeDIP/hMeDIP",
       "Bisulfite-PCR"}
DNAME = {"WGBS","RRBS","Methyl-array","MeDIP/hMeDIP","Bisulfite-PCR"}
REG   = {"ChIP-seq","CUT&RUN","CUT&Tag","ChIP-exo","ChIP-chip","ATAC-seq","scATAC-seq",
         "DNase-seq","FAIRE-seq","MNase-seq","Hi-C","HiChIP","Micro-C","Capture-HiC",
         "ChIA-PET","4C-seq"}
PAT = {"primary_tumor","metastasis","recurrence"}
# patient-DERIVED material: tumour tissue plus PDX and patient-derived organoids. A synovial
# sarcoma organoid is a patient's tumour grown out, not a decades-old cell line, so the two
# categories must be reported separately rather than collapsed.
PDER = PAT | {"PDX", "organoid"}

def main():
    rows = list(csv.DictReader(topen("T4_samples_atomic.tsv"), delimiter="\t"))
    epi = [r for r in rows if r["assay_class"] in EPI and r["is_duplicate"] != "Y"]
    mc = {r["disease"]: r for r in csv.DictReader(topen("T3_entity_counts.tsv"),
                                                  delimiter="\t")}
    by = collections.defaultdict(list)
    for r in epi: by[r["disease"]].append(r)

    out = []
    for key, disp, us, lo, hi, ac, basis in A:
        m = mc.get(key, {})
        g = by.get(key, [])
        n = len(g)
        studies = collections.Counter(r["gse"] for r in g)
        top = studies.most_common(1)[0] if studies else ("", 0)
        # ChIP/CUT&RUN input and IgG controls are not profiles of anything -- exclude them
        isreg = lambda r: (r["assay_class"] in REG
                           and r["epi_target_norm"] not in ("input/none","none"))
        nreg = sum(1 for r in g if isreg(r))
        nctl = sum(1 for r in g if r["assay_class"] in REG
                   and r["epi_target_norm"] in ("input/none","none"))
        ndna = sum(1 for r in g if r["assay_class"] in DNAME)
        k27 = sum(1 for r in g if r["epi_target_norm"] == "H3K27ac")
        acc = sum(1 for r in g if r["assay_class"] in ("ATAC-seq","scATAC-seq","DNase-seq",
                                                       "FAIRE-seq"))
        d3 = sum(1 for r in g if r["assay_class"] in ("Hi-C","HiChIP","Micro-C","Capture-HiC",
                                                      "ChIA-PET","4C-seq"))
        pt = sum(1 for r in g if r["sample_type"] in PAT)
        ptreg = sum(1 for r in g if r["sample_type"] in PAT and isreg(r))
        pdreg = sum(1 for r in g if r["sample_type"] in PDER and isreg(r))
        out.append({
            "atlas_disease": key, "display_name": disp, "age_class": ac,
            "US_cases_per_year_all_ages": us, "US_low": lo, "US_high": hi,
            "no_population_rate": "Y" if us == 0 else "",
            "incidence_basis": basis,
            "verification_caveat": CAVEAT.get(key, ""),
            "epigenomic_samples": n,
            "regulatory_epigenomic_samples": nreg,
            "ChIP_input_control_samples": nctl,
            "DNA_methylation_samples": ndna,
            "pct_DNA_methylation": round(100*ndna/n, 1) if n else "",
            "H3K27ac_samples": k27, "accessibility_samples": acc, "samples_3D": d3,
            "primary_tumour_epigenomic": pt,
            "primary_tumour_regulatory": ptreg,
            "patient_derived_regulatory": pdreg,
            "n_studies": len(studies),
            "largest_study": top[0], "pct_from_largest_study": round(100*top[1]/n, 1) if n else "",
            "distinct_models": int(m.get("model_units") or 0),
            "distinct_pdx": int(m.get("pdx_units") or 0),
            "models_with_H3K27ac": int(m.get("model_H3K27ac") or 0),
            "patients_with_H3K27ac": int(m.get("patient_H3K27ac") or 0),
            "epi_samples_per_case": round(n/us, 3) if us else "",
            "regulatory_per_case": round(nreg/us, 4) if us else "",
            "cases_per_epi_sample": round(us/n, 1) if (n and us) else "",
        })
    out.sort(key=lambda r: -r["US_cases_per_year_all_ages"])
    with twrite("T13_incidence_vs_data.tsv") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()), delimiter="\t")
        w.writeheader(); [w.writerow(r) for r in out]

    rated = [r for r in out if r["US_cases_per_year_all_ages"] > 0]
    def agg(g):
        c = sum(r["US_cases_per_year_all_ages"] for r in g)
        e = sum(r["epigenomic_samples"] for r in g)
        q = sum(r["regulatory_epigenomic_samples"] for r in g)
        return c, e, q, (e/c if c else 0), (q/c if c else 0)
    ped  = [r for r in rated if r["age_class"] == "paediatric"]
    both = [r for r in rated if r["age_class"] == "both"]
    adl  = [r for r in rated if r["age_class"] == "adult"]
    print(f"{'group':12s} {'n':>4s} {'US cases/yr':>12s} {'epi':>8s} {'regulatory':>11s} "
          f"{'epi/case':>9s} {'reg/case':>9s}")
    for lab, g in (("paediatric", ped), ("both", both), ("adult", adl), ("ALL", rated)):
        c, e, q, r1, r2 = agg(g)
        print(f"  {lab:12s} {len(g):3d} {c:12,d} {e:8,d} {q:11,d} {r1:9.2f} {r2:9.3f}")
    pc, pe, pq, pr1, pr2 = agg(ped); ac_, ae, aq, ar1, ar2 = agg(adl)
    print(f"\nadult US burden is {ac_/pc:.1f}x the paediatric burden of these entities")
    print(f"paediatric entities carry {pr1/ar1:.1f}x more epigenomic samples per case")
    print(f"paediatric entities carry {pr2/ar2:.1f}x more REGULATORY epigenomic samples per case")

    print("\n--- entities with ZERO regulatory (non-methylation) epigenomics ---")
    z = [r for r in out if r["regulatory_epigenomic_samples"] == 0]
    for r in sorted(z, key=lambda x: -x["US_cases_per_year_all_ages"]):
        print(f"  {r['display_name'][:44]:46s} {r['US_cases_per_year_all_ages']:6,d}/yr  "
              f"{r['epigenomic_samples']:4d} epi ({r['pct_DNA_methylation'] or 0}% methylation)  "
              f"studies={r['n_studies']}")
    print("\n--- entities with ZERO primary-tumour regulatory epigenomics ---")
    z2 = [r for r in out if r["primary_tumour_regulatory"] == 0]
    print(f"  {len(z2)} of {len(out)} entities")
    print("\n--- least regulatory epigenomics per case (rated entities only) ---")
    for r in sorted(rated, key=lambda x: x["regulatory_per_case"])[:16]:
        print(f"  {r['display_name'][:40]:42s} {r['US_cases_per_year_all_ages']:6,d}/yr  "
              f"reg={r['regulatory_epigenomic_samples']:5d}  {r['regulatory_per_case']:.4f}/case  "
              f"K27ac={r['H3K27ac_samples']:4d}")

if __name__ == "__main__":
    main()
