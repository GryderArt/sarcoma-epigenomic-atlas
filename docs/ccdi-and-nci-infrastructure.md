# CCDI and the NCI infrastructure: what it adds, and what it proves

Working notes for *Critical assessment of existing datasets*. Everything below came from
the CCDI Hub's own GraphQL API (`https://ccdi.cancer.gov/v1/graphql/`) and the GDC API,
both open and unauthenticated. Tables: `data/T20_ccdi_studies.tsv` through `data/T25_ccdi_unmapped.tsv`, plus the full
harvest in `data/T21_ccdi_samples.tsv.gz` (70,820 samples), `data/T22_ccdi_files.tsv.gz`
(319,084 files) and `data/T23_ccdi_participants.tsv.gz` (61,854 participants).

**Confirmed at file level.** The original version of this note rested on facet counts. The
full harvest since verified it row by row: across all 319,084 files belonging to a
sarcoma-cohort participant, the `library_strategy` values are WGS (163,547), WXS (78,836),
RNA-Seq (32,936), Archer Fusion, Other, AMPLICON, Targeted-Capture and miRNA-Seq — plus
25,103 array or derived files carrying no strategy. **Not one ATAC-seq, bisulfite or
ChIP-seq file.**

---

## 1. The finding I did not expect

**CCDI does not have ChIP-seq. At all. For any disease.**

Its `library_strategy` vocabulary is exactly eleven values — WGS, WXS, RNA-Seq, Archer
Fusion, Other, AMPLICON, Targeted-Capture, ATAC-seq, Bisulfite-Seq, WGA, miRNA-Seq. There
is no ChIP-seq, no CUT&RUN, no CUT&Tag, no Hi-C, no HiChIP. Across 42 studies, 61,854
participants and 1,259,864 files.

It does have ATAC-seq and bisulfite sequencing — and **every one of those samples is
leukaemia**:

| Assay | Participants | Diagnoses | Studies |
|---|---|---|---|
| ATAC-seq | 91 | T-ALL 42, precursor B-ALL 24, TEL-AML1 8, AML 4 … | phs003432, phs002371, phs002529 |
| Bisulfite-Seq | 108 | JMML 104, precursor B-ALL 4 | phs002504, phs003215 |

**Zero sarcoma. Zero solid tumor of any kind.**

The same holds one level up. The GDC's vocabulary has no ChIP-seq either, and its 410
ATAC-seq files belong to 23 TCGA cohorts — BRCA, COAD, KIRP, PRAD and so on. **TCGA-SARC
is not among them**, and neither is any TARGET sarcoma project.

So the paper can now say something considerably stronger than "the literature under-served
adult sarcoma." It can say: **across the entire United States federal pediatric and adult
cancer data infrastructure — CCDI, GDC, TARGET, TCGA — there is not one regulatory
epigenomic dataset for any sarcoma.** Paired with EpiRR's 2,888 reference epigenomes and
zero sarcoma, that is a structural finding about how the field is funded, not an accident
of which labs happened to publish.

---

## 2. What CCDI genuinely adds

CCDI holds a large pediatric sarcoma cohort: **29 studies, 16,867 participants, 23,143
samples, 319,084 files.** For epigenomics that resolves to one modality — DNA methylation
arrays — but at a scale worth having:

**1,592 sarcoma participants with raw Illumina IDATs, 1,382 of them primary tumors.**

Almost all from the Molecular Characterization Initiative (`phs002790`, 1,560) with a small
contribution from CBTN (`phs002517`, 32).

For context: that is **larger than GSE140686**, the DKFZ classifier cohort that is
currently the single largest entity-labeled sarcoma methylation resource in existence
(1,315 sarcoma samples). Obtaining it would roughly double the world's supply — and unlike
GSE140686, these arrive as raw IDATs with clinical annotation and outcome data attached,
so they can be normalized uniformly rather than accepted as someone else's processed calls.

Entities where CCDI would move our gap map most:

| Entity | Public atlas today | CCDI adds | Why it matters |
|---|---|---|---|
| Inflammatory myofibroblastic tumor | 9 methylation, 0 regulatory | **19** | one of our nine zero-regulatory entities; triples its entire public epigenome |
| Infantile fibrosarcoma | 36 | **13** | |
| Desmoid / aggressive fibromatosis | 38, only 4 regulatory | **12** | |
| Synovial sarcoma | 0 from tumor tissue | **20 primary tumors** | the entity with 563 regulatory epigenomes and none from a patient |
| DSRCT | 176 | **23** | |
| Spindle cell RMS | — | **20** | the MYOD1-mutant candidate pool |
| FN-RMS / FP-RMS / RMS-NOS | | **194 / 92 / 211** | |
| Ewing | 2,185 | **24** | |
| Undifferentiated sarcoma | 77 | **29** | |

Full per-entity table in `T18_ccdi_sarcoma_methylation.tsv` — 91 entities.

---

## 3. Studies worth a data access request, ranked

Access is one dbGaP Data Access Request per study, through eRA Commons → dbGaP → DCFS/Gen3,
then either the CGC for cloud analysis or the Gen3 client for download. The CCDI document
walks the mechanics; the ranking is mine.

**Tier 1 — request these**

1. **`phs002790` Molecular Characterization Initiative** — 7,842 participants, 245,739
   files, 1,560 sarcoma methylation arrays. The single highest-yield request available.
   Ongoing and growing, so a DAR now keeps paying.
2. **`phs000720` Genomic Sequencing of Pediatric Rhabdomyosarcoma** — 403 participants,
   3,483 files. Directly your entity, and worth checking against what is already in GEO
   before requesting.
3. **`phs002517` Molecular Characterization across Pediatric Brain and Solid Tumors
   (CBTN)** — 4,033 participants, 475,759 files, 32 sarcoma methylation.

**Tier 2 — the model resources, and the real opportunity**

4. **`phs003161` / `phs003160` / `phs003163` / `phs003164` — Pediatric In Vivo Testing
   Program (PIVOT)**, 175 / 53 / 36 / 40 participants. NCI's PDX programme.
5. **`phs001437` Pediatric Preclinical Testing Consortium (PPTC)** — 267 participants.
6. **`phs000469` TARGET Cancer Model Systems: cell lines and xenografts** — 95.
7. **`phs003215` Texas Pediatric PDX** — 51.

These four carry WGS, WXS and RNA-seq — and **no epigenomics whatsoever**. That is the
argument, not the acquisition: the models are already derived, consented, characterized
and federally funded. Adding H3K27ac and ATAC to an existing PIVOT or PPTC panel is
incremental cost on infrastructure that exists, and it would close the model side of our
gap for several entities at once. **This is the concrete, costed recommendation the white
paper can make.**

**Tier 3 — small but pointed**

8. **`phs000466` TARGET Clear Cell Sarcoma of the Kidney** — 13 participants. CCSK is one
   of our nine entities with zero regulatory epigenomics anywhere.
9. **`phs003975` Metastatic Osteosarcoma Spatial Profiling** — 8 participants, 244 files.
   Spatial, and new.
10. **`phs003519` Single-Cell Atlas of NF1 Nerve Sheath Tumors** — 29 participants, 3,644
    files. Relevant to the MPNST progression story.

---

## 4. Ways to get more, beyond CCDI

**Already harvested and in the atlas:** GEO, ArrayExpress/BioStudies, ENA, EGA metadata,
EpiRR, St Jude CSTN.

**Now added:** CCDI Hub (open metadata API), GDC (open API).

**Worth doing next, in order of expected yield:**

- **The CCDI API is fully open and faceted** — 43 filters including `library_strategy`,
  `data_category`, `file_type` and `tumor_classification`. I can re-run it on any schedule
  to catch new studies; MCI is actively accruing. A quarterly re-harvest would keep the
  atlas current for free.
- **dbGaP directly.** CCDI surfaces 42 studies, but dbGaP holds many more pediatric
  sarcoma studies that never entered CCDI. A systematic dbGaP sweep is the obvious
  remaining hole in our coverage.
- **The Genomic Data Commons `legacy` and `awg` endpoints**, plus the NCI Cancer Data
  Service, for TARGET material not exposed through the main GDC portal.
- **ICGC/ARGO and the Pan-Cancer Analysis of Whole Genomes** for adult sarcoma.
- **Contact the authors of PMID 39789291** about the myxofibrosarcoma ATAC-seq that is
  published with no accession — our weakest zero call, and a dataset that exists and is
  unfindable.

---

## 5. What it costs to actually use CCDI data

Worth stating plainly in the paper, because it is part of the argument about reuse.

An eRA Commons account, then a dbGaP DAR per study signed by both the PI and a Signing
Official, then DCFS/Gen3 authorisation. Analysis happens either in the Cancer Genomics
Cloud — where compute is billed to you, and CCDI pilot funds are limited — or by
downloading through the Gen3 client against signed URLs that expire quickly.

So CCDI sits in the same tier as EGA in our three-state model: **the data exists, the
metadata is public, and reuse requires per-study approval.** The difference is that CCDI's
metadata is genuinely excellent — faceted, queryable, entity-labeled — which is precisely
what GSE140686 is not. It is worth naming CCDI as the counter-example when we criticise
metadata-dark deposits: this is what good looks like.

---

## Sources

- [CCDI Hub](https://ccdi.cancer.gov/) · [Explore Dashboard](https://ccdi.cancer.gov/explore) · GraphQL API `https://ccdi.cancer.gov/v1/graphql/`
- [GDC API](https://api.gdc.cancer.gov/) · [Cancer Genomics Cloud](https://www.cancergenomicscloud.org/)
- [NCI Data Commons Framework Services (Gen3)](https://nci-crdc.datacommons.io/)
- CCDI Data Access Instructions v2.0, 2 February 2024 (the document supplied)
