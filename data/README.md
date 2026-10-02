# Data dictionary

All tables are tab-separated with a single header row, UTF-8, no embedded newlines.
Large tables are gzipped; the pipeline reads either form transparently (`pipeline/_paths.py`).

| Table | Size | Rows | What it is |
|---|---|---|---|
| `GSE140686_recovered_diagnoses.tsv` | 288 KB | 1,505 rows | The same 1,505 samples mapped onto atlas entities, with entity_kind (sarcoma / benign / non-sarcoma), the evidence level for each subtype call, and sample type from the manifestation field. |
| `GSE140686_sample_keys.tsv` | 199 KB | 1,505 rows | The recovered join key: GSM to the REFERENCE_SAMPLE / VALIDATION_SAMPLE id that GEO hides in !Sample_description, with the diagnosis and methylation class from the paper's supplement. |
| `T11_h3k27ac_fidelity.tsv` | 4 KB | 46 rows | H3K27ac fidelity analysis: how well models recapitulate patient tumors, quantified on a common reference peak set with genome builds harmonised and comparisons restricted within study to control for batch. |
| `samples/T4_samples_atomic.tsv.gz` | 2.0 MB | 92,138 rows | THE ATLAS. One row per GEO sample (GSM), 92,138 rows. Disease, age class, matched model, donor key, RRID, sample type, assay class, epigenetic target, antibody, duplicate flags, organism, platform, free-text title/source/characteristics, series title, PubMed, date, contact institute, and the evidence fields subtype_evidence / entity_kind / rms_call_basis. |
| `T11b_fidelity_by_study.tsv` | 1 KB | 7 rows | The batch-effect check behind T11 — same-study versus cross-study correlation. |
| `T12_ebi_all.tsv` | 2.3 MB | 9,589 rows | Full EBI sweep: ArrayExpress/BioStudies, ENA, EGA and EpiRR. Includes GEO mirrors, flagged as such. |
| `T12b_ebi_sarcoma_epigenomic.tsv` | 60 KB | 290 rows | The subset whose title is both sarcoma and epigenomic. |
| `T12c_ebi_unique.tsv` | 5 KB | 23 rows | The 23 studies a GEO-only search can never return (E-MTAB, PRJEB, PRJDB). |
| `T13_incidence_vs_data.tsv` | 11 KB | 45 rows | THE BURDEN JOIN. Curated US incidence per entity with its basis and low/high bracket, against epigenomic coverage: regulatory, methylation, H3K27ac, accessibility, 3D, primary-tumor and patient-derived counts, study concentration, and the adversarial verification caveat for every zero. |
| `T14_ega_sarcoma_datasets.tsv` | 335 KB | 560 rows | All 560 sarcoma datasets in EGA with sample counts, technologies, governing policy and data-access committee. All controlled access. |
| `T15_controlled_access.tsv` | 23 KB | 94 rows | The epigenomic subset of EGA plus St Jude, mapped onto atlas entities and split regulatory versus methylation. |
| `T16_stjude_cstn_inventory.tsv` | 12 KB | 50 rows | St Jude Childhood Solid Tumor Network resource inventory: accessions, repository, access model, assay, sample counts, diseases, PMIDs. |
| `T17_stjude_opdx_models.tsv` | 98 KB | 389 rows | 389 St Jude O-PDX and cell models with diagnosis, sample type, assays available and source publication. |
| `T1_taxonomy.tsv` | 64 KB | 142 rows | Sarcoma taxonomy. One row per subtype: WHO 2020 category, defining lesion, key epigenetic mechanism, markers, pediatric relevance, reference PMID. |
| `T20_ccdi_studies.tsv` | 4 KB | 42 rows | Every study in the NCI Childhood Cancer Data Initiative, with participant, sample and file counts. |
| `T21_ccdi_samples.tsv.gz` | 233 KB | 23,143 rows | The sarcoma-cohort samples. All 70,820 CCDI samples are fetched and the participant join runs over all of them — that is how the cohort is defined — but only sarcoma rows are kept. `sarcoma_call_basis` records which route made each call. All-disease totals survive in `T26`. |
| `T22_ccdi_files.tsv.gz` | 463 KB | 9,187 rows | The EPIGENOMIC files belonging to a sarcoma-cohort participant. All 319,084 are fetched and examined; the 97% that are WGS, WXS, RNA-seq, panels and their indexes are dropped — they say nothing about chromatin and cost 23 MB. |
| `T22b_ccdi_file_census.tsv` | 6 KB | 158 rows | `library_strategy` × `file_type` × `data_category` across ALL 319,084 files. This is what keeps "not one regulatory epigenomic file among 319,084" checkable without shipping the rows. |
| `T23_ccdi_participants.tsv.gz` | 128 KB | 8,727 rows | The sarcoma-cohort participants, with diagnosis, ICD-O category, anatomic site and survival status. The authoritative diagnosis source — sample rows are blank 76% of the time. |
| `T26_ccdi_disease_census.tsv` | 4 KB | 86 rows | Samples by ICD-O category and tumor status across ALL of CCDI, so the all-disease denominators — and the finding that CCDI's 91 ATAC-seq and 108 bisulfite-seq participants are entirely leukaemia — stay checkable after the sarcoma trim. |
| `T24_ccdi_entity_counts.tsv` | 2 KB | 36 rows | CCDI mapped onto atlas entities: participants, samples, tumor/normal split, and the methylation-array subset. All CONTROLLED access. |
| `T25_ccdi_unmapped.tsv` | 22 KB | 324 rows | CCDI diagnoses the mapping could not place, with the reason — mostly correctly-excluded non-sarcomas (neuroblastoma, Wilms) plus MCI's `see diagnosis_comment` placeholder. The residual is visible, not silently dropped. |
| `T27_access_routes.tsv` | 1 KB | 3 rows | Resources with a formal request route but no archive accession — St Jude's COMET methylation project and the St Jude Cloud CSTN dataset. They sit behind the same wall as EGA and dbGaP, and are counted with it; the missing accession is what makes them the weakest case within it. |
| `T28_stjude_viz_tracks.tsv` | 231 KB | 398 rows | The CSTN epigenetic browsers enumerated at sample × assay resolution, recovered from the ProteinPaint manifests the pages ship inline. 28 models, 317 regulatory tracks, the same 10-mark panel plus WGBS on every one. **Diagnoses are the St Jude data administrator's**, supplied by email and transcribed verbatim in stage 15; the fusion partner comes from the CSTN portal model table. `entity_call_basis`, `subtype_evidence`, `fusion_source` and `availability` carry the provenance per row — including the one model withdrawn from CSTN and the one model the administrator's list omits. NOT added to any total: the same material is already counted once under EGA. See `docs/methods.md` §3b and figure F26. |
| `T29_rna_assay_corrections.tsv` | 135 KB | 723 rows | Correction #15: samples whose `library_strategy` is RNA-based but which had inherited a DNA-based `assay_class` from their series context. One row per reassignment, with the old and new call and the evidence that forced it. |
| `T30_delta_samples.tsv.gz` | 2.4 MB | 89,475 rows | The most recent delta harvest: GEO series released after the previous census close, classified by the same `04_classify` rules as the rest of the atlas. Merged into T4 by `32_merge_delta.py`; kept so a re-harvest can be audited separately from the table it was folded into. |
| `T31_mixed_cohort_attribution.tsv` | 3 KB | 19 rows | Controlled-access cohorts deposited above the subtype, mapped to the named entities they bear on. `basis` is *verified* where the St Jude data administrator's per-model list names the entity and *parent* where the deposit states only the parent diagnosis. These are listed on an entity's page but added to no entity's totals: the same samples would otherwise be counted once per subtype. Correction #17. |
| `T32_administrator_corrections.tsv` | 35 KB | 265 rows | Every sample whose diagnosis was reassigned from the St Jude data administrator's per-model list, and every barcode dial-out demoted from epigenomic, with the model base and the basis for each. Corrections #18 and #18b. |
| `T2_models.tsv` | 164 KB | 550 rows | Cell line and PDX catalog. 550 models with aliases, type, disease, subtype or fusion, key alterations, RRID/Cellosaurus, source repository, origin PMID and a problematic_flag for the 95 with identity problems. |
| `T3_entity_counts.tsv` | 6 KB | 79 rows | Per-entity roll-up: epigenomic samples by assay family, distinct models, distinct PDX, patient epigenomes, and H3K27ac / accessibility / 3D / methylation counts split by model versus patient. |
| `T9b_rms_subtype_log.tsv` | 41 KB | 326 rows | Every RMS subtype reassignment with the rule that fired: mutant MYOD1 named, PAX fusion named, curated model identity, explicit fusion-negative wording, driver genotype annotated, or systematically-genotyped series. |
| `T9c_purity_log.tsv` | 60 KB | 296 rows | Samples removed from the bone-tumor bin. H3F3A/H3F3B capture pulls in the pediatric glioma residues K27M and G34R/V plus plant and mouse H3.3 model systems; GCTB uses G34W/L and chondroblastoma K36M. |
| `ega_datasets_raw.json` | 13.1 MB |  |  |

## `incidence/`

The curated evidence behind every rate in T13 — extracted rates with their source, population denominators, and the adult and pediatric literature sweeps. Kept separate so a reviewer can check an anchor without reading the pipeline.

| File | What it holds |
|---|---|
| `incidence/ADULT_incidence_A.tsv` | Adult sarcoma incidence, sweep A — registry and review series. |
| `incidence/ADULT_incidence_B.tsv` | Adult sarcoma incidence, sweep B — entity-specific literature. |
| `incidence/BONE_RMS_incidence.tsv` | Bone sarcoma and rhabdomyosarcoma rates. |
| `incidence/STS_incidence.tsv` | Soft tissue sarcoma rates by entity. |
| `incidence/population_denominators.tsv` | UN WPP 2024 and US census denominators used to convert rates to counts. |
| `incidence/reference_rates.tsv` | All extracted per-million rates with source and age scope. |
| `incidence/worldwide_estimates.tsv` | Global estimates where a US-only figure would mislead. |

## Key columns in T4 (the atomic table)

| Column | Meaning |
|---|---|
| `disease` | the entity bin this sample is assigned to |
| `entity_kind` | `sarcoma`, `benign` or `nonsarcoma` — controls and mimics are kept, not silently dropped |
| `age_class` | whether the ENTITY is pediatric-predominant, adult-predominant or both. A property of the disease, not the sample. |
| `model_matched` / `donor_key` / `model_rrid` | the model this sample came from, its de-duplication key, and its Cellosaurus RRID |
| `sample_type` | `primary_tumor`, `metastasis`, `recurrence`, `PDX`, `organoid`, `cell_line`, `xenograft(CDX)`, `mouse_model`, `normal/reference`, `unspecified` |
| `assay_class` | normalized assay, e.g. `ChIP-seq`, `CUT&RUN`, `ATAC-seq`, `Hi-C`, `Methyl-array` |
| `epi_target_norm` | normalized target: a histone mark, a factor, `5mC`, or `input/none` for controls |
| `is_duplicate` / `duplicate_of_gsm` / `series_redundant` | re-deposits of the same data, flagged rather than deleted |
| `subtype_evidence` | for GSE140686: whether the subtype came from a stated fusion, the methylation class, or histology alone |
| `rms_call_basis` | for RMS: which rule assigned the subtype |
| `disease_source` | how the entity was determined |

## A CCDI counting trap worth knowing

76% of CCDI **sample** rows carry no `diagnosis_category` of their own. The diagnosis lives
on the **participant**, and the Explore Dashboard's sample counts come from a server-side
join. Flagging sarcoma from the sample row alone undercounts roughly 17-fold — 1,350
against a true 23,143. `T21` therefore carries `sarcoma_call_basis` showing which route
made each call.

A second trap: CCDI returns file-level ids wrapped in brackets (`[1794211]`) while
sample-level ids are bare (`003-1`), and several studies use different id spaces in the two
tables entirely. Joins between `T21` and `T22` are unreliable at sample level and are done
participant-to-file instead.

## Counting conventions

**Regulatory** = ChIP-seq, ChIP-exo, ChIP-chip, CUT&RUN, CUT&Tag, ATAC-seq, scATAC-seq, DNase-seq, FAIRE-seq, MNase-seq, Hi-C, HiChIP, Micro-C, Capture-HiC, ChIA-PET, 4C-seq — **excluding input and IgG controls**, which are not profiles of anything.

**DNA methylation** = WGBS, RRBS, methylation array, MeDIP/hMeDIP, bisulfite PCR.

**Patient-derived** = tumor tissue, metastasis, recurrence, PDX or patient-derived organoid. Reported separately from tumor tissue alone, because an organoid grown from a patient's tumor is not a decades-old cell line.

**Epigenomic totals** include input controls and duplicates are removed first (`is_duplicate != Y`).
