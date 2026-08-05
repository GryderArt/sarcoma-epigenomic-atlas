# Data dictionary

All tables are tab-separated with a single header row, UTF-8, no embedded newlines.
Large tables are gzipped; the pipeline reads either form transparently (`pipeline/_paths.py`).

| Table | Size | Rows | What it is |
|---|---|---|---|
| `GSE140686_recovered_diagnoses.tsv` | 282 KB | 1,505 rows | The same 1,505 samples mapped onto atlas entities, with entity_kind (sarcoma / benign / non-sarcoma), the evidence level for each subtype call, and sample type from the manifestation field. |
| `GSE140686_sample_keys.tsv` | 196 KB | 1,505 rows | The recovered join key: GSM to the REFERENCE_SAMPLE / VALIDATION_SAMPLE id that GEO hides in !Sample_description, with the diagnosis and methylation class from the paper's supplement. |
| `T11_h3k27ac_fidelity.tsv` | 4 KB | 46 rows | H3K27ac fidelity analysis: how well models recapitulate patient tumours, quantified on a common reference peak set with genome builds harmonised and comparisons restricted within study to control for batch. |
| `samples/T4_samples_atomic.tsv.gz` | 2.0 MB | 85,698 rows | THE ATLAS. One row per GEO sample (GSM), 85,698 rows. Disease, age class, matched model, donor key, RRID, sample type, assay class, epigenetic target, antibody, duplicate flags, organism, platform, free-text title/source/characteristics, series title, PubMed, date, contact institute, and the evidence fields subtype_evidence / entity_kind / rms_call_basis. |
| `T11b_fidelity_by_study.tsv` | 0 KB | 7 rows | The batch-effect check behind T11 — same-study versus cross-study correlation. |
| `T12_ebi_all.tsv` | 2.3 MB | 9,589 rows | Full EBI sweep: ArrayExpress/BioStudies, ENA, EGA and EpiRR. Includes GEO mirrors, flagged as such. |
| `T12b_ebi_sarcoma_epigenomic.tsv` | 58 KB | 290 rows | The subset whose title is both sarcoma and epigenomic. |
| `T12c_ebi_unique.tsv` | 4 KB | 23 rows | The 23 studies a GEO-only search can never return (E-MTAB, PRJEB, PRJDB). |
| `T13_incidence_vs_data.tsv` | 10 KB | 45 rows | THE BURDEN JOIN. Curated US incidence per entity with its basis and low/high bracket, against epigenomic coverage: regulatory, methylation, H3K27ac, accessibility, 3D, primary-tumour and patient-derived counts, study concentration, and the adversarial verification caveat for every zero. |
| `T14_ega_sarcoma_datasets.tsv` | 328 KB | 560 rows | All 560 sarcoma datasets in EGA with sample counts, technologies, governing policy and data-access committee. All controlled access. |
| `T15_controlled_access.tsv` | 17 KB | 79 rows | The epigenomic subset of EGA plus St Jude, mapped onto atlas entities and split regulatory versus methylation. |
| `T16_stjude_cstn_inventory.tsv` | 11 KB | 50 rows | St Jude Childhood Solid Tumor Network resource inventory: accessions, repository, access model, assay, sample counts, diseases, PMIDs. |
| `T17_stjude_opdx_models.tsv` | 95 KB | 389 rows | 389 St Jude O-PDX and cell models with diagnosis, sample type, assays available and source publication. |
| `T1_taxonomy.tsv` | 63 KB | 142 rows | Sarcoma taxonomy. One row per subtype: WHO 2020 category, defining lesion, key epigenetic mechanism, markers, paediatric relevance, reference PMID. |
| `T2_models.tsv` | 160 KB | 550 rows | Cell line and PDX catalogue. 550 models with aliases, type, disease, subtype or fusion, key alterations, RRID/Cellosaurus, source repository, origin PMID and a problematic_flag for the 95 with identity problems. |
| `T3_entity_counts.tsv` | 5 KB | 72 rows | Per-entity roll-up: epigenomic samples by assay family, distinct models, distinct PDX, patient epigenomes, and H3K27ac / accessibility / 3D / methylation counts split by model versus patient. |
| `T9b_rms_subtype_log.tsv` | 203 KB | 1,705 rows | Every RMS subtype reassignment with the rule that fired: mutant MYOD1 named, PAX fusion named, curated model identity, explicit fusion-negative wording, driver genotype annotated, or systematically-genotyped series. |
| `T9c_purity_log.tsv` | 54 KB | 288 rows | Samples removed from the bone-tumour bin. H3F3A/H3F3B capture pulls in the paediatric glioma residues K27M and G34R/V plus plant and mouse H3.3 model systems; GCTB uses G34W/L and chondroblastoma K36M. |
| `ega_datasets_raw.json` | 13.1 MB |  |  |

## `incidence/`

The curated evidence behind every rate in T13 — extracted rates with their source, population denominators, and the adult and paediatric literature sweeps. Kept separate so a reviewer can check an anchor without reading the pipeline.

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
| `age_class` | whether the ENTITY is paediatric-predominant, adult-predominant or both. A property of the disease, not the sample. |
| `model_matched` / `donor_key` / `model_rrid` | the model this sample came from, its de-duplication key, and its Cellosaurus RRID |
| `sample_type` | `primary_tumor`, `metastasis`, `recurrence`, `PDX`, `organoid`, `cell_line`, `xenograft(CDX)`, `mouse_model`, `normal/reference`, `unspecified` |
| `assay_class` | normalised assay, e.g. `ChIP-seq`, `CUT&RUN`, `ATAC-seq`, `Hi-C`, `Methyl-array` |
| `epi_target_norm` | normalised target: a histone mark, a factor, `5mC`, or `input/none` for controls |
| `is_duplicate` / `duplicate_of_gsm` / `series_redundant` | re-deposits of the same data, flagged rather than deleted |
| `subtype_evidence` | for GSE140686: whether the subtype came from a stated fusion, the methylation class, or histology alone |
| `rms_call_basis` | for RMS: which rule assigned the subtype |
| `disease_source` | how the entity was determined |

## Counting conventions

**Regulatory** = ChIP-seq, ChIP-exo, ChIP-chip, CUT&RUN, CUT&Tag, ATAC-seq, scATAC-seq, DNase-seq, FAIRE-seq, MNase-seq, Hi-C, HiChIP, Micro-C, Capture-HiC, ChIA-PET, 4C-seq — **excluding input and IgG controls**, which are not profiles of anything.

**DNA methylation** = WGBS, RRBS, methylation array, MeDIP/hMeDIP, bisulfite PCR.

**Patient-derived** = tumour tissue, metastasis, recurrence, PDX or patient-derived organoid. Reported separately from tumour tissue alone, because an organoid grown from a patient's tumour is not a decades-old cell line.

**Epigenomic totals** include input controls and duplicates are removed first (`is_duplicate != Y`).
