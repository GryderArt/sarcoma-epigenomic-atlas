# Adult sarcoma burden, the EBI archives, and a metadata failure worth naming

Working notes for *Critical assessment of existing datasets* — SASS Epigenetics Working Group.
Atlas version 3 (all ages). Every number below is reproducible from
`SarcomaEpigenomicDataAtlas_v3.xlsx`, tab `T13_Incidence_vs_Data`.

---

## 1. Kevin Jones's point, with numbers behind it

**The headline.** Adult sarcoma carries roughly eight times the case burden of pediatric
sarcoma (89.4% of US cases) and 30% of the regulatory epigenomic data. Per new case diagnosed,
pediatric-predominant entities carry **~24× more** ChIP-seq / CUT&RUN / ATAC / Hi-C data
than adult-predominant entities.

**National figures (verified against primary sources).**

| | Value | Source |
|---|---|---|
| US sarcoma cases, 2026 | **18,020** (4,110 bone + 13,910 soft tissue) | Siegel et al., *Cancer Statistics 2026*, PMID 41528114 |
| Aged 0–19 | **~1,900/yr** site-coded; ~1,700/yr by ICCC histology | SEER 21 age distribution; Siegel DA et al., *JNCI* 2023, PMID 37433078 |
| Adult fraction | **89.4%** | derived from the two rows above |
| Sarcoma as a share of all cancers, age 0–19 | **11.5%** (ICCC groups VIII + IX) | Siegel DA et al., PMID 37433078 |
| Sarcoma as a share of all cancers, age 20+ | **0.8%** site-coded, ~1.1% histology-based | derived |

Two corrections to figures that circulate widely and that we should not repeat:

- Pediatric sarcoma is **~1,700–1,900/yr**, not ~1,500. The commonly quoted "18,000 total,
  1,500 pediatric" pair is internally inconsistent — it implies 92% adult, not 89%.
- Sarcoma is **~11.5%** of pediatric malignancy, not the often-quoted 15–20%. I could find
  no source stratifying ICCC VIII/IX for ages 15–19 alone, which may be where the higher
  number came from; we should not cite it without one.

**One caveat to state in the paper.** The 18,020 figure is topography-coded, so it excludes
sarcomas arising in viscera and uterus (GIST, uterine leiomyosarcoma) that are coded to
those organs. Extrapolating the French central-review rate of 70.7 sarcomas per million to
the US population gives **~24,000/yr**. The undercount is in the adult direction, so the
gap we are describing is if anything larger than stated.

**Where the asymmetry actually bites.** It is not that adult entities have *less* data —
it is that the data they have cannot answer regulatory questions. **13 of the 45 named
entities in the atlas have no regulatory epigenomics at all** — no ChIP-seq, no CUT&RUN,
no ATAC, no Hi-C, in any sample type, open or controlled. Nine of those 13 have a
published incidence rate, so we can say what the gap costs per patient:

| Entity | US cases/yr | Epigenomic samples | Regulatory |
|---|---|---|---|
| Dermatofibrosarcoma protuberans | 1,400 | 40 | **0** |
| Angiosarcoma | 1,312 | 46 | **0** |
| Myxofibrosarcoma | 820 | 17 | **0** |
| Pleomorphic liposarcoma | 185 | 28 | **0** |
| Inflammatory myofibroblastic tumor | 175 | 9 | **0** |
| Extraskeletal myxoid chondrosarcoma | 75 | 11 | **0** |
| Intimal sarcoma | 60 | 1 | **0** |
| Clear cell sarcoma of the kidney | 20 | 12 | **0** |
| BCOR-CCNB3 sarcoma | 15 | 9 | **0** |

The other four have no published population rate anywhere *and* no regulatory
epigenomics: myoepithelial carcinoma of soft tissue, ossifying fibromyxoid tumor,
angiomatoid fibrous histiocytoma, and atypical fibroxanthoma / pleomorphic dermal sarcoma.
For these we cannot even state what the gap costs — the denominator does not exist either.

**A note on denominators, since every count below depends on it.** An *entity* is one
diagnostic category in the atlas taxonomy: a WHO 2020 diagnosis, except that RMS and
liposarcoma are split to molecular subtype because that is the resolution this section
argues at. There are **45 named entities**. The atlas also carries four *unresolved
residual bins* — "Sarcoma NOS", "RMS-NOS" and two one-off cases, holding samples whose
diagnosis the deposit never states — and one *control tissue* category. Those five are not
diseases and are excluded from every "of N" figure in this memo, in the figures and in the
gap map. They still contribute to sample and record totals, because the samples exist; we
just cannot say what they are. The residual bins are themselves a metadata-quality
measure: 535 epigenomic samples sit in them.

**And the primary-tumor gap you flagged is broader than EHE.** **25 of 45** entities have
never had a single piece of patient-derived material — tumor tissue, PDX or
patient-derived organoid — profiled for active chromatin, accessibility or 3D architecture — including leiomyosarcoma (3,310 US cases/yr),
desmoid (1,740), DFSP (1,400), angiosarcoma (1,312), chondrosarcoma (1,280), Kaposi
sarcoma (1,268) and solitary fibrous tumor (1,200). Synovial sarcoma is the single entity
rescued by widening the definition from tumor tissue to patient-derived material: its 563
regulatory epigenomes are cell lines and mouse models except for one organoid series
(GSE148722), and it still has nothing from tumor tissue.

**Each zero was independently verified.** A separate adversarial sweep of GEO (series and
sample level), SRA across all library strategies, ENCODE, ChIP-Atlas (845,824 experiments),
ArrayExpress/BioStudies, GDC/TARGET, dbGaP and Cellosaurus cell-line rosters, plus PubMed
and Europe PMC full text, tried to falsify each call. Caveats that survived are recorded
per entity in the `verification_caveat` column of T13. Three matter for wording:

1. **Scope angiosarcoma to human.** Canine hemangiosarcoma — the accepted spontaneous
   model — does have CUT&Tag (GSE304509, H3K27ac/H3K4me3/Kla/RNAPII-Ser5, PMID 41767677)
   and ChRO-seq (GSE150705).
2. **Soften myxofibrosarcoma to "no retrievable data."** PMID 39789291 reports ATAC-seq on
   an MFS PDX, but states no accession and no matching series exists in GEO or SRA. Worth
   emailing the authors — this is a dataset that exists and is unfindable.
3. **Say "DNA methylation profiling (predominantly array-based)"** rather than "methylation
   arrays". DFSP has one nanopore-methylation tumor (GSE320108); CCSK has cfDNA RRBS
   (E-MTAB-8770).

A quotable confirmation from the field itself, for BCOR-CCNB3: *"For BCOR-rearranged fusion
proteins, no fusion-specific dataset was available"* — PMID 40841360, *Nat Commun* 2025,
whose authors built BCOR::CCNB3 tumoroids and had to fall back on ENCODE BCOR ChIP-seq.

---

## 2. How we get EBI data — and what it adds

All four EBI resources are reachable programmatically without credentials. I ran them; the
results are in tabs `T12_EBI_all`, `T12b_EBI_sarcoma_epi` and `T12c_EBI_unique`.

| Resource | Endpoint | Result for sarcoma |
|---|---|---|
| ArrayExpress / BioStudies | `https://www.ebi.ac.uk/biostudies/api/v1/arrayexpress/search?query=…` | 7,856 hits; 2,033 native (`E-MTAB`/`E-MEXP`/`E-TABM`), the rest `E-GEOD` GEO mirrors |
| ENA portal | `https://www.ebi.ac.uk/ena/portal/api/search?result=study&query=study_title="*…*"` | 1,732 studies; mostly `PRJNA` SRA-side mirrors of GEO |
| EGA | `https://metadata.ega-archive.org/studies?query=…` | metadata open, **data controlled** — DAC approval per study |
| EpiRR (IHEC registry) | `https://www.ebi.ac.uk/vg/epirr/view/all?format=tsv` (returns JSON) | 2,888 reference epigenomes, **0 sarcoma** |

**The EpiRR result is the finding.** The International Human Epigenome Consortium's
registry of reference epigenomes contains 2,888 entries and not one is a sarcoma. Whatever
else the white paper says about gaps, sarcoma has never been part of the reference
epigenome effort at all.

**What EBI actually adds beyond GEO: 23 studies** that a GEO-only search can never return
(tab `T12c_EBI_unique`). Most are methylation, but three are regulatory and directly
relevant to entities we score as sparse:

- `PRJEB74484` — **ChIP-seq of myxoid liposarcoma PDX models**
- `PRJDB12634` — **Hi-C of sarcoma** and renal cell carcinoma
- `PRJDB8273` — ZNF217 ChIP-seq in the leiomyosarcoma line UT1
- `E-MTAB-9875` — methylation (450K + EPIC) for sarcoma classification
- `E-MTAB-6961` — methylation array profiling of **undifferentiated sarcomas of adults**
- `E-MTAB-8864` — MPNST methylation · `E-MTAB-11031` — chondrosarcoma methylation
- `E-MTAB-6708` — malignant rhabdoid tumor methylation
- `PRJEB104098` — histiocytic sarcoma vs UPS · `PRJEB112272` — iPSC model of MPNST progression

Practically: `E-GEOD-*` and `PRJNA*` are mirrors and can be dropped; `E-MTAB-*`, `PRJEB*`,
`PRJDB*` and EGA are the genuinely additional records. That is the rule the harvester
(`harvest_ebi.py`) applies.


---

## 2b. EGA and St Jude: what exists but cannot be reused

I enumerated the **entire EGA dataset catalog — 21,279 datasets** — via the metadata API
(`https://metadata.ega-archive.org/datasets?limit=200&offset=N`; the apparent `query=`
parameter is silently ignored, so the catalog has to be pulled whole and filtered
locally). Results in tabs `T14_EGA_sarcoma` and `T15_Controlled_epigenomics`.

**560 sarcoma datasets, 22,832 samples, and every single one is controlled access.** Of
those, **65 datasets / 2,775 samples are epigenomic**. Folding in the St Jude holdings
gives **1,207 regulatory epigenomes** (ChIP-seq, ATAC, CUT&RUN, Hi-C, 4C) and **7,464
methylation samples** that exist and cannot be downloaded.

This forces a distinction the paper should make explicitly, because a GEO-only survey
cannot see it:

| State | Meaning | Scale (regulatory samples) |
|---|---|---|
| **Open** | downloadable and reanalysable today | 10,233 |
| **Controlled** | the experiment was done, the metadata is public, the data needs a per-dataset DAA | 1,207 |
| **Absent** | no record in any archive | 9 entities at zero |

The controlled tier is concentrated in a handful of gatekeepers: BC Cancer governs 204 of
the 560 datasets, DKFZ-HIPO 40, ICGC DACO 30, Sanger CGP 50, Princess Máxima 28, and the
St Jude–WashU Pediatric Cancer Genome Project 20.

**The point that matters for the gap map: none of the controlled-access holdings cover the
nine entities with no regulatory epigenomics.** I checked each explicitly — angiosarcoma,
DFSP, myxofibrosarcoma, pleomorphic liposarcoma, IMT, EMC, intimal sarcoma, CCSK and
BCOR-CCNB3 appear nowhere in EGA either. Those zeros survive.

Datasets in the controlled tier that fill entity gaps our open-data map shows as thin:
`EGAD00001004135` (SS18-SSX BAF hijacking, synovial sarcoma, 85), `EGAD00001006253`
(MPNST, 98), `EGAD50000002533` (MPNST iPSC progression model with ATAC-seq, 56),
`EGAD00001015649` (bulk ATAC-seq of 42 malignant rhabdoid tumors), `EGAD00001011820`
(ATAC / Hi-C / 4C, malignant rhabdoid), `EGAD00001005109` (GCTB genomic-epigenomic, 29),
`EGAD00010002571` (epithelioid sarcoma 850K, 32), `EGAD00010002338` (chordoma
methylation, 68) and `EGAD00003133` (Ewing RRBS via ICGC, 86).

**St Jude CSTN** (tabs `T16_StJude_CSTN`, `T17_StJude_OPDX_models`) — the outstanding
non-GEO source from the Aug 3 notes. The disease-cohort browser is an R Shiny app that
returns 403 to non-browser clients, but it is backed by an undocumented JSON API
(`cstn-gateway-prod.azurewebsites.net/p/samplesearch/search`) which enumerates **375 live
models across 20 diagnoses**: osteosarcoma 71, rhabdomyosarcoma 66, Wilms 63,
retinoblastoma 56, neuroblastoma 34, Ewing 27, then a long tail including synovial sarcoma,
DSRCT, rhabdoid tumor, clear cell sarcoma, epithelioid sarcoma, MPNST, GIST and
liposarcoma as singletons. 389 models once the 14 Nature-2017-only models are added back.

Three findings there belong in the paper:

1. **The entire CSTN epigenome is EGA-only.** The 756 ChIP-seq libraries (19 antibodies)
   and the WGBS from Stewart et al., *Cancer Cell* 2018 (**PMID 30146332**) sit in
   `EGAD00001004312`, `EGAD00001006398` and `EGAD00001004315`. **Zero CSTN ChIP-seq and
   zero CSTN WGBS are in GEO** — and they are not in St Jude Cloud either: the CSTN dataset
   there (`SJC-DS-1008`) is 143 samples of WGS/WES/RNA-seq only.
2. **The portal displays epigenetic browsers for six cohorts — Ewing, neuroblastoma,
   osteosarcoma, retinoblastoma, RMS and rare tumors — but only RMS, plus four
   RMS/OS/NB models from the Nature 2017 paper, are deposited anywhere.** Ewing,
   retinoblastoma and rare-tumor ChIP-seq are viewable in the browser with no deposit in
   EGA, GEO or dbGaP that an exhaustive search can find. That is undeposited data, not
   merely controlled data.
3. **COMET**, the ">4,700 sample" methylation project the Data Types page advertises, has
   **no public accession anywhere**, and the portal shows it active for 66 models, all
   rhabdomyosarcoma. Treat the 4,700 as an internal figure, not an accessible resource.

Access, for the record: models are free to academics under an MTA (max 5 samples per
request, no onward sharing, mandated by the MAST consent, NCT01050296). Sequencing data is
a separate agreement — WGS/WES/RNA-seq via the St Jude Cloud DAA; ChIP-seq and WGBS via a
per-`EGAD` data-access agreement signed by the PI *and* an authorised organizational
representative, roughly 15–20 business days. Browsing the chromHMM tracks is open;
downloading the data is not.

---

## 3. The metadata failure we should name in the paper

The single largest sarcoma epigenomic resource in existence is **GSE140686** — 1,505 DNA
methylation arrays from the DKFZ sarcoma classifier (Koelsche et al., *Nat Commun* 2021,
PMID 33479225). Every sample is deposited as `sarcoma classifier reference case N`, with
characteristics reading `tissue: sarcoma`. **No diagnosis appears anywhere in the GEO
record.** The per-case entity label exists only in the paper's Supplementary Data, joined
through an identifier GEO buries in `!Sample_description` — a field most harvesting
pipelines never read.

I recovered the join. Those 1,505 samples resolve to **46 entities**: 1,315 sarcoma, 154
benign or non-neoplastic controls, 36 non-sarcoma. The recovery is in tab
`GSE140686_recovered` and can be reused by anyone.

The consequence is not academic. Before recovery, angiosarcoma, DFSP, IMT, EMC, CCSK,
BCOR-CCNB3, intimal sarcoma, AFH, OFMT and AFX/PDS each appeared to have **zero epigenomic
data of any kind**. All of them in fact have methylation arrays — sitting in a public,
open-access deposit, invisible to every entity-level query. For ten entities, **100% of all
public epigenomic data comes from this one metadata-dark deposit** (figure F22).

This is a concrete, fixable recommendation for the white paper: **entity labels belong in
`Sample_characteristics`, not in a supplementary spreadsheet.** If a consortium deposit
cannot be found by disease, it cannot be reused, and it does not count as shared data.

---

## 4. Corrections applied to the atlas this round

| # | What was wrong | Effect |
|---|---|---|
| 1 | GSE140686 unlabeled | 1,505 samples recovered into 46 entities |
| 2 | RMS subtype rule lost in the v1→v2 rebuild | Re-applied. The Gryder/Yohe primary tumors now read 5 FN-RMS / 3 FP-RMS, matching your figure |
| 3 | Driver-genotype inference over-firing on cell lines | Curated cell line identity now outranks series context; RH4/RH30/RH41/CW9019 no longer swept into FN-RMS, RD/SMS-CTR/RH36 no longer sitting in FP-RMS |
| 4 | H3F3A/H3F3B capture pulling glioma into the bone bin | 288 rows removed (K27M, G34R/V, plant and mouse H3.3 systems). GCTB uses G34W/L; chondroblastoma uses K36M |
| 5 | ChIP inputs counted as profiles | Input and IgG controls excluded from every "regulatory" count |

Atlas totals after correction: **85,698 in-scope samples · 16,787 unique epigenomic ·
4,157 series · 550 cataloged models (95 with identity problems)**.

---

## 5. Figures

All seven are single-page PDFs, Arial, TrueType-embedded — every label is live text in
Illustrator.

| Figure | Point it makes |
|---|---|
| **F18** | Adult sarcoma: ~5× the burden of these entities, a fraction of the data, 24× gap per case |
| **F23** | Sarcoma is 1% of adult cancer and 11% of childhood cancer — 11% of cases carry 70% of the data |
| **F19** | Burden vs regulatory epigenomics, log–log, with the zero floor called out in red |
| **F20** | The methylation-only tier, ordered by incidence (the nine with a published rate) |
| **F21** | 25 of 45 entities have never had patient-derived material profiled — the version of your EHE figure that covers the whole field |
| **F24** | Open vs controlled vs absent — separating unmeasured biology from data that exists and cannot be reused |
| **F22** | Single-study dependence, and how much of it is one metadata-dark deposit |

---

## Sources

- [Siegel RL et al. Cancer statistics, 2026. *CA Cancer J Clin* (PMID 41528114)](https://pubmed.ncbi.nlm.nih.gov/41528114/)
- [Siegel DA et al. Pediatric cancer in the United States, 2003–2019. *JNCI* 2023 (PMID 37433078)](https://pubmed.ncbi.nlm.nih.gov/37433078/)
- [Koelsche C et al. Sarcoma classification by DNA methylation profiling. *Nat Commun* 2021 (PMID 33479225)](https://pubmed.ncbi.nlm.nih.gov/33479225/) · [GSE140686](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140686)
- [Wagner MJ et al. Incidence and presenting characteristics of angiosarcoma in the US, 2001–2020. *JAMA Netw Open* 2024 (PMID 38607625)](https://pubmed.ncbi.nlm.nih.gov/38607625/)
- [de Pinieux G et al. Nationwide incidence of sarcomas using an expert pathology review network. *PLoS One* 2021 (PMID 33630918)](https://pubmed.ncbi.nlm.nih.gov/33630918/)
- [Ray-Coquard I et al. Concordance between initial diagnosis and centralized expert review. *Ann Oncol* 2012 (PMID 22331640)](https://pubmed.ncbi.nlm.nih.gov/22331640/) · [Lurkin A et al. *BMC Cancer* 2010 (PMID 20403160)](https://pubmed.ncbi.nlm.nih.gov/20403160/)
- [Criscione VD & Weinstock MA. DFSP epidemiology. *J Am Acad Dermatol* 2007 (PMID 17141362)](https://pubmed.ncbi.nlm.nih.gov/17141362/) · [Kreicher KL et al. *Dermatol Surg* 2016 (PMID 26730971)](https://pubmed.ncbi.nlm.nih.gov/26730971/)
- [Ma J et al. Incidence of undifferentiated pleomorphic sarcoma in the United States. *Sarcoma* 2024 (PMID 39502684)](https://pubmed.ncbi.nlm.nih.gov/39502684/)
- [O'Donnell JS et al. *npj Precis Oncol* 2025 (PMID 39789291)](https://pubmed.ncbi.nlm.nih.gov/39789291/) — the unretrievable myxofibrosarcoma ATAC-seq
- [BCOR::CCNB3 tumoroids, *Nat Commun* 2025 (PMID 40841360)](https://pubmed.ncbi.nlm.nih.gov/40841360/)
- [Stewart E et al. Orthotopic patient-derived xenografts of pediatric solid tumors. *Nature* 2017 (PMID 28854174)](https://pubmed.ncbi.nlm.nih.gov/28854174/)
- [Stewart E et al. Integrated genomic, epigenomic and proteomic analyses of rhabdomyosarcoma. *Cancer Cell* 2018 (PMID 30146332)](https://pubmed.ncbi.nlm.nih.gov/30146332/)
- [St Jude Childhood Solid Tumor Network portal](https://cstn.stjude.cloud/) · [EGA dataset catalog](https://ega-archive.org/datasets/)
- [EpiRR / IHEC reference epigenome registry](https://www.ebi.ac.uk/vg/epirr/) · [ArrayExpress](https://www.ebi.ac.uk/biostudies/arrayexpress) · [ENA portal API](https://www.ebi.ac.uk/ena/portal/api/) · [EGA](https://ega-archive.org/)
