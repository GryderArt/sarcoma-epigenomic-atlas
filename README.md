# Sarcoma Epigenomic Data Atlas

A sample-level census of every publicly recorded epigenomic experiment on sarcoma — all
ages, bone and soft tissue — built to answer one question: **where are the holes?**

Assembled for the *Critical assessment of existing datasets* section of the SASS
Epigenetics Working Group white paper.

**[→ Open the interactive gap map](https://gryderart.github.io/sarcoma-epigenomic-atlas/)**

---

## What's here

| | |
|---|---|
| **16,394** | sarcoma epigenomic samples, de-duplicated |
| **4,416** | assay × model × study records — the "H3K27ac ChIP-seq for RH4" grain |
| **943** | studies across GEO, ArrayExpress, ENA, EGA, St Jude CSTN and NCI CCDI |
| **45** | named diagnostic entities (plus 4 unresolved bins and 1 control tissue) |
| **550** | catalogued cell lines and PDX models, 95 flagged for identity problems |

## Five findings the atlas was built to support

**1. Adult sarcoma carries ~8× the burden and a third of the data.** 89.4% of the 18,020
US sarcoma cases per year occur at age 20+, but paediatric-predominant entities hold 70%
of all regulatory epigenomic data ever generated — **~24× more per new case diagnosed**.

**2. Thirteen of 45 named entities have no regulatory epigenomics at all.** No ChIP-seq,
no CUT&RUN, no ATAC, no Hi-C, in any sample type, open or controlled. Nine have a
published incidence rate, so the gap can be costed per patient: dermatofibrosarcoma
protuberans (1,400 US cases/yr), angiosarcoma (1,312), myxofibrosarcoma (820), pleomorphic
liposarcoma, inflammatory myofibroblastic tumour, extraskeletal myxoid chondrosarcoma,
intimal sarcoma, clear cell sarcoma of the kidney, BCOR-CCNB3 sarcoma. Every one of these
zeros was independently verified against GEO, SRA, ENCODE, ChIP-Atlas, ArrayExpress, GDC,
dbGaP and Europe PMC full text; surviving caveats are recorded per entity in
[`data/T13_incidence_vs_data.tsv`](data/T13_incidence_vs_data.tsv).

**3. Twenty-five of 45 entities have never had a patient's tumour profiled** for active
chromatin, accessibility or 3D architecture — counting tumour tissue, PDX and
patient-derived organoids together. Synovial sarcoma is the sharpest case: 563 regulatory
epigenomes exist and all but one organoid series are cell lines and mouse models.

**4. The largest sarcoma epigenomic resource in existence is unfindable by disease.**
[GSE140686](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140686) — the DKFZ
sarcoma methylation classifier, 1,505 arrays — is deposited with every sample labelled
`sarcoma classifier reference case N` and characteristics reading `tissue: sarcoma`. No
diagnosis appears anywhere in the GEO record. **This repository recovers the mapping**
(see [`data/GSE140686_recovered_diagnoses.tsv`](data/GSE140686_recovered_diagnoses.tsv)):
1,505 samples across 46 entities. Ten entities depend on that one metadata-dark deposit
for 100% of their public epigenome.

**5. The US national cancer data infrastructure has no regulatory epigenomics for sarcoma
at all.** The NCI Childhood Cancer Data Initiative — 42 studies, 61,854 participants,
1,259,864 files — has **no ChIP-seq in its assay vocabulary**, for any disease. Its 91
ATAC-seq and 108 bisulfite-seq participants are entirely leukaemia. Across all **319,084
files belonging to a sarcoma-cohort participant there is not one regulatory epigenomic
file**. The GDC is the same: no ChIP-seq in the vocabulary, and its 410 ATAC-seq files
belong to 23 TCGA cohorts, none of them SARC or any TARGET sarcoma project.

What CCDI does hold is worth having: **1,592 sarcoma participants with raw methylation
IDATs, 1,382 of them primary tumours** — larger than GSE140686, and with clinical
annotation attached. `data/T24_ccdi_entity_counts.tsv`.

The CCDI tables ship trimmed to what an epigenomics atlas needs: the sarcoma cohort, and
within it the epigenomic files. The full harvest is run in every case — that is how the
cohort and the absence are established — but 163,547 WGS and 78,836 WXS rows are not
carried. `T22b` and `T26` keep the full denominators, so every claim above remains
checkable from this repository.

**6. A third of the field's regulatory epigenomics is locked, not missing.** 560 sarcoma
datasets in EGA, all controlled access — 1,207 regulatory epigenomes behind data-access
committees. The entire St Jude CSTN epigenome (756 ChIP-seq libraries plus WGBS) is
EGA-only; none of it is in GEO, and it is not in St Jude Cloud either. Crucially, **none
of the controlled holdings cover the thirteen entities with nothing** — those zeros
survive.

And the model resources are the opportunity: PIVOT, PPTC, TARGET Model Systems and the
Texas PDX bank hold roughly 500 paediatric sarcoma models with WGS, WXS and RNA-seq and
**no epigenomics whatsoever**. Those models are already derived, consented and federally
funded — adding H3K27ac and ATAC to an existing panel is incremental cost on
infrastructure that exists.

Separating *open*, *controlled* and *absent* is the point of
[`figures/F24_open_vs_controlled_vs_absent.pdf`](figures/F24_open_vs_controlled_vs_absent.pdf).
A GEO-only survey cannot tell the three apart.

## What counts as an "entity"

One diagnostic category in the atlas taxonomy: a WHO 2020 diagnosis, except that
rhabdomyosarcoma and liposarcoma are split to molecular subtype (FP-RMS / FN-RMS /
MYOD1-mutant; WD / DD / myxoid / pleomorphic / NOS), because that is the resolution the
paper argues at. There are **45**.

The atlas also carries bins that are *not* entities and are excluded from every
denominator: four **unresolved residual bins** (`Sarcoma NOS`, `RMS-NOS` and two one-off
cases — 535 samples whose diagnosis the deposit never states, a metadata-quality measure
rather than a disease) and one **control tissue** category.

RMS subtype is assigned by the rule that fusion-negative means *neither a PAX3/PAX7 fusion
nor a mutant MYOD1* — a RAS mutation is explicitly **not** the criterion, since some
FN-RMS carry none. Every call carries its evidence in `rms_call_basis`.

## Layout

```
data/          the atlas tables, T1–T25, tab-separated (see data/README.md)
  samples/     the atomic sample table — one row per GSM, 85,698 rows, gzipped
  incidence/   the curated incidence evidence behind every rate in T13
pipeline/      the harvest → classify → aggregate → publish stages, numbered in order
figures/       F18–F24 as editable PDFs (Arial, TrueType-embedded)
workbooks/     the same tables as Excel, built for clean Google Sheets import
docs/          the interactive gap map, and the working briefing note
```

## Publishing your own copy

```powershell
.\setup-github.ps1 -User <your-github-username>     # Windows
./setup-github.sh <your-github-username>            # macOS / Linux
```

Installs the GitHub CLI if needed, signs in through the browser once, creates the
repository, pushes, and enables Pages so the gap map is live.

## Reproducing it

```bash
git clone https://github.com/<you>/sarcoma-epigenomic-atlas
cd sarcoma-epigenomic-atlas
pip install openpyxl matplotlib

# everything downstream of the harvest, from the shipped tables — about a minute
cd pipeline
python3 08_aggregate_entities.py      # entity-level counts        -> T3
python3 13_join_incidence.py          # incidence join             -> T13
python3 12_merge_controlled_access.py # EGA + St Jude mapping      -> T15
python3 20_make_figures.py            # F18–F24
python3 21_build_workbook.py          # the atlas workbook
python3 23_build_gap_map.py           # the interactive gap map
```

Stages `01`–`04` re-run the GEO sweep and classification from scratch. That is several
hours and ~450,000 sample records against NCBI E-utilities; the classified output is
shipped so you do not have to. `09`–`11` re-harvest EBI and EGA (minutes).

Paths resolve relative to the checkout. Override with `SASS_ROOT`, `SASS_DATA`,
`SASS_FIGURES`, `SASS_WORK`.

## Known limits

- **Absence means absence of a public, entity-labelled deposit** — not proof no experiment
  was done. GSE140686 is the standing proof that this failure mode is real and large.
- **A sample count is not an information count.** 4,511 osteosarcoma epigenomes are not
  4,511 independent observations; `data/T3_entity_counts.tsv` carries distinct models,
  distinct PDX and distinct patients alongside the raw totals.
- **Incidence anchors are curated, not computed.** Each carries its basis in T13, and
  `US_low`/`US_high` bracket the disagreement between US registry coding and European
  central-pathology-review series, which for several entities differ two-fold.
- **Classification is rule-based over free text.** The correction logs
  (`T9b`, `T9c`) record every reassignment and why, so you can audit the calls rather than
  trust them. Corrections applied to date are listed in [`docs/methods.md`](docs/methods.md).

## Citing

Please cite the white paper when it appears. Until then, cite this repository — see
[`CITATION.cff`](CITATION.cff).

## Licence

Code (`pipeline/`) under the [MIT licence](LICENSE). Data, figures and documentation
(`data/`, `figures/`, `workbooks/`, `docs/`) under
[CC BY 4.0](LICENSE-DATA) — reuse freely, with attribution.

Underlying records remain governed by their source archives. EGA datasets listed in
`T14`/`T15` are metadata only; the data itself requires a data-access agreement with the
committee named in each row.
