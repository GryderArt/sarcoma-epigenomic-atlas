# Methods

How the atlas is built, what each decision assumes, and every correction applied — so a
reader can audit the calls rather than trust them.

## 1. Finding the data

**Two sweeps, because one is not enough.** A disease-term sweep of GEO
(`01_sweep_disease_terms.py`) misses series whose record never names the disease — a study
titled *"Chromatin landscape of RH4 cells"* is invisible to a query for
"rhabdomyosarcoma". So a second sweep queries by **model name** instead
(`02_sweep_model_names.py`, 1,069 queries built from the model catalogue), which recovered
1,644 series the disease sweep had not seen.

**A third blind spot: the bare word.** GSE140686 — the largest sarcoma epigenomic deposit
in existence — describes itself only as "sarcoma". Neither an entity query nor a model
query returns it. The `PANSARC` term group (`"sarcoma"`, `"sarcomas"`, `"soft tissue
sarcoma"`, `"bone sarcoma"`) exists solely because a coauthor's list surfaced a dataset
the pipeline had missed. Any survey of this kind should assume it has a blind spot it has
not found yet.

**Beyond GEO.** `09_harvest_ebi.py` covers ArrayExpress/BioStudies, ENA, EGA metadata and
EpiRR; `10_harvest_ega.py` enumerates the whole EGA dataset catalogue, because EGA's
apparent `query=` parameter is silently ignored and the catalogue must be pulled whole and
filtered locally. The St Jude CSTN inventory came from the portal's undocumented JSON
gateway plus the source publications.

Harvest is at **sample level**, not series level. A series is a bag of heterogeneous
samples; only `!Sample_characteristics`, `!Sample_title`, `!Sample_source_name` and
`!Sample_description` say what a given file actually is.

## 2. Classification

Rule-based over free text, with the evidence recorded per call. Five things are assigned:
assay class, epigenetic target, sample type, model identity, disease entity.

**Model matching** uses a normalised n-gram token index rather than a regex alternation
over ~1,000 names (17 ms → 0.07 ms per sample). Antibody catalogue numbers are stripped
before matching, because Novus `NB100-…` otherwise matches the CHP-100 alias `NB-100`.

**RMS subtype** follows the rule that fusion-negative means *neither a PAX3/PAX7 fusion
nor a mutant MYOD1*. A RAS mutation is explicitly not the criterion — some FN-RMS carry
none. GEO almost never writes "fusion-negative"; it writes the driver
(`RMS primary tumor, NRAS`), so within a series that systematically genotypes its cases,
the absence of a named PAX fusion is the usable signal. Precedence
(`06_rms_subtype_rule.py`):

1. mutant MYOD1 named in the sample record → **RMS-MYOD1**
2. a PAX3/PAX7 fusion named (any partner, plus `t(2;13)` / `t(1;13)`) → **FP-RMS**
3. for models, the curated identity in `T2_models.tsv` → that subtype
4. a non-PAX RMS fusion named (VGLL2, SRF::NCOA2, TEAD1::NCOA2) → **FN-RMS**
5. explicit fusion-negative or embryonal wording → **FN-RMS**
6. patient material with a named driver gene, or in a systematically-genotyped series →
   **FN-RMS** by absence
7. otherwise → **RMS-NOS**

Two guard rails matter. Only the RMS-NOS bin is re-derived, so a call already made is
never demoted. And **curated cell line identity outranks series context** — RH4, RH30,
RH41 and CW9019 are PAX-fusion lines that carry TP53 and other mutations, and rule 6 would
otherwise sweep them into FN-RMS.

**Subtype calls use sample-level text only.** Adding the series title once made
GSE140686's headline ("PAX3-FOXO1 establishes myogenic super enhancers") mark all eight of
its tumours fusion-positive, including the fusion-negative ones.

## 3. Recovering GSE140686

Koelsche et al., *Nat Commun* 2021 (PMID 33479225) deposited 1,505 methylation arrays as
the DKFZ sarcoma classifier reference and validation sets. Every GEO sample is titled
`sarcoma classifier reference case N` with characteristics reading `tissue: sarcoma`. **No
diagnosis appears anywhere in the GEO record.**

The join exists but is hidden in two places at once: GEO buries a `REFERENCE_SAMPLE N` /
`VALIDATION_SAMPLE N` identifier in `!Sample_description` (a field most harvesters never
read), and the diagnosis for that identifier lives only in the paper's Supplementary Data
1 and 3. `05_recover_GSE140686.py` fetches the supplement from Europe PMC, joins on that
identifier, and maps 1,505 samples onto 46 entities: 1,315 sarcoma, 154 benign or
non-neoplastic control, 36 non-sarcoma.

Subtype evidence is graded, because the sources differ in strength: a stated fusion beats
a methylation class, which beats histology alone. The grade is carried in
`subtype_evidence` for every recovered sample. Benign mimics and non-sarcoma controls are
retained and labelled rather than deleted — they are part of what the deposit contains.

## 4. Counting

**Regulatory** epigenomics means ChIP-seq, ChIP-exo, ChIP-chip, CUT&RUN, CUT&Tag,
ATAC-seq, scATAC-seq, DNase-seq, FAIRE-seq, MNase-seq, Hi-C, HiChIP, Micro-C, Capture-HiC,
ChIA-PET and 4C-seq, **excluding input and IgG controls** — a control is not a profile of
anything, and counting them inflates sparse entities most.

**Patient-derived** means tumour tissue, metastasis, recurrence, PDX or patient-derived
organoid, reported separately from tumour tissue alone. A synovial sarcoma organoid grown
from a patient's tumour is not a decades-old cell line, and collapsing the two would have
called synovial sarcoma a hard zero when it is not.

**Duplicates** are flagged, not deleted (`is_duplicate`, `duplicate_of_gsm`,
`series_redundant`), so a reader can choose the denominator. Headline counts use the
de-duplicated set.

**Entities versus bins.** 45 named diagnostic entities carry every denominator. Four
unresolved residual bins (`Sarcoma NOS`, `RMS-NOS` and two one-off cases, 535 samples) and
one control tissue category are excluded from all "of N" figures — they are a
metadata-quality measure, not diseases.

## 5. Incidence

Anchors are **curated, not computed**. Each carries its basis in `T13_incidence_vs_data.tsv`,
and `US_low`/`US_high` bracket a real disagreement: US registry coding and European
central-pathology-review series (NETSARC) differ roughly two-fold for several entities,
because central review reclassifies what registries code by default. Entities with no
published population rate anywhere are flagged `no_population_rate=Y` and excluded from
every per-case ratio rather than given a guessed denominator.

The 18,020 US figure is topography-coded and therefore excludes sarcomas coded to the
organ they arise in — GIST, uterine leiomyosarcoma. Extrapolating the French central-review
rate of 70.7 per million gives ~24,000/yr. The undercount runs in the adult direction, so
the burden gap is if anything larger than reported.

## 6. Verifying the zeros

Every entity called zero was subjected to an adversarial sweep designed to falsify it:
GEO at series and sample level, SRA across all library strategies, ENCODE, ChIP-Atlas
(845,824 experiments), ArrayExpress/BioStudies, GDC/TARGET, dbGaP, Cellosaurus cell-line
rosters, and PubMed plus Europe PMC full text including preprints. Caveats that survived
are recorded per entity in the `verification_caveat` column of T13. Three shape the
wording:

- **Angiosarcoma is a human zero, not an absolute one.** Canine haemangiosarcoma — the
  accepted spontaneous model — has CUT&Tag (GSE304509) and ChRO-seq (GSE150705).
- **Myxofibrosarcoma is the weakest call.** PMID 39789291 reports ATAC-seq on an MFS PDX
  but states no accession, and no matching series exists in GEO or SRA. It is a dataset
  that exists and is unfindable — reported as "no retrievable data".
- **"DNA methylation profiling", not "methylation arrays".** DFSP has one nanopore
  methylation tumour; CCSK has cfDNA RRBS.

## 7. Corrections applied

Each was a real error caught against ground truth, and each is logged.

| # | What was wrong | How it was found | Fix |
|---|---|---|---|
| 1 | GSE140686's 1,505 samples had no entity label | a coauthor's dataset list | recovered the `!Sample_description` join; `05` |
| 2 | The RMS subtype rule was lost in a pipeline rebuild — the Gryder/Yohe tumours had drifted back to RMS-NOS | checking a known result | re-applied; the five fusion-negative tumours now read FN-RMS; `T9b` |
| 3 | Driver-genotype inference swept PAX-fusion cell lines into FN-RMS | auditing which models the rule touched | curated model identity now outranks series context; `06` |
| 4 | Series titles leaked their headline subtype onto every sample | all eight tumours in one series called fusion-positive | subtype calls restricted to sample-level text |
| 5 | H3F3A/H3F3B capture pulled paediatric glioma into the bone-tumour bin | 865 K27M/G34R samples in GCTB | residue disambiguation — GCTB is G34W/L, chondroblastoma K36M; `T9c` |
| 6 | Antibody catalogue numbers matched model aliases | Novus `NB100-…` → CHP-100 | antibody fields stripped before model matching |
| 7 | Short model aliases matched patient labels | `OS9` → `OS9-1` | length and context guards |
| 8 | ChIP input controls counted as profiles | inflated sparse entities most | excluded from every regulatory count |
| 9 | Mixed genome builds and mouse bigWigs compared at human coordinates | models beat the patient-patient ceiling | build detection, liftOver, mouse dropped; then restricted to within-study comparisons after same-study ρ 0.714 vs cross-study 0.452 |
| 10 | "Entities catalogued" mixed diseases with residual bins and a control tissue | inconsistent denominators across figures | `entity_kind`; 45 named entities carry every denominator |

## 8. What this cannot tell you

**Absence in this atlas means absence of a public, entity-labelled deposit.** It is not
proof that no experiment was done. Data under controlled access, or deposited without an
entity label, is invisible to any search of this kind — and GSE140686 is the standing proof
that the second failure mode is real and large. `T14`/`T15` quantify the first: 1,207
regulatory epigenomes exist behind data-access committees.

**A sample count is not an information count.** Sample totals are inflated by treatment
arms, replicates and re-deposits of the same material. `T3_entity_counts.tsv` carries
distinct models, distinct PDX and distinct patients alongside them.

**The classifier is rules over free text**, and free text is written by humans in a hurry.
The correction table above is not a list of problems that have been solved — it is a
demonstration of the error rate, and the eleventh error has not been found yet.
