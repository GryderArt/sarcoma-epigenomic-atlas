# Methods

How the atlas is built, what each decision assumes, and every correction applied — so a
reader can audit the calls rather than trust them.

## 1. Finding the data

**Two sweeps, because one is not enough.** A disease-term sweep of GEO
(`01_sweep_disease_terms.py`) misses series whose record never names the disease — a study
titled *"Chromatin landscape of RH4 cells"* is invisible to a query for
"rhabdomyosarcoma". So a second sweep queries by **model name** instead
(`02_sweep_model_names.py`, 1,069 queries built from the model catalog), which recovered
1,644 series the disease sweep had not seen.

**A third blind spot: the bare word.** GSE140686 — the largest sarcoma epigenomic deposit
in existence — describes itself only as "sarcoma". Neither an entity query nor a model
query returns it. The `PANSARC` term group (`"sarcoma"`, `"sarcomas"`, `"soft tissue
sarcoma"`, `"bone sarcoma"`) exists solely because a coauthor's list surfaced a dataset
the pipeline had missed. Any survey of this kind should assume it has a blind spot it has
not found yet.

**Beyond GEO.** `09_harvest_ebi.py` covers ArrayExpress/BioStudies, ENA, EGA metadata and
EpiRR; `10_harvest_ega.py` enumerates the whole EGA dataset catalog, because EGA's
apparent `query=` parameter is silently ignored and the catalog must be pulled whole and
filtered locally. The St Jude CSTN inventory came from the portal's undocumented JSON
gateway plus the source publications, and `15_harvest_stjude_viz.py` resolves its two
epigenetic browsers to sample x assay grain (§3b).

Harvest is at **sample level**, not series level. A series is a bag of heterogeneous
samples; only `!Sample_characteristics`, `!Sample_title`, `!Sample_source_name` and
`!Sample_description` say what a given file actually is.

## 2. Classification

Rule-based over free text, with the evidence recorded per call. Five things are assigned:
assay class, epigenetic target, sample type, model identity, disease entity.

**Model matching** uses a normalized n-gram token index rather than a regex alternation
over ~1,000 names (17 ms → 0.07 ms per sample). Antibody catalog numbers are stripped
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
its tumors fusion-positive, including the fusion-negative ones.

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
retained and labeled rather than deleted — they are part of what the deposit contains.

## 3b. The St Jude CSTN browsers, and who says what a sample is

The two CSTN "epigenetic landscape" pages on `viz.stjude.cloud` are ProteinPaint embeds
that ship their whole track manifest inline, so the browser contents can be enumerated
exactly — 398 tracks on 28 models — at a resolution EGA's own description ("ChIP-Seq files
for RMS, 242 samples") does not reach. `15_harvest_stjude_viz.py` does that.

**The entity calls are curated, not inferred.** The first version of that stage read the
subtype off the label suffix — `(ERMS)`, `(ARMS)`, `(SCLEROS)`. The St Jude data
administrator then supplied the definitive per-sample list, which is transcribed verbatim
in the `ADMIN_EMAIL` block of the stage and is the primary source for every call. The
CSTN portal's own model table (`T17`, harvested independently) is read as a second source
and supplies the fusion partner, which the administrator's list does not give: two models
are PAX3::FOXO1 and two are PAX7::FOXO1. Every row of `T28` carries `entity_call_basis`,
`subtype_evidence`, `fusion_or_driver` and `fusion_source`, so no call in that table has to
be taken on trust.

Applying the administrator's list confirmed eleven fusion-negative and three
fusion-positive models against the label suffix and corrected two calls the suffix had got
wrong:

- `SJHGS015726_X2` is **"high grade sarcoma"** to both the administrator and the portal.
  The label-based rule had mapped `SJHGS` to UPS/MFH. Undifferentiated pleomorphic sarcoma
  is a specific diagnosis; "high grade sarcoma" is the absence of one, so it now sits in
  the `Sarcoma NOS` residual bin and is out of every "of 45 entities" denominator.
- `SJRHB015720_X1` is **"sclerosal RMS"** — a histology. The suffix rule had asserted
  `RMS-MYOD1`, which is the WHO entity's defining genotype. Neither the administrator, the
  portal, nor the literature states MYOD1 status for this model; the portal reports PIK3CA,
  which co-occurs in roughly a third of MYOD1 L122R cases (PMID 24793135) and is
  suggestive, not confirmatory. The atlas bin is kept but `subtype_evidence` reads
  `histology_spindle_sclerosing`, not a genotype.

**Three discrepancies are reported and left standing**, printed by the stage on every run:

1. `SJRHB010463_X16` is on the browser page and in the portal (PAX3::FOXO1, with ChIP-seq
   deposited as EGAD00001003432) but does **not** appear in the administrator's list. It is
   kept, with its entity taken from the browser label and flagged as such.
2. `SJOS010930_X1` has been **withdrawn from CSTN and is no longer available** — a fact
   available from nowhere but the administrator. Its 16 tracks are still drawn on the
   browser page. They are kept, marked `withdrawn`, and excluded from the enumerated
   resource count in `T27`. The withdrawal is independently corroborated: the live CSTN
   portal API does not return this model either.
3. `SJRHB013757_X1` is on the browser page; the portal carries only `SJRHB013757_X2`. A
   fusion is a property of the patient's tumor rather than of the passage, so the sibling
   row supplies PAX7::FOXO1 — recorded as `fusion_source = sibling passage`, not as an
   exact match.

One near miss is worth recording because it would trip a naive parser: the portal lists
`PAX3` in the mutation field for `SJRHB012405_X1`, which both the administrator and the
portal's own subtype field call **fusion negative**. A *PAX3 point mutation* is not a
*PAX3::FOXO1 fusion*, and the RMS rule in §2 turns on the fusion.

None of this material is added to any sample total. It is already counted once inside the
EGA bar of F25; what the enumeration adds is the assay-level resolution EGA does not
publish, shown in **F26**.

## 4. Counting

**Regulatory** epigenomics means ChIP-seq, ChIP-exo, ChIP-chip, CUT&RUN, CUT&Tag,
ATAC-seq, scATAC-seq, DNase-seq, FAIRE-seq, MNase-seq, Hi-C, HiChIP, Micro-C, Capture-HiC,
ChIA-PET and 4C-seq, **excluding input and IgG controls** — a control is not a profile of
anything, and counting them inflates sparse entities most.

**Patient-derived** means tumor tissue, metastasis, recurrence, PDX or patient-derived
organoid, reported separately from tumor tissue alone. A synovial sarcoma organoid grown
from a patient's tumor is not a decades-old cell line, and collapsing the two would have
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

- **Angiosarcoma is a human zero, not an absolute one.** Canine hemangiosarcoma — the
  accepted spontaneous model — has CUT&Tag (GSE304509) and ChRO-seq (GSE150705).
- **Myxofibrosarcoma is the weakest call.** PMID 39789291 reports ATAC-seq on an MFS PDX
  but states no accession, and no matching series exists in GEO or SRA. It is a dataset
  that exists and is unfindable — reported as "no retrievable data".
- **"DNA methylation profiling", not "methylation arrays".** DFSP has one nanopore
  methylation tumor; CCSK has cfDNA RRBS.

## 7. Corrections applied

Each was a real error caught against ground truth, and each is logged.

| # | What was wrong | How it was found | Fix |
|---|---|---|---|
| 1 | GSE140686's 1,505 samples had no entity label | a coauthor's dataset list | recovered the `!Sample_description` join; `05` |
| 2 | The RMS subtype rule was lost in a pipeline rebuild — the Gryder/Yohe tumors had drifted back to RMS-NOS | checking a known result | re-applied; the five fusion-negative tumors now read FN-RMS; `T9b` |
| 3 | Driver-genotype inference swept PAX-fusion cell lines into FN-RMS | auditing which models the rule touched | curated model identity now outranks series context; `06` |
| 4 | Series titles leaked their headline subtype onto every sample | all eight tumors in one series called fusion-positive | subtype calls restricted to sample-level text |
| 5 | H3F3A/H3F3B capture pulled pediatric glioma into the bone-tumor bin | 865 K27M/G34R samples in GCTB | residue disambiguation — GCTB is G34W/L, chondroblastoma K36M; `T9c` |
| 6 | Antibody catalog numbers matched model aliases | Novus `NB100-…` → CHP-100 | antibody fields stripped before model matching |
| 7 | Short model aliases matched patient labels | `OS9` → `OS9-1` | length and context guards |
| 8 | ChIP input controls counted as profiles | inflated sparse entities most | excluded from every regulatory count |
| 9 | Mixed genome builds and mouse bigWigs compared at human coordinates | models beat the patient-patient ceiling | build detection, liftOver, mouse dropped; then restricted to within-study comparisons after same-study ρ 0.714 vs cross-study 0.452 |
| 10 | "Entities cataloged" mixed diseases with residual bins and a control tissue | inconsistent denominators across figures | `entity_kind`; 45 named entities carry every denominator |
| 18 | The St Jude data administrator's per-model diagnoses governed the browser tracks they were obtained for, but not the GEO samples generated on the same models — which carry no diagnosis of their own, so their subtype came from the curated catalogue instead | the same student asked why GSE174376's single-cell ATAC-seq was not under FN-RMS | 250 samples across 6 series were filed under the wrong RMS subtype, 186 of them belonging to FN-RMS. Precedence follows correction #3 one level out: the institution that derived, holds and diagnosed the model outranks our catalogue's reading of it. Matching is on the model base, asserted safe because no base carries conflicting calls across passages. Five epigenomic samples move; the rest are transcriptomic and fix the tables rather than the counts. `37`, `T32` |
| 18b | Fifteen `*_barcode` samples in GSE174376 — `protocol: Barcode dialout PCR`, `library_strategy: OTHER` — were classified scATAC-seq because the series is a single-cell ATAC study | counting GSE174376's scATAC and getting 22 where the deposit describes 7 | a barcode dial-out is a PCR readout of lineage labels, not a chromatin profile. The series-context leak of corrections #4 and #15 once more, in a record whose own words were clear enough to prevent it. Demoted to non-epigenomic; single-cell epigenomes atlas-wide fall 246 to 231. `37` |
| 17 | Controlled-access cohorts deposited above the subtype reached no entity page at all — the St Jude RMS ChIP-seq, 400 samples in EGAD00001004312 and EGAD00001006398, was filed under the bin `RMS (any subtype)` and therefore absent from FN-RMS, FP-RMS and RMS-MYOD1 alike | a PhD student looked for the St Jude PDX ChIP-seq under FN-RMS and could not find it | refusing to guess a subtype was right; refusing to say anything was not. 413 regulatory and 467 methylation samples sat in cohort-level bins, 48% of the controlled regulatory record. Stage `36` attaches each cohort to the entities it bears on at two evidence tiers — *verified*, where the St Jude data administrator's per-model list (correction #11) names them, and *parent*, where the deposit states only the parent diagnosis. Neither tier is added to any entity's totals, because the same samples cannot be counted once per subtype; they are listed with their size and their basis. No headline number moves. `36`, `T31` |
| 16 | A delta harvest adds samples that the cross-series duplicate pass never sees, so new rows carry no `is_duplicate` flag | re-harvesting to 2026-09-01 and asking what stage `32` could not reproduce | the rule that produced the existing 197 flags could not be reconstructed from the table it wrote, and a guessed replacement over-flagged by an order of magnitude when tested against those 197. Rather than substitute a different rule, new rows are left unflagged and the bound is stated: 197 of 85,698 pre-delta rows carry the flag (0.23%), so the de-duplicated denominators move by less than their rounding. `32` |
| 15 | RNA-seq samples in multi-assay series inherited the series' epigenomic assay class | asking what LGFMS/SEF's 20 "regulatory" samples actually were — all 20 were a HUVEC fusion-expression experiment, 12 of them titled `RNA_seq_*` | `library_strategy` is now authoritative over any assay word in the text; 723 samples reassigned, 688 of which had no epigenomic assay word in their own record. Regulatory total falls 6.2%, 10,556 → 9,898. Every structural claim is unchanged: still 13 entities with none, 20 with patient-derived, 11 unvalidatable. Stage `30`, log in `T29`; `04` fixed for future runs |
| 14 | `distinct_models` was read as "models that exist"; it counts models name-matched to the curated catalog, so three entities with unmatched cell-line data were miscounted | cross-checking the model bullet against `sample_type` | claim now computed from `sample_type` directly — 14 entities have never had a cell line, PDX or organoid profiled, a different set from `distinct_models = 0`; §7b |
| 13 | Normal and reference tissue counted toward an entity's regulatory total | LGFMS / SEF appeared to have 20 regulatory epigenomes while dropping out of the cell-line-vs-tissue figure | 3% of regulatory samples atlas-wide; decisive only for LGFMS / SEF (100% control tissue) and endometrial stromal sarcoma (75%). Reported rather than silently re-binned |
| 12 | Two scatter labels sat on top of each other in WP1-D and F19, so `FP-RMS` appeared to name the `Rhabdoid tumor / ATRT` point | a reader asked which dot each label meant, and whether the two were pooled | they were never pooled — separate rows in every table. Placement rewritten in `_labels.py`: measured boxes, an ownership test, leader lines; `28`, `20` |
| 11 | St Jude CSTN entity calls were read off the browser label suffix — `SJHGS` became UPS/MFH, `(SCLEROS)` asserted a MYOD1 genotype | the data administrator supplied the definitive per-sample list | curated calls keyed on model ID, with `entity_call_basis` per row; three discrepancies reported rather than resolved; `15`, §3b |

## 7b. Figures for the white-paper section

`28_whitepaper_figures.py` builds the four figures used in the "Critical assessment of
existing datasets" section: WP1 accumulation and modality mix, WP2 burden versus data,
WP3 what the atlas is made of, WP4 whether the models can be validated. Figure numbers
follow the order they are cited in the section, and the section opens on growth rather
than on deficit — the census documents a record that has expanded 87% in a decade, and a
reading that reports only the gaps misrepresents it. An earlier pair on access tiering (where the data
lives; the St Jude CSTN panel) was cut when the section moved to mentioning controlled
access only in passing; the underlying tables T15, T27 and T28 are unchanged. They are composed from T3/T4/T13/T15/T27/T28 at
render time rather than assembled from the exploratory figures, so every number on every
panel is recomputed from the tables on each run and cannot drift from them. Four panels
maximum per figure; vector PDF, Arial, `pdf.fonttype 42` so the text stays editable.

**Scatter labels carry leader lines, and placement is measured rather than guessed.**
`_labels.py` holds the shared placer used by WP1-D and F19. The first version tested a
guessed 52x11 pixel box; `Rhabdoid tumor / ATRT` (95 US cases per year, 1,052 regulatory
samples) and `FP-RMS` (110, 640) both took the slot above themselves, their guessed boxes
cleared each other by a pixel, and the pair read as one two-line label over the upper
point. Four rules now prevent that: candidate boxes are drawn and measured with the real
renderer (noting that `Annotation.get_window_extent` includes the leader, so the search
uses a bare label); a candidate is rejected if it sits nearer to a foreign marker than to
its own; candidates are tried in order of how well they point away from the nearest
neighbour; and every label gets a leader line, so the reading never rests on proximity.
Placement is done only after the scales, limits and legend are final, because it is
measured in display pixels. A label that cannot be placed cleanly is reported on stdout
rather than drawn silently.

Two claims in that section needed checking against the tables rather than against the
exploratory figures, and one of them was wrong on the first pass:

- `distinct_models = 0` means **no model has ever contributed an epigenomic sample** for
  that entity (14 of 45). It does *not* mean no model exists — 13 of those 14 do have
  epigenomic data, from patient material. The first draft read it the second way.
- Model-to-model H3K27ac similarity exceeds patient-to-patient similarity in **six of
  seven** study-internal comparisons (`T11b`), not all seven. Clear cell sarcoma by
  signal Spearman is the exception.

## 7bis. Assay class and the series-context leak

The assay detector searched each sample's own text first and fell back to series-level
context. For a multi-assay series — "ChIP-seq and RNA-seq of X", or a 10x multiome whose
ATAC and GEX halves sit side by side — that fallback stamps the series' epigenomic assay
onto its RNA samples. It is the assay-level form of the series-title leak already logged
as correction 4, and it was found by asking what a single entity's twenty "regulatory"
samples actually were.

`library_strategy` is a structured field set by the depositor and is now authoritative: a
sample GEO calls RNA-Seq is not ChIP-seq, whatever its series is titled. 723 samples were
reassigned, 688 of them with no epigenomic assay word anywhere in their own record; the
remaining 35 are almost all the GEX half of a multiome. `30_fix_rna_assay_leak.py` repairs
the shipped table, which cannot be regenerated without re-running the multi-hour GEO
sweep, and logs every reassignment to `T29_rna_assay_corrections.tsv`. Stage `04` carries
the same rule for any future full run.

The correction is one of magnitude, not of structure. The regulatory total falls 6.2%
(10,556 → 9,898) and three entities move materially — ASPS −29%, FN-RMS −20%, LGFMS/SEF
20 → 8. No entity crosses zero: still 13 with no regulatory data, 20 with patient-derived,
31 in the validation figure, 11 of them unvalidatable.

## 7c. Spelling

The atlas reports US incidence for a US-authored working group, so its own text — entity
labels, derived column values, figure captions, documentation — is US English, applied by
`29_normalize_spelling.py` using the rules in `_enus.py`.

Two categories are deliberately left as deposited. The **verbatim contents of third-party
records** keep their original spelling: a GEO `source_name` reading "Malignant peripheral
nerve sheath tumor" and a DKFZ classifier label reading "Malignant rhabdoid tumor" are
quoted, not corrected, and the gap map displays them as the depositor wrote them. So do
the **regular expressions that match those records** — European depositors write "tumor"
and "hemangioendothelioma", and Americanising the patterns would stop them matching. The
separation is enforced mechanically: TSV rewriting is scoped by column against a
deny-list of verbatim fields, and Python rewriting masks every raw-string literal before
substituting.

The rename touched 53,824 table cells and every count was verified unchanged afterwards:
re-running stages 08, 13 and 12 reproduced all 45 entity rows identically.

## 7d. Re-harvest and the delta

The full sweep (`01`–`04`) queries roughly 450,000 records and takes hours. Bringing the
atlas forward does not need it: every GEO series carries a release date, so the same term
sweep bounded to dates after the latest series already held reaches only what is new.
`31_reharvest_delta.py` runs that bounded sweep and classifies the result by **importing**
`classify` from `04`, not by reimplementing it, so a sample deposited last month is
annotated by exactly the rules that annotated one deposited in 2012. `32_merge_delta.py`
appends the delta, then `06`, `07` and `30` are re-run over the whole merged table.

The 2026-09-01 re-harvest swept 658 queries, fetched 98,638 sample records from 4,357
series the atlas did not hold, and kept 6,440 in scope. Stage `07` then moved 296 of those
out of the bone-tumor bin — 276 non-human H3.3 model systems and 20 H3.3 glioma-residue
rows — which is correction #5 operating on new material rather than a new error. Stage
`30` found nothing to correct, because the `library_strategy` guard added to `04` after
correction #15 now rejects the RNA-assay leak at classification time rather than
afterwards. Net: +681 epigenomic samples, 15,671 → 16,352.

Every structural claim in the white-paper section survived the re-harvest unchanged:
entity coverage (8/21/44 of 45 by 2010/2015/2020), 13 entities with no
regulatory epigenomics, 25 without patient-derived material, 11 of 31 not
validatable, and a 23x per-case gap between pediatric and adult disease. What moved
were magnitudes, which is what a re-harvest should move.

## 7e. Numbers in the prose

`33_whitepaper_facts.py` computes every quantity the white-paper section asserts and
writes `docs/whitepaper_facts.json`; the document build reads that file and interpolates.
Nothing numeric is typed into the document. `33` imports `28_whitepaper_figures.py` and
reuses its loaders, its `REG` / `MODERN` / `MODALITY` definitions and its bootstrap, so a
sentence and the panel it describes cannot print different numbers — they did, for six
quantities, when `33` first carried its own copies of those definitions.

## 8. What this cannot tell you

**Absence in this atlas means absence of a public, entity-labeled deposit.** It is not
proof that no experiment was done. Data under controlled access, or deposited without an
entity label, is invisible to any search of this kind — and GSE140686 is the standing proof
that the second failure mode is real and large. `T14`/`T15` quantify the first: 1,207
regulatory epigenomes exist behind data-access committees.

**A sample count is not an information count.** Sample totals are inflated by treatment
arms, replicates and re-deposits of the same material. `T3_entity_counts.tsv` carries
distinct models, distinct PDX and distinct patients alongside them.

**The classifier is rules over free text**, and free text is written by humans in a hurry.
The correction table above is not a list of problems that have been solved — it is a
demonstration of the error rate, and the sixteenth error has not been found yet.
