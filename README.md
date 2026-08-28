# AfroTB Phylogeny-Aware GNN — Project Handoff README

Phylogeny-aware GNN for joint prediction of TB drug resistance and
*M. tuberculosis* lineage on the Afro-TB dataset. See `CLAUDE.md` for full
project rules, team ownership, and coding conventions — this file is the
practical entry point for regenerating and using Person 1's processed data.

## 1. Where the raw data is

`data/raw/Afro_TB/` (gitignored — not committed, must exist locally):

- `0-StartHERE_Afro-TB.xlsx` — primary source. Sheet `AfroTB`: 13,753 isolates
  × 157 mutation columns (row 4 = mutation names, row 5 = `Name, Country,
  Lineage, Drug` headers, data from row 6). Each mutation cell holds either a
  placeholder (`_`, `-`, blank) or a drug code (see §6).
- `Lineage-drug-resitance-classifiation.xlsx` — source of the `Country`,
  `Lineage`, `Drug` columns (same content as in `labels.csv`).
- `WHO-resistance-associated-mutations.xlsx` — global WHO mutation→drug
  catalog (12 drug codes). Reference only; **not** the source of `y_amr`
  (see §6).
- `Undescribed-mutations.xlsx`, `Validation-strains.xlsx`,
  `Acession-Numbers.xlsx` — supporting reference sheets, not yet consumed by
  any script.
- `AFRO_TB_VCF/` — 13,753 per-isolate VCFs, one file per isolate, no
  subdirectories. Filename convention (confirmed by
  `scripts/validate_vcf_mapping.py` against every file in the directory):
  `<ID>_MT.vcf.gz`, e.g. `ERR036186_MT.vcf.gz` → canonical ID `ERR036186`.
  Verified exact 1:1 with `sample_ids.csv` (see §4.1). **Not yet parsed for
  variants — identity mapping only so far.**
- `AFRO_TB_dataset/` — per-isolate annotation `.txt` files. **Not yet
  processed by anything**, including ID-mapping validation.

## 2. How to regenerate processed data — Start Here

`data/processed/` is gitignored — a fresh clone has none of it. Regenerate
everything in this exact order from the repo root, using `.venv`. Each
script is standalone, reads only the raw source and/or already-generated
files it needs, and has its own strict validation (asserts on failure).

```
.venv/Scripts/python.exe scripts/prepare_afrotb_matrix.py   # raw xlsx -> features.csv, labels.csv, dataset_metadata.json
.venv/Scripts/python.exe scripts/create_sample_ids.py       # features.csv, labels.csv -> sample_ids.csv
.venv/Scripts/python.exe scripts/create_y_amr.py            # raw xlsx, sample_ids.csv, labels.csv -> y_amr.csv, y_amr_metadata.json
.venv/Scripts/python.exe scripts/create_splits.py           # sample_ids.csv, labels.csv -> splits.csv, splits_metadata.json
```

This order matters: `create_sample_ids.py` must run after
`prepare_afrotb_matrix.py` (it reads `features.csv`/`labels.csv`), and both
`create_y_amr.py` and `create_splits.py` must run after
`create_sample_ids.py` (they read `sample_ids.csv` as the canonical row
order). `create_y_amr.py` and `create_splits.py` are independent of each
other and can run in either order relative to one another.

All four scripts are deterministic and reproduce the currently-committed
outputs exactly (verified byte-for-byte: row count, column count, IDs,
ordering, values, and all metadata-relevant statistics). Re-run them
whenever the raw workbook changes; don't hand-edit any file in
`data/processed/`.

## 3. What each processed file means

All files live in `data/processed/` and share the row order defined by
`sample_ids.csv` (§4) unless noted otherwise.

| File | Contents |
|---|---|
| `features.csv` | `Name` + 157 binary columns — one per catalogued resistance mutation. `1` = mutation detected (any drug code), `0` = placeholder/absent. This is `X_mutations`. |
| `labels.csv` | `Name, Country, Lineage, Drug` — raw, unmodified. `Drug` is an aggregate phenotype-style category (`Sensitive/MDR/Mono/Pre-XDR/Other/Other*`), not per-drug. `Lineage` has 10 rows with value `"-"` (missing) and an unmerged `BOV_AFRI`/`BOV-AFRI` spelling split — not yet cleaned into a `y_lineage` artifact. |
| `dataset_metadata.json` | Feature-matrix metadata: mutation names, recognized presence codes, placeholder-handling assumption (§7), source row/column layout. |
| `sample_ids.csv` | The canonical ID list (§4). |
| `y_amr.csv` | `Name` + 9 binary drug columns — the multi-label AMR target (§6). This is `y_amr`. |
| `y_amr_metadata.json` | Derivation method, D=9 decision rationale, strict validation results, and a QC comparison against `labels.csv`'s `Drug` column (§8). |
| `splits.csv` | `Name, split` (`train`/`val`/`test`) — fixed partition (§5). |
| `splits_metadata.json` | Seed, proportions, stratification method, singleton handling, per-split Drug distribution. |
| `vcf_mapping_report.json` | Identity/file-mapping audit between `sample_ids.csv` and `AFRO_TB_VCF/` (§4.1) — counts, discrepancies (none found), confirmed filename convention. Not variant data. |
| `vcf_structural_qc_report.json` | Structural integrity audit of all 13,753 VCFs (§4.2) — gzip/header/record/genotype checks, sample-ID header comparison, CHROM consistency, variant-count stats. Structural only, not biological validation. |
| `vcf_content_qc_report.json` | Content/biological-sanity audit of all 13,753 VCFs (§4.3) — reference/CHROM consistency, SNP/indel representation, REF/ALT sanity, FILTER/QUAL/genotype distributions, duplicate positions, per-isolate variant burden. Nothing filtered or corrected. |
| `mutation_matrix_vcf_crosscheck_report.json` | Cross-check of the 157-mutation catalog (`features.csv`) against SnpEff-annotated VCF calls (§4.4) — per-mutation concordance rates, documented matching method, and full transparency on approximations. |

`edge_index` / `edge_weight` (the graph component) do not exist yet — that
is Person 2's work (§9).

## 4. Canonical ID / order convention

`data/processed/sample_ids.csv` is the **single source of truth for row
order**. It was derived from `features.csv`'s row order and validated to be:
13,753 IDs, all unique, and identical (name-for-name, in order) to
`labels.csv`'s order. `y_amr.csv` and `splits.csv` were both built and
written in this exact order. **Any new artifact (graph node order,
predictions, additional label files) must be indexed against
`sample_ids.csv`'s order** — join on `Name`, don't assume row-position
equality with a differently-sourced file.

### 4.1 VCF filename → canonical ID mapping (verified)

`scripts/validate_vcf_mapping.py` audits `sample_ids.csv` against every
file in `AFRO_TB_VCF/` (identity/file mapping only — no variant parsing).
Regenerate with:

```
.venv/Scripts/python.exe scripts/validate_vcf_mapping.py   # sample_ids.csv, AFRO_TB_VCF/ -> vcf_mapping_report.json
```

**Confirmed convention**: `<ID>_MT.vcf.gz`, e.g. `ERR036186_MT.vcf.gz` →
canonical ID `ERR036186`. All 13,753 files sit directly in `AFRO_TB_VCF/`
(no subdirectories) and match this pattern exactly.

**Result (see `vcf_mapping_report.json` for full detail)**: exact 1:1
mapping confirmed — 13,753 canonical IDs, 13,753 VCF files, 13,753 unique
VCF-derived IDs, 0 missing, 0 extra, 0 duplicates on either side. **Person
2 can safely use `sample_ids.csv`'s canonical order to look up each
isolate's VCF path** (join on `Name` → `<Name>_MT.vcf.gz`), consistent with
`features.csv`, `labels.csv`, `y_amr.csv`, and `splits.csv`.

### 4.2 VCF structural integrity audit

`scripts/audit_vcf_structure.py` performs a **read-only, structural-only**
audit of all 13,753 VCFs (Step 2, after Step 1's filename/ID mapping in
§4.1) — this is not biological validation and does not interpret variant
calls. Regenerate with:

```
.venv/Scripts/python.exe scripts/audit_vcf_structure.py   # AFRO_TB_VCF/ -> vcf_structural_qc_report.json
```

**What was checked, per file**: gzip validity and non-emptiness; presence
of `##` header lines and a `#CHROM` header with the required fixed columns
(`#CHROM POS ID REF ALT QUAL FILTER INFO`, plus `FORMAT`/sample columns);
every variant record's field count against the header, `POS` integer
validity, non-empty `REF`/`ALT`; genotype-field column-count consistency
against `FORMAT` (structural only — no biological interpretation); the
header's sample-column name compared against the Step 1 canonical
filename-ID; CHROM-value consistency across the whole collection; and
per-file variant-count statistics (dataset-wide min/max/mean/median plus
an explicit Tukey-fence outlier rule — files are reported, never
discarded).

**Result** (full detail in `data/processed/vcf_structural_qc_report.json`):
all 13,753 files are valid gzip, non-empty, correctly headered, and have
zero malformed variant records or malformed genotype fields. Every file
uses the same single CHROM convention (`M.tuberculosis_H37Rv`). Variant
counts range 4–197,179 per file (median 1,150); 69 files are statistically
high outliers by the Tukey-fence rule, 0 are structurally invalid, 0 have
zero variants.

The one anomaly category found: **1,604 of 13,753 files (~11.7%) have a
VCF-internal sample-column name that does not exactly equal the
filename-derived ID** — e.g. filename `ERR171145_MT.vcf.gz` but header
sample column `ERR171145/ERR171145.sorted.rmdup.bam`. In every one of
these 1,604 cases the header string **starts with** the correct ID; the
mismatch is a BAM-path-derived naming artifact from the original
variant-calling pipeline (some isolates' headers carry
`_library1.sorted/..._library1.sorted.sorted.rmdup.bam` or
`/....sorted.rmdup.bam` suffixes), not a different identity. Not corrected
here — reported only, with full file list in the JSON report.

**Person 2 can proceed with VCF parsing**, using the filename-derived ID
(already verified 1:1 against `sample_ids.csv` in §4.1) as the sample
identity — do not rely on the VCF's internal sample-column header for ID
matching, since ~11.7% of files don't carry a bare accession there.

### 4.3 VCF content / biological-sanity audit

`scripts/audit_vcf_content.py` goes beyond §4.2's structural checks to
audit **content** that could materially affect SNP-distance / phylogeny
construction — still **read-only**, still no filtering or "cleaning" of
any record. Regenerate with:

```
.venv/Scripts/python.exe scripts/audit_vcf_content.py   # AFRO_TB_VCF/ -> vcf_content_qc_report.json
```

**What was checked** (aggregated across all 13,753 files, ~20.8M variant
records): reference-assembly consistency; SNP vs. indel vs. MNP
representation and multiallelic records; REF/ALT allele-character sanity;
FILTER-value distribution; QUAL distribution; genotype-call distribution
and missingness; duplicate genomic positions; per-isolate variant burden
(total and SNP-only), with an explicit Tukey-fence outlier rule. Every
unusual record is counted and a bounded set of examples is kept for
investigation — nothing is discarded.

**Results** (full detail in `vcf_content_qc_report.json`):

- **Reference consistency**: the `##reference=` header path differs
  across 3 values, but this is purely a processing-batch directory
  difference (`Algeria/...`, `MTBseq_source-master/...`,
  `TB_africa/...`) — all 3 point to the **same reference FASTA filename**
  (`M._tuberculosis_H37Rv_2015-11-13.fasta`), and every file shares one
  `##contig=` identifier/length and one `CHROM` value
  (`M.tuberculosis_H37Rv`). Reference assembly is consistent.
- **Variant representation**: 19,164,310 SNPs, 1,670,709 indels; 4,527
  multiallelic records (ALT with >1 allele in one line).
- **REF/ALT sanity**: 0 truly invalid (non-ACGTN) alleles anywhere in the
  dataset. 108,314 REF and 109,140 ALT values use **lowercase** base
  letters instead of uppercase — a formatting-convention difference, not
  a validity problem (all are valid nucleotides case-insensitively).
- **FILTER distribution**: 100% of records carry `.` (no `PASS`/failed
  distinction was set by the calling pipeline) — worth knowing before
  assuming FILTER can be used to subset variants.
- **QUAL distribution**: 0 missing/unparseable; mean 203.6 (histogram
  range ~3–225).
- **Genotype distribution**: haploid calls (ploidy 1, consistent with the
  `bcftools call --ploidy 1` command recorded in each VCF header) —
  observed tokens are `1` (19,126,773) and missing `./.` (1,703,454); 0
  genotype tokens in an unexpected shape. **Missingness rate: 8.18%** of
  genotype calls.
- **Duplicate positions**: 1,070 positions repeated within a file with
  *differing* REF/ALT (legitimate split multiallelic representation); **0**
  exact duplicate (same CHROM+POS+REF+ALT) records.
- **Per-isolate variant burden**: total variants min 4 / median 1,150 /
  max 197,179; SNP-only min 4 / median 1,050 / max 196,960. High-end
  outliers exist (consistent with §4.2's 69 flagged files) — retained, not
  removed.

Nothing here blocks Person 2 from parsing variants; the FILTER-field and
lowercase-REF/ALT observations are worth being aware of when writing a
VCF parser (e.g. don't filter on `FILTER=="PASS"` — nothing would survive;
uppercase alleles before comparing sequences).

### 4.4 Mutation-matrix ↔ VCF cross-check

`scripts/crosscheck_mutations_vcf.py` checks whether the 157 resistance
mutations encoded in `features.csv` (`X_mutations`) are independently
observable in the SnpEff-annotated VCFs in
`data/raw/Afro_TB/AFRO_TB_ANNOTATION_VCF/` — a consistency check between
`X_mutations`, `y_amr`, and the VCF-derived genomic data Person 2 will use
for phylogeny. Read-only; no VCF, `features.csv`, or `y_amr.csv` change.
Regenerate with:

```
.venv/Scripts/python.exe scripts/crosscheck_mutations_vcf.py   # AFRO_TB_ANNOTATION_VCF/, features.csv, sample_ids.csv -> mutation_matrix_vcf_crosscheck_report.json
```

**Matching method** (stated explicitly so the result isn't overclaimed):
149 of 157 catalog entries (substitutions, synonymous changes,
single-residue del/dup, range del/ins) are matched by **exact string
equality** against the annotation's `HGVS.p` field, after normalizing the
AfroTB name into the same `p.<Xxx><pos><...>` shape SnpEff uses (1-letter
codes converted via the standard universal amino-acid table — not
dataset-specific information). The remaining 8 are frameshift (`...fs`)
entries, matched only by **gene + amino-acid position** (SnpEff's
frameshift notation carries extra detail — e.g. a terminal-stop offset —
the catalog name doesn't encode), so these are reported as an
**approximate, position-based match**, not exact. 0 of the 157 entries
were unparseable.

**Result** (full per-mutation detail in
`mutation_matrix_vcf_crosscheck_report.json`, aggregated across all
13,753 isolates — not a per-isolate dump):

- **Mean concordance rate across the 150 mutations with ≥1 positive
  isolate: 99.81%.**
- **145 of 150 mutations have 100% concordance** (every isolate marked
  positive in `features.csv` has a matching variant call in its
  annotated VCF).
- **0 mutations have 0% concordance.**
- 5 mutations have partial concordance (90–99.8%): `pncA V180F` (10/11),
  `katG S315N` (68/74), `rpsL K88M` (15/16), `embB M306L` (37/39),
  `rpsL K43R` (1031/1033) — small counts of isolates marked positive in
  the XLSX without a matching VCF-annotated call. Not investigated
  further or corrected here; flagged for whoever owns label QA.
- 15 mutations have a small number of isolates (1–7 each) marked
  *negative* in `features.csv` where the VCF annotation *does* show the
  variant (e.g. `katG S315T`: 7 isolates) — again small in absolute
  terms and reported, not corrected.
- 7 catalog mutations have **zero** positive isolates in `features.csv`
  dataset-wide (e.g. `pncA Val130Val`, `rpoB Gln432His`) — nothing to
  cross-check for these; consistent with them simply being rare/absent in
  this cohort, not a VCF-side problem.

**Bottom line**: `X_mutations`, `y_amr`, and the VCF-derived genomic data
are demonstrably compatible — 99.81% mean concordance with 0 mutations at
zero-concordance is strong evidence the three sources describe the same
underlying calls, safely citable as "we audited the complete VCF
collection, cross-checked it against the mutation matrix, and confirmed
compatibility" rather than needing any correction before Person 2 starts.

## 5. Train / validation / test split

`data/processed/splits.csv` — 70% / 15% / 15% (train=9,629 / val=2,062 /
test=2,062), stratified on `labels.csv`'s `Drug` column, seed `42`.

- Every ID appears in exactly one split; order matches `sample_ids.csv`.
- The single `Drug="Other"` isolate (`SRR1577832`) can't be stratified
  across 3 splits and is deterministically placed in `train` (documented in
  `splits_metadata.json`).
- This is a **per-sample** split only — it does not yet account for
  phylogenetic adjacency/leakage across the graph Person 2 will build (see
  §9's leakage note).
- Full per-split Drug distribution is in `splits_metadata.json`.

## 6. y_amr = 9 drugs

`y_amr.csv` has exactly **9** drug columns: `RIF, INH, EMB, PZA, STM, LEV,
CAP, ETH, LZD`. These are re-derived directly from the raw workbook by
preserving the drug code in each of the 157 mutation cells (rather than the
flat 0/1 in `features.csv`, which discards which drug each mutation
confers resistance to): `y_amr[isolate, drug] = 1` if any mutation cell for
that isolate carried that drug's code, else `0`.

**D=9, not 12** — deliberately. The WHO catalog (`WHO-resistance-
associated-mutations.xlsx`) lists 12 drug codes (adds `AMI, KAN, MXF`), but
those 3 never appear as a value in this dataset's mutation cells, so there
is no ground truth to populate them from — adding them as all-zero columns
would fabricate "tested and negative" where the data actually says
"untested." Full rationale in `y_amr_metadata.json`.

## 7. Placeholder assumption

In the raw workbook's 157 mutation columns, three placeholder values
(`_`, `-`, blank/`None`) all appear. The source publication (Laamarti et
al., *Scientific Data*, 2023) and the workbook do not define whether these
distinguish "not detected" from "not tested" or "not applicable." Both
`features.csv` and `y_amr.csv` treat **all three as "mutation absent" (0)**.
This is a provisional simplification — revisit if authoritative
documentation (e.g. UM6P Afro-TB database) becomes available. Recorded in
`dataset_metadata.json.preprocessing_assumption`.

## 8. Known QC issues

- **`Lineage` column**: 10 isolates have `Lineage = "-"` (missing); label
  spelling is split across `BOV_AFRI` and `BOV-AFRI` (likely the same
  class). Not yet cleaned — do not assume `labels.csv`'s `Lineage` is
  ready to use as `y_lineage` without addressing this first.
- **`Drug` vs. `y_amr` inconsistencies** (from `y_amr_metadata.json`'s QC
  pass, non-destructive — nothing was corrected):
  - 1 isolate labeled `Sensitive` has a non-zero `y_amr` row.
  - 904 isolates labeled `Mono` have `y_amr` drug-counts ≠ 1 (up to 5
    drugs flagged) — `Drug`'s categories do not map to a literal count of
    `y_amr` positives; do not assume they're reconcilable without further
    investigation.
  - 0 isolates in a resistant category (`MDR/Mono/Pre-XDR/Other/Other*`)
    have an all-zero `y_amr` row (consistent).
  - Full per-category `y_amr` positive-count stats and the exact
    inconsistent IDs are in `y_amr_metadata.json`.
- **Sample-ID namespace vs. VCFs**: RESOLVED — `sample_ids.csv` was
  cross-validated against `AFRO_TB_VCF/` filenames and confirmed exact 1:1
  (§4.1, `vcf_mapping_report.json`). Not yet checked against
  `AFRO_TB_dataset/*.txt`.
- **VCF sample-column header naming** (§4.2): 1,604/13,753 files (~11.7%)
  have a VCF-internal sample-column name that is a BAM-path-derived
  string rather than a bare accession (though it always starts with the
  correct ID). Use the filename, not the header, for sample identity.
- **VCF FILTER field is uninformative**: 100% of ~20.8M variant records
  carry `FILTER="."` (§4.3) — the calling pipeline never set `PASS` vs. a
  failure reason, so `FILTER` cannot be used to subset variants.
- **REF/ALT case inconsistency**: 108,314 REF and 109,140 ALT values use
  lowercase base letters instead of uppercase (§4.3) — all are valid
  nucleotides, but any VCF parser should uppercase alleles before string
  comparison rather than assume a single case convention.
- **Genotype missingness**: 8.18% of genotype calls are missing (`./.`)
  across the dataset (§4.3) — expected for sequencing data, but worth
  accounting for in any per-position or per-isolate completeness
  assumption.
- **`features.csv` vs. VCF-annotation concordance** (§4.4, non-destructive
  QC, nothing corrected): 99.81% mean concordance overall; 5 mutations
  have 90–99.8% concordance (small numbers of XLSX-positive isolates
  without a matching VCF call) and 15 mutations have 1–7 isolates marked
  negative in the XLSX where the VCF annotation shows the variant. Full
  per-mutation detail in `mutation_matrix_vcf_crosscheck_report.json`.

## 9. What Person 2 owns

Per `CLAUDE.md`: Graph Neural Network architecture, multi-task learning,
model training, validation, model evaluation. Concretely, next steps on top
of what Person 1 has delivered:

- Build the SNP-distance / phylogenetic representation from
  `AFRO_TB_VCF/*.vcf.gz`, keyed to `sample_ids.csv`'s canonical IDs. The
  filename→ID mapping is already verified exact 1:1 (§4.1) — safe to use
  directly, no further ID reconciliation needed for the VCF directory. The
  VCF collection has also passed a structural (§4.2) and content/
  biological-sanity (§4.3) audit, plus a mutation-matrix cross-check
  (§4.4) — safe to parse directly, with the caveats noted in §8 (FILTER
  is uninformative, REF/ALT case varies, ~8% genotype missingness).
- Produce `edge_index` (2×E) and `edge_weight` (E) from that representation.
- Consume `features.csv` as node features (`X_mutations`), `y_amr.csv` as
  the multi-label AMR target, and (once cleaned) a `y_lineage` target from
  `labels.csv`'s `Lineage` column.
- Use `splits.csv` for train/val/test — do not re-split.
- Consider whether the per-sample split in `splits.csv` needs
  graph-aware adjustment (e.g. isolates connected by a phylogenetic edge
  landing in different splits) before training — this was explicitly out of
  scope for Person 1's split design.

## 10. What Person 3 owns

Per `CLAUDE.md`: Random Forest baseline, XGBoost baseline, global/non-
African comparison where data permits, evaluation metrics, explainability,
result visualization. Concretely:

- Baselines should train on `features.csv` (`X_mutations`) against
  `y_amr.csv` (multi-label) and/or `labels.csv`'s `Drug` column (categorical
  baseline task, kept intact per §6/§8) using the same `splits.csv`.
- Evaluation metrics should be computed per-drug (from `y_amr`) and
  optionally against the aggregate `Drug` category, given the two are not
  fully consistent (§8) — report both rather than silently picking one.
- Explainability/visualization work should trace feature importance back to
  the 157 mutation names in `dataset_metadata.json.feature_names`.
