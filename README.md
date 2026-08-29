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
| `snp_matrix.npz` | Core-genome SNP matrix (§12), sparse `(13753, 465257)` binary, row order = `sample_ids.csv`. Confident ALT calls only; see §12 for the missingness/reference-collapse caveat. |
| `snp_matrix_sites.csv` | `pos, alt, global_call_count, missing_rate, kept_in_core` for every variant site observed dataset-wide (§12) — the column index for `snp_matrix.npz`. |
| `pyg_graph.pt` | `edge_index` (2×E) / `edge_weight` (E) as a `torch_geometric.data.Data` object, plus `sample_ids` for node-order alignment (§12). This is the graph component of the data contract. |
| `graph_construction_report.json` | Full graph-build metadata (§12): site counts, distance metric, k/mode, degree distribution, connected-component sizes, cross-split near-duplicate diagnostic, timings, and the correctness-check results. |
| `results/baseline_metrics.json` | Phase 1 tabular-baseline results (§11): RF / XGBoost / multi-task MLP, default and per-drug-tuned thresholds, val + test metrics. |

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
model training, validation, model evaluation. Tabular baselines (§11) and
graph construction (§12) have been completed cross-role in this session;
what remains is the GNN itself:

- **Done**: `edge_index` (2×E) / `edge_weight` (E) exist in
  `data/processed/pyg_graph.pt` (§12) — built from `AFRO_TB_VCF/*.vcf.gz`
  keyed to `sample_ids.csv`'s canonical order. Read §12 before using it,
  especially the missingness/reference-collapse limitation and the 98
  cross-split near-duplicate pairs it flags (not removed).
- **Remaining**: design and train the actual GNN. Consume `features.csv` as
  node features (`X_mutations`), `y_amr.csv` as the multi-label AMR target,
  and (once cleaned) a `y_lineage` target from `labels.csv`'s `Lineage`
  column (§8's `Lineage` caveat still applies — not yet cleaned).
- Use `splits.csv` for train/val/test — do not re-split. §12's cross-split
  near-duplicate diagnostic (98 val/test isolates with a near-identical
  train neighbor) is worth reading before treating test performance as a
  fully independent estimate in a transductive setting.
- `results/baseline_metrics.json` (§11) is the number to beat.

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

## 11. Phase 1 — Tabular baselines (completed)

`src/baselines/` — Random Forest, XGBoost, and a hard-parameter-sharing
multi-task MLP (`sklearn.neural_network.MLPClassifier`, shared hidden
layers + one sigmoid output per drug — no `torch` needed for this baseline),
all trained on `features.csv` → `y_amr.csv` using `splits.csv` exactly
(train only; `val` used for the MLP's early stopping and for per-drug
threshold tuning, never for model selection beyond that; `test` scored once).
Run with:

```
.venv/Scripts/python.exe -m src.baselines.run_all
```

Per-drug decision thresholds are tuned on `val` only (F1-maximizing via
`precision_recall_curve`) and applied unchanged to `test` — see
`src/baselines/metrics.py:tune_per_drug_thresholds`. Both the default-0.5
and tuned-threshold results are kept in `results/baseline_metrics.json` for
comparison. Test macro-F1: **RF 0.826 → 0.869 tuned**, **XGBoost 0.656 →
0.657 tuned** (its per-drug models have no real signal for the rarest drugs
— CAP/LZD ROC-AUC ≈ 0.5–0.71 — so no threshold fixes that), **MLP 0.762 →
0.850 tuned**. `ETH` has 0 test positives, so its AUC is `null`/excluded
from macro-AUC by design, not a bug (`metrics.py`'s documented convention).

## 12. Phase 2 — Genomic distance & graph construction (completed)

`src/phylogeny/` builds the SNP matrix, pairwise distance, and k-NN graph
from `AFRO_TB_VCF/*.vcf.gz`, reusing the filename↔ID mapping already
verified in §4.1 and the stdlib gzip-streaming parse pattern from
`scripts/audit_vcf_content.py` (no VCF library added). Regenerate with:

```
.venv/Scripts/python.exe -m src.phylogeny.run_build_graph
```

**Pipeline** (two streaming passes over all 13,753 VCFs, SNP-only —
indels excluded, standard practice for bacterial SNP-distance phylogenetics):
pass 1 tallies every observed `(pos, alt)` site's global call count and
missingness rate; sites with missingness >10% are dropped from the **core**
set (465,257 of 465,370 total observed sites kept — the filter barely
trims anything, i.e. most sites have low missingness). Pass 2 builds a
sparse `(13753, 465257)` binary matrix (`snp_matrix.npz`). Distance and
graph construction never read `y_amr.csv`, `labels.csv`, or the `split`
column of `splits.csv` — edges are purely a function of genomic content;
`splits.csv` is read only afterward, to compute an informational diagnostic
(see below), and is never modified.

**Distance metric — Jaccard, not raw Hamming (deliberate, evidence-based
choice, not the original default)**: `d(i,j) = 1 - |shared calls| /
|union of calls|`. Two failure modes were found and fixed during
development, in order:
1. **Raw Hamming distance** (`row_sum[i] + row_sum[j] - 2*shared`) let
   isolates with an unusually *low* total call count become spurious hubs —
   one isolate (`ERR181826`, 314 calls vs. a population median of 945)
   reached degree 11,978 (87% of the population), because its distance to
   *anyone* is bounded above by its own tiny call count regardless of true
   relatedness. Jaccard normalizes by each isolate's own call count, which
   removes this artifact.
2. A **`uint8` overflow bug**: the sparse matmul used to compute shared-call
   counts accumulated in the SNP matrix's storage dtype (`uint8`, max 255)
   *before* any cast to float, so any pair sharing more than 255 calls
   silently wrapped around and produced a wrong distance. This corrupted
   the very first full-scale run (caught by re-running the correctness
   spot-check with a wider dtype and finding mismatches, *not* by luck) —
   fixed by widening to `float32`/`float64` before the dot product
   (`distance.py`). The `graph_construction_report.json.correctness_check`
   block documents a passing spot-check (0/190 mismatches, max shared-call
   count 2,218) as a permanent record that this was actually verified, not
   just patched.

**k-NN mode — union, not mutual**: strict mutual-kNN (edge kept only if
each isolate is in the other's k-nearest list) left 98.4% of nodes isolated
at full population scale — most isolates have a clear best candidate
neighbor, just not a *reciprocal* one, in a population this diverse (union
mode keeps an edge if either side nominates it).

**Result** (`k=20`, `graph_construction_report.json` has full detail):
404,944 directed edges (202,472 undirected); degree min 20 / median 26 /
mean 29.4 / max 195 — no isolated nodes, no runaway hub; 6 connected
components — one giant component covering 98.3% of isolates (13,520/13,753)
plus 5 small satellite clusters (77/56/39/33/28 nodes), a biologically
plausible structure (a few isolates/clusters more distant from the main
population than any of their `k=20` neighbors would place them in it).

**Known limitations, stated rather than hidden**:
- These are per-isolate variant-call VCFs, not joint/gVCF all-sites calls —
  a site absent from an isolate's record cannot be distinguished from
  "matches reference" vs. "no confident call" without depth data this
  dataset doesn't provide. The 10%-missingness core-site filter bounds
  this rather than eliminating it.
- **Cross-split near-duplicates**: 98 val/test isolates have a train-set
  isolate at Jaccard distance ≤0.01 (≥99% shared calls) — flagged in
  `graph_construction_report.json.cross_split_near_duplicate_diagnostic`,
  not removed or acted on. `splits.csv` stays fixed per project rules;
  whoever trains the GNN should read this before treating `test` metrics
  as a fully independent estimate in a transductive setting.
- k=20 and the 10% missingness threshold are documented defaults validated
  at ~100-isolate and full scale, not exhaustively tuned — reasonable
  starting points, not claimed-optimal.

## 13. Phase 3 — Multi-task GNN (first pass, AMR-only)

`src/gnn/` — a hard-parameter-sharing multi-task GCN: 2×`GCNConv` (157→64→64,
using `pyg_graph.pt`'s Jaccard-derived `edge_weight`) + one shared linear head
(9 logits, one per drug) — the graph analogue of `train_mlp.py`'s baseline.
AMR-only for this first pass (no `y_lineage` head yet — see §9). Transductive:
every node/edge is visible during the forward pass regardless of split (that
is standard, not leakage — edges/features never encode labels); only the
training loss is masked to `train`. Regenerate with:

```
.venv/Scripts/python.exe -m src.gnn.run_train
```

Per-drug BCE `pos_weight` from train-split class balance (same imbalance
rationale as `train_random_forest.py`'s `class_weight="balanced"`); Adam,
early stopping on val macro-F1 (patience 20). Same evaluation protocol as
Phase 1 — thresholds tuned on `val` only via
`src/baselines/metrics.py:tune_per_drug_thresholds`, applied unchanged to
`test`. Full results in `results/gnn_metrics.json`.

**Result: this first-pass GNN does not beat the Phase 1 tabular baselines
on macro-F1** — worth stating plainly rather than as a success. Test
macro-F1: GNN 0.604 → 0.637 tuned, vs. RF 0.826 → 0.869, MLP 0.762 → 0.850,
XGBoost 0.656 → 0.657 (§11). Trained in 61s (81 epochs, early-stopped ~100).

The interesting part is *where* it underperforms: per-drug **ROC-AUC is
strong and close to RF's** (test macro-ROC-AUC 0.961 vs. RF's 0.998,
clearly ahead of XGBoost's 0.900 — including on CAP/LZD, where XGBoost's
AUC was only 0.5–0.71) — the model ranks isolates well. The gap is in
**F1 at a fixed threshold** on the majority drugs (RIF/INH/EMB/PZA/STM/LEV),
where RF gets near-1.0 F1 by exploiting an almost-deterministic
mutation→resistance mapping directly, while the GCN's graph convolution
mixes in neighbors' features every layer, which plausibly smooths out
exactly the sharp per-isolate signal those drugs don't need help with.
Consistent with, not contradicted by, strong AUC: the model still separates
the classes, it's just less confidently calibrated at a single threshold.
Not chased further in this pass — candidates for a next iteration: a
skip/residual connection from raw features to the head (let the model use
graph context only where it helps, not replace direct signal), fewer/more
GCN layers, or GraphSAGE/GAT instead of GCN.

**Edge-weight ablation** (`src/gnn/run_ablation.py`, same seed/architecture/
hyperparameters, only `edge_weight` differs — `None` treats every kept edge
as weight 1 instead of its Jaccard similarity): the Jaccard weighting is
carrying real signal, not noise. Test macro-F1 **0.604→0.637 with weights vs.
0.582→0.602 without** (0.5 / tuned), macro-ROC-AUC **0.961 vs. 0.950**,
macro-PR-AUC **0.776 vs. 0.742** — weighted wins on every metric. The
unweighted run also converged much earlier (epoch 41 vs. 81) and to a lower
val macro-F1, i.e. it's not close-but-slower, it plateaus at a worse optimum.
Full detail (including per-drug) in `results/gnn_edge_weight_ablation.json`.

**Skip connection** (`src/gnn/run_skip_ablation.py`, `MultiTaskGCN(...,
skip_connection=True)` in `model.py` — concatenates each isolate's raw
157-feature vector onto the final GCN hidden state before the head, both
variants trained with `use_edge_weight=True`): confirms the oversmoothing
hypothesis above and **closes most of the gap to the tabular baselines**.
Test macro-F1 (0.5 / tuned): **0.604→0.682 without skip vs. 0.637→0.707
with skip**; macro-ROC-AUC **0.961→0.991**. Per-drug detail
(`results/gnn_skip_connection_ablation.json`) shows the gain concentrated
exactly where predicted — the majority drugs the plain GCN struggled with:
RIF 0.833→0.936, INH 0.819→0.929, EMB 0.770→0.902, PZA 0.784→0.881,
STM 0.709→0.870, LEV 0.685→0.840. One regression: LZD 0.133→0.000 (n=2
test positives — a single-example threshold-tuning artifact at this sample
size, not treated as a real signal either way, same caveat Phase 1's
baselines hit on the same rare drugs). CAP stays 1.0, ETH stays undefined
(0 test positives). Trained longer before early-stopping (epoch 181 vs. 81,
167s vs. 83s) — direct feature access gives the optimizer more useful
gradient to keep improving on.

**Current best GNN vs. Phase 1** (test macro-F1, tuned threshold): **GCN +
skip 0.707** vs. RF 0.869, MLP 0.850, XGBoost 0.657 — narrowed from a ~23-point
gap to ~16 points, and now clearly ahead of XGBoost. `src/gnn/run_train.py`
still defaults to `skip_connection=False`; adopting the skip connection as
the default (and re-running the primary `results/gnn_metrics.json`) is a
reasonable next step but wasn't done automatically here since it changes
what "the" GNN result means going forward.
