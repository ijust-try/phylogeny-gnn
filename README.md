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
- `AFRO_TB_VCF/` — 13,753 per-isolate VCFs. **Not yet processed by anything.**
- `AFRO_TB_dataset/` — per-isolate annotation `.txt` files. **Not yet
  processed by anything.**

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
- **Sample-ID namespace**: `features.csv`/`labels.csv` IDs (e.g.
  `ERR6397155`) have not yet been explicitly cross-validated against
  `AFRO_TB_VCF/` filenames (e.g. `ERR036186_MT.vcf.gz`) or
  `AFRO_TB_dataset/*.txt` — confirm this mapping before building the
  phylogenetic distance matrix from VCFs.

## 9. What Person 2 owns

Per `CLAUDE.md`: Graph Neural Network architecture, multi-task learning,
model training, validation, model evaluation. Concretely, next steps on top
of what Person 1 has delivered:

- Build the SNP-distance / phylogenetic representation from
  `AFRO_TB_VCF/*.vcf.gz`, keyed to `sample_ids.csv`'s canonical IDs (confirm
  the ID-namespace mapping first, §8).
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
