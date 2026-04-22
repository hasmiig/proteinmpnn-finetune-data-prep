# Training Data Pipeline — `filter_map_train_prep.py`

This directory contains the main pipeline script for preparing ProteinMPNN fine-tuning
training data from PMGen pMHC structure predictions.

---

## Overview

`filter_map_train_prep.py` runs four sequential stages:

1. **FILTER** — Select the best structure per (allele, peptide) pair based on pLDDT; apply confidence threshold
2. **MAP** — Join filtered structures back to raw IEDB parquet to recover peptide metadata and MHC sequences
3. **RESAMPLE** — Re-run iterative median sampling to correct residual allele bias introduced during structure prediction
4. **SPLIT** — Produce train/val/test splits in HLA-stratified or anchor-stratified mode

> **Dependency:** This script imports plotting and sampling utilities from `pmhc_sampling.py`
> (located in `scripts/01_data_preparation/pmhc_sampling.py` directory). Ensure that script is accessible before running.

---

## Usage

### Binders — HLA-stratified split

```bash
python filter_map_train_prep.py \
    --plddt_csv      path/to/plddt_means_binder.csv \
    --parquet        path/to/PMDb_class1.parquet \
    --mhc_encodings  path/to/mhc1_encodings.csv \
    --pdb_base_dir   path/to/outputs/binder \
    --output_dir     path/to/trainprep/binder/ \
    --mode           binder \
    --plddt_threshold 80 \
    --k              5 \
    --split_mode     hla
```

### Non-binders — Anchor-stratified split (keep both structures)

```bash
python filter_map_train_prep.py \
    --plddt_csv      path/to/plddt_means_nonbinder.csv \
    --parquet        path/to/PMDb_class1.parquet \
    --mhc_encodings  path/to/mhc1_encodings.csv \
    --pdb_base_dir   path/to/outputs/nonbinder \
    --output_dir     path/to/trainprep/nonbinder/ \
    --mode           nonbinder \
    --keep_all \
    --plddt_threshold 0 \
    --k              5 \
    --split_mode     anchor
```

---

## Key Options

| Flag | Required | Description |
|------|----------|-------------|
| `--plddt_csv` | Yes | CSV of per-structure pLDDT scores (output of pLDDT QC step) |
| `--parquet` | Yes | Raw IEDB parquet file with peptide and binding metadata |
| `--mhc_encodings` | Yes | CSV mapping allele names to MHC protein sequences |
| `--pdb_base_dir` | No | Base directory of PMGen PDB outputs; used to construct `pdb_path` column |
| `--output_dir` | Yes | Directory where all outputs will be saved |
| `--mode` | Yes | `binder` or `nonbinder` — controls filtering and resampling behavior |
| `--keep_all` | No | *(nonbinder only)* Retain both predicted structures per pair instead of selecting the best |
| `--plddt_threshold` | No | Minimum `pep_mean_plddt` to retain a structure (default: `80`; set `0` to disable) |
| `--k` | No | Number of cross-validation folds (default: `5`) |
| `--split_mode` | Yes | `hla`, `anchor`, or `both` |

---

## Inputs Required

| File | Description |
|------|-------------|
| `plddt_means_{mode}.csv` | Output of pLDDT QC step; must contain columns `allele`, `peptide`, `chunk`, `struct_idx`, `pep_mean_plddt`, `anchor_mean_plddt` |
| `PMDb_class1.parquet` | Raw IEDB export; must contain `long_mer`, `allele`, `assigned_label`, `mhc_class` |
| `mhc1_encodings.csv` | Allele-to-sequence mapping; must contain `key` and `mhc_sequence` columns |

---

## Binder vs Non-binder Mode

| Behaviour | Binder | Non-binder |
|-----------|--------|------------|
| Structure selection | Best per pair (highest mean pLDDT) | Best per pair, or both if `--keep_all` |
| pLDDT threshold | 80 (recommended) | 0 (disabled, recommended) |
| Resampling | Phase 1 only (allele balancing) | Phase 1 + Phase 2 (allele + anchor balancing) |

Anchor residue balancing (Phase 2) is skipped for binders because anchor preferences
reflect genuine MHC groove chemistry — homogenising them would remove biologically
meaningful signal.

---

## Split Modes

### `hla` — Allele-based split
- Holds out the rarest 20% of alleles as an independent test set
- Remaining alleles are divided into `k` folds for cross-validation
- **Use case:** Evaluate generalisation to unseen or rare MHC alleles

Output structure:
```
splits/hla/
    test.parquet
    test_alleles.csv
    fold_summary.csv
    split_meta.json
    fold_1/train.parquet, val.parquet
    fold_2/train.parquet, val.parquet
    ...
```

### `anchor` — Anchor-combination-based split
- Splits on unique anchor combinations (P2 × C-terminal residue pairs)
- Rotating fold assignment: test = chunk i, val = chunk i+1, train = remainder
- Every anchor combination is held out as test exactly once across `k` folds
- **Use case:** Evaluate generalisation to unseen anchor contexts within known alleles

Output structure:
```
splits/anchor/
    fold_summary.csv
    fold_1/train.parquet, val.parquet, test.parquet, test_combos.csv
    fold_2/...
    ...
```

### `both`
Runs both split strategies in a single execution.

---

## Outputs

All outputs are saved under `--output_dir`:

```
output_dir/
    mapped.parquet              # After Stage 2 (map)
    resampled.parquet           # After Stage 3 (resample) — input to Stage 4
    pipeline.log                # Full run log with per-stage counts
    pipeline_summary.txt        # Human-readable summary of args and row counts
    plots/
        plddt_distribution_before_after.png
        allele_distribution_after_mapping.png
        allele_distribution_after_resampling.png
        comparison_allele_distribution.png
        comparison_peptide_lengths.png
        comparison_anchor1_residues.png
        comparison_anchor2_residues.png
        anchor_combo_heatmap.png
        anchor_kl_divergence_boxplot.png
        stats_report.txt
        stats_summary.csv
        anchor_combo_stats.csv
    splits/
        hla/   ...              # If --split_mode hla or both
        anchor/ ...             # If --split_mode anchor or both
```

---

## Dependencies

```
pandas
numpy
matplotlib
pyarrow
```

Also requires `pmhc_sampling.py` from the parent `scripts/01_data_preparation/pmhc_sampling.py` directory
for diagnostic plots. If not found, the pipeline still runs but skips
the comparison plots and stats reports.

Install with:
```bash
pip install pandas numpy matplotlib pyarrow
```

---

## Next Step

The `splits/` directory produced by this script is the direct input to
`prepare_pmhc_data.py` (in `03_filtering_analysis/`), which converts
the parquet splits into ProteinMPNN-ready `data.jsonl` and
`fixed_positions.json` files.