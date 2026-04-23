# ProteinMPNN Fine-tuning — pMHC Data Preparation

## Citation

This data preparation pipeline was developed for **PMGen: From Peptide-MHC Structure Prediction to Peptide Generation**.


Pipeline repository: https://github.com/hasmiig/proteinmpnn-finetune-data-prep.git

---

End-to-end pipeline for curating a high-quality dataset of pMHC class I binder structures for fine-tuning ProteinMPNN. Starting from raw IEDB exports, the pipeline produces balanced, structure-predicted, QC-filtered datasets in ProteinMPNN-ready JSONL format with MHC chain fixed and peptide chain designable.

For full methodology and design decisions, see [docs/METHODOLOGY.md](docs/METHODOLOGY.md).

---

## Overview

```
IEDB export (binders only)
    │
    ▼
01_data_preparation/      ← filter to binders, deduplicate, balance alleles
    │
    ▼
02_structure_prediction/  ← prepare PMGen input, run structure prediction
    │
    ▼
03_filtering_analysis/    ← pLDDT QC, re-balance, produce ProteinMPNN JSONL splits
    │
    ▼
data.jsonl + fixed_positions.json  (per train/val/test fold)
```

---

## Quick Start

If you have a clean IEDB export or similar pMHC class I dataset, here's the fastest path to ProteinMPNN-ready fine-tuning data:

1. **Stage 1: Prepare & Balance** — Filter binders, balance alleles
   ```bash
   cd scripts/01_data_preparation/
   python pmhc_sampling.py --input your_data.parquet --mode binder --phases only_phase1 --output binders_sampled.parquet --plots results/
   ```
   **Output:** Balanced dataset (~100–200k peptides), ready for structure prediction

2. **Stage 2: Predict Structures** — Run PMGen (on a GPU cluster)
   - Prepare input: `python prepare_pmgen_input.py`
   - Chunk for parallel jobs: `python chunk_tsv.py --input pmgen_input.tsv --size 4000`
   - Submit jobs: `for chunk in chunks/*.tsv; do sbatch run_pmgen.sh $chunk; done`
   - **Output:** PDB files in `outputs/`

3. **Stage 3: Filter & Split** — Apply pLDDT QC, rebalance, generate train/val/test
   ```bash
   cd scripts/03_filtering_analysis/
   python compute_plddt_means.py --mode binder --base /path/to/outputs
   python analyse_plddt.py  # Edit MODE, BASE, PLOTS_BASE constants at top of script
   python filter_map_train_prep.py --plddt_csv outputs/binder/plddt_means_binder.csv \
     --parquet binders_sampled.parquet --pdb_base_dir outputs/ \
     --output_dir splits/ --split_mode hla
   python prepare_pmhc_data.py --splits_dir splits/hla --output_dir proteinmpnn_input/
   ```
   **Output:** ProteinMPNN-ready JSONL files in `proteinmpnn_input/`

For detailed options and troubleshooting, see the stage-specific READMEs in `scripts/`.

---

## Repository Structure

```
proteinmpnn-finetune-data-prep/
├── scripts/
│   ├── 01_data_preparation/
│   │   ├── pmhc_sampling.py         # Iterative median sampling (allele balancing)
│   ├── 02_structure_prediction/
│   │   ├── prepare_pmgen_input.py   # Build TSV input for PMGen
│   │   └── chunk_tsv.py             # Split large files for parallel processing
│   │   └── run_pmgen.sh             # SLURM job script for PMGen structure prediction
│   ├── 03_filtering_analysis/
│   │   ├── compute_plddt_means.py   # Extract mean pLDDT per structure
│   │   ├── analyse_plddt.py         # Visualize pLDDT distributions, apply threshold
│   │   ├── filter_map_train_prep.py # Filter, re-sample, produce train/val/test splits
│   │   └── prepare_pmhc_data.py     # Convert split parquets → ProteinMPNN JSONL
│   └── utils/
│       ├── setup_env.sh             # Environment setup for cluster jobs
│       └── submit_prepare_pmhc.sh   # Batch submission helper
├── docs/
│   ├── METHODOLOGY.md               # Full design rationale and statistics
│   └── <dated run logs>/            # Per-run exploration outputs
└── README.md
```

---

## Requirements

```
pandas
numpy
matplotlib
scipy
pyarrow     # for .parquet files
```

```bash
pip install pandas numpy matplotlib scipy pyarrow
```

Structure prediction requires PMGen and a GPU node (see `scripts/utils/setup_env.sh`).

---

## Input Data Setup

### Directory Structure

Before starting, organize your input files as follows:

```
<your_project_root>/
├── data/
│   ├── raw/
│   │   └── your_iedb_export.parquet          # Stage 1 input
│   └── mhc1_encodings.csv                    # Stage 2 requirement: MHC sequences
└── proteinmpnn-finetune-data-prep/
    └── scripts/
        ├── 01_data_preparation/
        ├── 02_structure_prediction/
        └── 03_filtering_analysis/
```

### Configuring Script Paths

Several scripts have hardcoded paths at the top. Before first use, update these templates:

**`scripts/02_structure_prediction/prepare_pmgen_input.py`**
```python
# Line ~10–15 (update these paths)
PARQUET_PATH = "/path/to/binders_sampled.parquet"
MHC_ENCODINGS_PATH = "/path/to/data/mhc1_encodings.csv"
OUTPUT_PATH = "pmgen_input.tsv"
```

**`scripts/03_filtering_analysis/filter_map_train_prep.py`**
```python
# Line ~20–30 (set defaults or pass via CLI args)
PLDDT_CSV = "outputs/binder/plddt_means_binder.csv"
PARQUET_PATH = "binders_sampled.parquet"
MHC_ENCODINGS_PATH = "data/mhc1_encodings.csv"
PDB_BASE_DIR = "outputs/binder/"
```

Alternatively, all parameters can be passed as command-line arguments — see `python script_name.py --help` for options.

### Required Columns in Input Data

If using a different data source, ensure your `.parquet` or `.csv` has:

| Column | Type | Example |
|--------|------|----------|
| `long_mer` | string | `LLFGYTWP` |
| `allele` | string | `HLA-A*02:01` |
| `assigned_label` | float | `1.0` (binder) or `0.0` (non-binder) |
| `mhc_class` | float | `1.0` (MHC-I) |

If your source uses different column names, update the constants at the top of `scripts/01_data_preparation/pmhc_sampling.py`:

```python
COL_PEPTIDE = "long_mer"           # your peptide sequence column
COL_MHC     = "allele"             # your MHC allele column
COL_LABEL   = "assigned_label"     # your binding label column
MHC_CLASS   = "mhc_class"          # your MHC class column
```

---

## Scope

This project focuses on **MHC-I binder** data for ProteinMPNN fine-tuning.
Non-binder support is present throughout the codebase (via `--mode nonbinder`)
but was developed for a separate project and is not the focus of this pipeline.

---

## Stage 1 — Data Preparation

**Directory:** `scripts/01_data_preparation/`

**Methodology:** See [docs/2026-03-25_subsampling/README.md](docs/2026-03-25_subsampling/README.md) for detailed rationale on two-phase sampling strategy.

### 1a. Inspect raw IEDB export

Verify data structure and columns before processing:

```bash
python pmhc_sampling.py --input iedb_export.csv --inspect
```

Prints schema and first rows. Verify column names match the constants at the top of `pmhc_sampling.py` (see [Input column requirements](#input-column-requirements) below).

### 1b. Explore raw binder distributions (optional)

Understand allele and anchor composition before sampling:

```bash
python pmhc_sampling.py \
    --input iedb_export.csv \
    --mode binder \
    --explore \
    --plots exploration_raw/
```

### 1c. Filter, deduplicate, balance, and generate plots (**REQUIRED**)

After inspection/exploration, perform the full data preparation:

```bash
python pmhc_sampling.py \
    --input iedb_export.csv \
    --mode binder \
    --phases only_phase1 \
    --output binders_sampled.parquet \
    --plots results/binders/
```

This single command performs:

1. **Filter** — Keep MHC class I binders only (`assigned_label == 1.0`, `mhc_class == 1.0`)
2. **Deduplicate** — Remove exact duplicates by (peptide, allele) pair
3. **Validate** — Remove peptides with non-standard amino acids or length ≥ 15 residues
4. **Phase 1 Balancing** — Iterative median sampling to balance allele representation (cap: 1000 peptides per allele)
5. **Generate Plots** — Summary statistics and allele/anchor distributions to `results/binders/`

**Phase 1 only (not Phase 2):** Anchor residue preferences in binders reflect biological MHC-peptide constraints and should not be down-sampled. Phase 2 anchor balancing is only applied to non-binder datasets.

**Output:** `binders_sampled.parquet` — ready for Stage 2 (structure prediction).

### Input column requirements

| Column | Description |
|---|---|
| `long_mer` | Peptide amino acid sequence |
| `allele` | MHC allele name (e.g. `HLA-A*02:01`) |
| `assigned_label` | Binding label: `1.0` = binder, `0.0` = non-binder |
| `mhc_class` | MHC class: `1.0` = class I, `2.0` = class II |

If your file uses different column names, update the constants at the top of `pmhc_sampling.py`:

```python
COL_PEPTIDE = "long_mer"
COL_MHC     = "allele"
COL_LABEL   = "assigned_label"
MHC_CLASS   = "mhc_class"
SAMPLE_CAP  = 1000
```

---

## Stage 2 — Structure Prediction

**Directory:** `scripts/02_structure_prediction/`

### 2a. Prepare PMGen input

```bash
python prepare_pmgen_input.py
```

Merges sampled peptides with MHC sequences from `mhc1_encodings.csv` and writes a TSV for PMGen. Update the hardcoded paths at the top of the script to point to your parquet and MHC encoding files.

### 2b. Chunk TSV for parallel processing

Split the large TSV into smaller chunks for parallel structure prediction:

```bash
python chunk_tsv.py --input pmgen_input.tsv --size 4000 --output chunks/
```

Creates multiple `chunk_*.tsv` files in the `chunks/` directory (adjust `--size` based on available GPU memory and job duration preferences).

### 2c. Run PMGen on the cluster

```bash
for chunk in chunks/*.tsv; do
    sbatch run_pmgen.sh $chunk
done
```

Each job runs PMGen in `--initial_guess --multiple_anchors` mode, generating two structures per pMHC pair. PDB outputs are written to `outputs/<chunk>/`.

---

## Stage 3 — Filtering, Splitting, and ProteinMPNN Format Conversion

**Directory:** `scripts/03_filtering_analysis/`

### 3a. Extract pLDDT means

```bash
python compute_plddt_means.py --mode binder --base /path/to/outputs
```

Reads PMGen PDB outputs and writes `plddt_means_binder.csv` with per-structure mean peptide pLDDT and anchor pLDDT.

**Flags:**
- `--mode` — `binder` or `nonbinder` (required)
- `--base` — Path to outputs folder containing chunk_* directories (required)
- `--output` — Directory to save output CSV (optional; defaults to `base/mode/`)

### 3b. Analyse pLDDT distributions

```bash
python analyse_plddt.py
```

Generates pLDDT histograms, applies the 80-threshold filter, and saves the best structure per (allele, peptide) pair to `plddt_means_binder_best.csv`.

**Note:** This script uses hardcoded constants at the top. Before running, update:
```python
MODE = "binder"                          # or "nonbinder"
BASE = Path("/path/to/outputs")          # Base outputs directory
PLOTS_BASE = Path("/path/to/plots")      # Directory to save plots
```

### 3c. Filter, re-sample, and split

```bash
python filter_map_train_prep.py \
    --plddt_csv      outputs/binder/plddt_means_binder.csv \
    --parquet        iedb_mhc1_binders.parquet \
    --mhc_encodings  data/mhc1_encodings.csv \
    --pdb_base_dir   outputs/binder \
    --output_dir     trainprep/binder/ \
    --mode           binder \
    --plddt_threshold 80 \
    --k              5 \
    --split_mode     hla
```

Four sequential stages:

1. **Filter** — select best structure per (allele, peptide); apply pLDDT threshold
2. **Map** — join back to original parquet to recover metadata, MHC sequences, and paths to PDB files
3. **Resample** — re-run iterative median sampling to correct any allele bias from prediction
4. **Split** — train/val/test splits in one of two modes:
   - `hla` — rare alleles (bottom 20% by frequency) held out as test; k-fold CV on remainder
   - `anchor` — k-fold CV stratified by anchor residue combinations (P2 × C-terminal)

### 3d. Convert to ProteinMPNN JSONL format

```bash
python prepare_pmhc_data.py \
    --splits_dir trainprep/binder/splits/hla \
    --output_dir proteinmpnn_input/binder_hla \
    --split_mode hla
```

Reads each split parquet, finds the corresponding PDB files, and writes:
- `data.jsonl` — ProteinMPNN-format structure records
- `fixed_positions.json` — chain A (MHC) fixed, chain P (peptide) designable

Output structure (HLA mode):
```
proteinmpnn_input/binder_hla/
  test/
    data.jsonl
    fixed_positions.json
  fold_1/
    train/  val/
  fold_2/ ... fold_5/
```

---

## Key Numbers

| Stage | Count |
|---|---|
| Initial IEDB MHC-I binders | ~1,100,000 |
| After allele balancing | ~163,000 |
| After structure prediction (pLDDT ≥ 80) | 87,187 |
| After final re-balancing | **63,817** |
| Unique MHC-I alleles | 426 |

---

## Validating Your Run

### Stage 1 — Data Preparation

**Expected outputs** in the `--plots` directory:
```
results/binders/
├── allele_distribution.pdf          # Bar chart: # peptides per allele
├── anchor_heatmaps_*.pdf            # P2 × C-terminal anchor combinations
├── sequence_length_dist.pdf         # Peptide length distribution
└── summary_stats.txt                # Row counts before/after filtering
```

**Expected parquet size:**
- Input: ~1M rows (raw IEDB)
- Output: ~100–200k rows (after balancing)

**Quick check:**
```python
import pandas as pd
df = pd.read_parquet("binders_sampled.parquet")
print(f"Rows: {len(df)}, Alleles: {df['allele'].nunique()}")
# Expected: Rows: ~163000, Alleles: 459
```

---

### Stage 2 — Structure Prediction

**Expected outputs** in `outputs/`:
```
outputs/chunk_1/
├── <peptide>_<allele>_0.pdb        # Structure 1
├── <peptide>_<allele>_1.pdb        # Structure 2 (best structure kept later)
└── ... (many more .pdb files)
```

**Expected counts:**
- Input TSV: ~163k peptide–allele pairs (Stage 1 output)
- Output PDBs: ~326k files (2 per pair)
- Expected runtime: ~1–3 days on GPU cluster (size + hardware dependent)

**Failure indicators:**
- Fewer than ~10k PDBs across all chunks → check PMGen logs, GPU availability, or anchor enumeration

---

### Stage 3 — Filtering & Splitting

**Expected outputs** in `proteinmpnn_input/`:

```
proteinmpnn_input/binder_hla/
├── test/
│   ├── data.jsonl                  # Hold-out test (rare alleles)
│   └── fixed_positions.json
├── fold_1/
│   ├── train/
│   │   ├── data.jsonl              # 60% of data
│   │   └── fixed_positions.json
│   └── val/
│       ├── data.jsonl              # 20% of data
│       └── fixed_positions.json
├── fold_2/ ... fold_5/             # Repeat for 5-fold CV
```

**Expected size:**
- Total records across train/val/test: **~60–65k** (after pLDDT filtering + re-balancing)
- Per-fold train size: ~24–26k records
- Per-fold val size: ~8–10k records

**Quick check:**
```bash
# Count records in train JSONL (should be ~24k)
wc -l proteinmpnn_input/binder_hla/fold_1/train/data.jsonl

# Verify fixed_positions.json exists and contains chain A only
python -c "import json; print(json.load(open('proteinmpnn_input/binder_hla/fold_1/train/fixed_positions.json')))"
```

**Common issue:** If counts are much lower than expected, check pLDDT threshold and allele filtering; see `filter_map_train_prep.py` logs.

---

## Documentation

- [docs/METHODOLOGY.md](docs/METHODOLOGY.md) — rationale for binders-only approach,
  sampling strategy, pLDDT thresholds, splitting design
- [`docs/<date>/`](docs/) — per-run exploration plots and stats
- [`scripts/01_data_preparation/README.md`](scripts/01_data_preparation/README.md) — `pmhc_sampling.py` usage and workflow
- [`scripts/03_filtering_analysis/README.md`](scripts/03_filtering_analysis/README.md) — `filter_map_train_prep.py` usage and workflow