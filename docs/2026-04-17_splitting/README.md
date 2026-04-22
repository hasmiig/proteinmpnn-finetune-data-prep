# Training Data Preparation & Splitting — Binder Dataset

## Date
2026-04-17

## Analyst
Hasmig Aintablian

## What we did
Prepared final training data from filtered binder structures using a four-stage pipeline:

1. **FILTER** — Selected best structure per (allele, peptide) pair based on pLDDT quality
   (threshold: 80 for peptide mean pLDDT)

2. **MAP** — Joined filtered structures back to raw IEDB parquet data to recover:
   - Peptide sequences and binding labels
   - MHC allele annotations
   - MHC protein sequences from encodings

3. **RESAMPLE** — Applied Phase 1 iterative median sampling to balance MHC allele 
   representation. Anchor dominance in binders is biologically constrained and was 
   not further balanced in Phase 2.

4. **SPLIT** — Analyzed split strategies (HLA-stratified and anchor-stratified) 
   to understand optimal cross-validation design for downstream ProteinMPNN fine-tuning.

## Why we did it
The raw structures had extreme allele imbalance — some alleles dominated the dataset. 
Resampling ensures fair representation across all MHC alleles for model training. 
Analysis of splitting strategies informs which validation scheme best represents 
unseen data diversity.

## Results

### Stage-by-Stage Data Flow

| Stage | Rows | Alleles | Median rows/allele | Max rows/allele |
|-------|------|---------|-------------------|-----------------|
| After pLDDT filter | 163,948 | 459 | — | — |
| After map + filter | 87,187 | 426 | 53 | 876 |
| After Phase 1 resample | 63,817 | 426 | 53 | 376 |

### Key Findings

#### Allele Representation
- **Mapped dataset**: Median 53 peptides/allele; 426 alleles present (none at cap of 1,000)
- **After Phase 1**: Median maintained at 53; max capped at 376 (from 876)
- **Resampling effect**: Removed 23,370 rows (27% reduction); no allele loss

#### Peptide Length Distribution (Post-Phase 1)
- 9-mer: 49,608 peptides (77.7%) — most abundant
- 8-mer: 6,530 (10.2%)
- 10-mer: 5,192 (8.1%)
- 11-15-mer: <500 each
- All peptides < 15 residues (valid)

#### Anchor Residue Composition

**Position 2 (P2) frequencies post-Phase 1:**
- Most common: P (proline, 7,749), L (leucine, 7,250), E (glutamate, 6,329)
- Balanced across 20 amino acids with biological skew toward P, L, E, R, T

**C-terminal (last position) frequencies post-Phase 1:**
- Strongly biased: L (leucine, 19,675), F (phenylalanine, 8,809), Y (tyrosine, 7,265)
- MHC class I anchor requirement: aliphatic or aromatic (L, F, V, I, M, Y, W preferred)
- Only 87/63,817 (0.1%) non-canonical C-terminals (D, N, S, etc.)

### Diagnostic Plots Generated

- `plddt_distribution_before_after.png` — Peptide and anchor pLDDT histograms showing quality improvement
- `allele_distribution_after_mapping.png` — Per-allele peptide counts (before resampling)
- `allele_distribution_after_resampling.png` — Per-allele peptide counts (after Phase 1)
- `comparison_allele_distribution.png` — Allele counts across stages (mapped vs post-Phase 1)
- `comparison_peptide_lengths.png` — Peptide length distribution across stages
- `comparison_anchor1_residues.png` — P2 amino acid frequency comparison
- `comparison_anchor2_residues.png` — C-terminal amino acid frequency comparison
- `per_allele_anchor1_postphase_1.png` — P2 residue diversity within each of 426 alleles
- `per_allele_anchor2_postphase_1.png` — C-terminal residue diversity within each of 426 alleles
- `anchor_combo_heatmap.png` — Heatmap of (P2 + C-term) anchor pair frequencies
- `anchor_kl_divergence_boxplot.png` — KL divergence of anchor distributions (stability metric)
- `anchor_combo_stats.csv` — Per-allele unique anchor combinations and top combo ratio

## Conclusion
Binder training data successfully prepared with:
- **63,817 peptides** across **426 alleles** (balanced representation)
- **Strong anchor constraints** preserved (biologically appropriate)
- **High pLDDT quality** post-filtering (≥80)
- **Ready for ProteinMPNN fine-tuning** with train/val/test splits

## Files
- `plots/binder/` — All diagnostic plots and CSV summaries
  - `stats_report.txt` — Detailed per-stage statistics
  - `stats_summary.csv` — Condensed summary table
  - `anchor_combo_stats.csv` — Per-allele anchor pair analysis
