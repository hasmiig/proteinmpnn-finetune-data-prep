# Subsampling — pMHC Binder Dataset

## Date
2026-03-25

## Analyst
Hasmig Aintablian

## What we did
Applied a two-phase iterative median sampling pipeline (`pmhc_sampling.py`) to balance 
the raw pMHC class I dataset from Amirreza's database (PMDb_2025_11_18_class1.parquet).
Evaluated both Phase 1 (MHC allele balancing) and Phase 2 (anchor residue diversity balancing),
but proceeded with **Phase 1 only** for fine-tuning.

## Why we did it
The raw dataset is heavily imbalanced. Some alleles have hundreds of thousands of 
peptides while others have only a handful. This would introduce bias into any downstream 
model. Phase 1 balances MHC allele representation to ensure fair training data.

## Results

- Raw: 1,137,039 peptides across 459 alleles
- Post-Phase 1: 163,948 peptides across 459 alleles
- Cap: 1000 peptides per allele
- 128 alleles at cap, 331 below cap (too few peptides)

## Conclusion
Binder dataset reduced from 1.1M to 163k peptides. Dataset is now balanced 
across alleles and ready for PMGen input preparation.

## GitHub
Pipeline and full documentation available in the scripts/ directory.