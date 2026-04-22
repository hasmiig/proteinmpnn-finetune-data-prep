# PMGen Input Preparation

## Date
2026-03-27

## Analyst
Hasmig Aintablian

## What we did
Prepared input TSV files for PMGen predictions from the sampled binder 
parquet file. Two steps:

1. Merged the sampled peptides with MHC sequences from `mhc1_encodings.csv` to 
   produce a PMGen-compatible TSV with columns: peptide, mhc_seq, mhc_type, anchors, id.
2. Split the TSV into chunks of 4000 rows each for parallelized GPU job submission.

## Why we did it
PMGen requires a specific TSV format as input. Chunking allows each chunk to be 
submitted as an independent SLURM job, enabling parallelization and easier 
resubmission of failed jobs.

## Results
- 163,948 rows → 41 chunks of 4000 rows each


## Conclusion
Input chunks ready for PMGen structure prediction.