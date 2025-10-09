# 01_GetCounts

**Purpose**  
Combine per-sample `.featureCounts` outputs into a single count matrix for downstream analysis (DE, modeling, figures).

## Expected inputs
- One `.featureCounts` file per sample produced by Subread **featureCounts v2.0.6**.
- All files must be generated with the *same* reference and counting options.
- An experiment design file `Input/ExperimentDesign_batch_Z01-F001_F002_combined_update.csv`.

### Reference & counting contract (required)
- Genome: **Nfu_20140520**
- Annotation: **GCF_001465895.1** 


## Outputs (gitignored)
- `Output/Counts_Atlas_Plate01-02_240319.csv` — raw count matrix
- `Output/TPM_Atlas_Plate01-02_240319.csv` — tpm
- `Output/GeneLengths_Atlas_Plate01-02_240319.csv` — gene lengths

## Run
- Run 1_getCounts_platescombo.R in RStudio.
- To run 2_getCombinedCountMatrix_220805_codechk_platescombo.pl: 
```bash
perl 2_getCombinedCountMatrix_220805_codechk_platescombo.pl
```
