# 03_DESeq2

**Purpose**  
DESeq2 analyses

## Expected inputs
- 'Counts_Atlas_allbatches_merged_v3.csv': raw counts
- 'ExperimentDesign_allbatches_combined_v7.csv': metadata
- 'TPM_Atlas_allbatches_merged_v3.csv': TPM


## Outputs (gitignored)
- This is where all of the output tables and some plots will be stored.


## Run
- 01_DESeq2_EC.R - this is where we generate the list of Deseq2 objects by tissue and perform sex combined age-related DE analysis
- 02_DESeq2_mvsf_eachtimepoint.R - perform DE analysis in each time bin between the sexes, for each tissue
- 03_DESeq2_sepbysex.R - perform DE analysis across time, separated by sex
- Fig1e_DESeq2_mvsf_eachtimepoint.R - plotting for Fig1e

