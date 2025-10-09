# 010_Batch_effects

**Purpose**  
Assessing batch effects.

## Expected inputs
- 'dds_TPS_allsamples_Gonadcombo_240714.bin': List of DESeq2 objects for all tissues
- 'TPM_Atlas_allbatches_merged_v3.csv': TPM
- 'ExperimentDesign_allbatches_combined_v7.csv': metadata
- '241121_Alltissue_variancePartition_results_Gonadcombo.csv': tissue-specific results from variancePartition


## Outputs (gitignored)
- This is where all of the output tables and some plots will be stored.


## Run
- 01_PCA_batchanalysis.R - plot tissue PCA colored by age, sex, cohort, RNA  batch
- 02_VariancePartition_assessbatcheffects.R - compute mediance variance explained by factor across tissues; bar plots
