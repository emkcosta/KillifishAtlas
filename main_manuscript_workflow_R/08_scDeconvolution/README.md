# 08_scDeconvolution

**Purpose**  
Single cell deconvolution for kidney

## Expected inputs
- 'kidneyOnly_Teefy_SeurObj_Fordeconvolution.Rdata': Seurat object from single cell reference dataset Teefy et al
- 'GeneLengths_Atlas_Plate02_Lane03-4_BatchF004_240417.csv': gene lengths for killifish genes



## Outputs (gitignored)
- This is where all of the output tables and some plots will be stored.


## Run
- 01_singlecell_mat_granulator.R - Creates a single cell matrix of the kidney data in TPM
- 02_Kidney_deconvolution_granulator.R - Performs the scdeconvolution for the kidney using granulator
- 03_Kidney_deconvolution_granulator_plotting.R - Plots scdecon results 
- 04_singlecell_mat_granulator_grouped.R - Creates a single cell matrix of the kidney data in TPM, using cell type groups
- 05_Kidney_deconvolution_granulator_grouped_plotting_.R - Performs the scdeconvolution for the kidney using granulator, using cell type groups
