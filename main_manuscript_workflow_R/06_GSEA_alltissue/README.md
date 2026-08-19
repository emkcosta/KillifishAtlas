# 06_GSEA_alltissue

**Purpose**  
GSEA GO enrichment analysis all tissues

## Expected inputs
- "241208_allCorrresults", t ,"_", sex,"Only.csv": age-correlated genes for each tissue and sex
- 'NCBI-Human-orthologs.txt': orthologs in human for killifish genes


## Outputs (gitignored)
- This is where all of the output tables and some plots will be stored.


## Run
- 4_GSEA_AllTissue_GOterms_SexSplitCorrelation_241209.R - GSEA 
- 5_GSEA_AllTissue_SexDivergentTerms_250108.R - Bubble plots for selected GO terms that show opposite directions of change for M/F
- 5_GSEA_ExampleTissue_BubblePlot_250108.R - Plot the bubble plots for selected GO terms
- 6_GSEA_ExampleTissue_Heatmap_SelectGOterms_250108.R - Plot the heatmap for selected GO terms