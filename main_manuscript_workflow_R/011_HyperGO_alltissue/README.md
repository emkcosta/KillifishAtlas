# 011_HyperGO_alltissue

**Purpose**  
Hypergeometric GO enrichment analysis all tissues

## Expected inputs
- "241208_allCorrresults", t ,"_", sex,"Only.csv": age-correlated genes for each tissue and sex
- 'GO_killifish-human_best-hits.txt': orthologs in human for killifish genes

## Outputs (gitignored)
- This is where all of the output tables and some plots will be stored.


## Run
- 4_HyperGO_AllTissue_GOterms_SexSplitCorrelation_forSexDivergentTerms_250521.R - Hypergeometric GO enrichment analysis focusing on sex-divergent terms
- 4_HyperGO_AllTissue_GOterms_SexSplitCorrelation_forSexSharedTerms_250504.R - Hypergeometric GO enrichment analysis focusing on sex-shared terms
- 5_HyperGO_BubblePlot_AllTissue_SexDivergentTerms_250507.R - Plot the bubble plots for selected GO terms, focusing on sex-divergent terms
- 5_HyperGO_BubblePlot_AllTissue_SexSharedTerms_250520.R - Plot the bubble plots for selected GO terms, focusing on sex-shared terms