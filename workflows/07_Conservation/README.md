# 07_Conservation

**Purpose**  
Conservation analysis

## Expected inputs
- 'nfur-ncbi_orthologs_20170302_pps.csv': orthologs
- 'SupFile_4_CorrResults_SexCombined.xlsx': supplemental file 4 from this publication
- 'dds_DESeq2_TPS_allsamples_agesbinned_240805.bin': DESeq2 results for killifish by tissue
- 'CAS_geneset_Ohahn.xlsx': common aging score genes from Hahn et al
- 'GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.csv': from GTEx
- GTEx v10 '.gct' files for all tissues
- 'dds_TPS_alltissue_femalesOnly_DESeq2results_241106.bin': DESeq2 results for killifish by tissue, female
- 'dds_TPS_alltissue_malesOnly_DESeq2results_241106.bin': DESeq2 results for killifish by tissue, male
- 'SupFile_3_CorrResults_SexSplit.xlsx': supplemental file 3 from this publication
- 'NCBI-Human-orthologs.txt': human orthologs


## Outputs (gitignored)
- This is where all of the output tables and some plots will be stored.


## Run
- 01_mouse_xspec.R - Compare DEGs in mouse + killifish dataset
- 02_mouse_casagingscore.R - CAS aging score genes 
- 03_mouse_DEGdist_plot.R - Plot gene boxplots by physiological system
- 04_mouse_GSEA_killifishDEGs.R - Performs GSEA using mouse terms for killifish DEGs
- 05_mouse_GSEA_mouseDEGs.R - Performs GSEA using mouse terms for mouse DEGs
- 06_mousevskillifish_compareGOterms.R - Find overlapping GO terms
- 07_mousevskillifish_analyzecommonGO.R - Plot the bubble plots for the overlapping terms
- 08_human_DESeq2_gtex.R - Get DESeq2 results for gtex samples
- 09_human_xspec.R - Compare DEGs in human + killifish dataset
- 10_gtex_forplotting.R - Plot gene boxplots by physiological system
- 11_human_GSEA_killifishDEGs.R - Performs GSEA using human terms for killifish DEGs
- 12_human_GSEA_humanDEGs.R - Performs GSEA using human terms for human DEGs
- 13_humanvskillifish_compareGOterms.R - Find overlapping GO terms
- 14_humanvskillifish_analyzecommonGOterms.R - Plot the bubble plots for the overlapping terms