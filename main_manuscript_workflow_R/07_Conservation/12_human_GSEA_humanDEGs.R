# ------------------------------------------------------------------------------
# Title:          GO GSEA Analysis for Human GTEx Data using ClusterProfiler
# Author:         Emma Costa
# Date:           Code compiled on 20250808
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description:    Performs Gene Ontology enrichment using GSEA (ClusterProfiler)
#                 on DESeq2 results from GTEx tissues (age 70–79 vs 20–29).
#                 Generates ranked gene lists and computes GO enrichment scores 
#                 for each tissue using human Ensembl IDs.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------
# Load packages
# ------------------------------------------------------------------

library("clusterProfiler") 
library("org.Hs.eg.db") # orgDB for human

# ------------------------------------------------------------------
# Set working directory and input location
# ------------------------------------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
outdir = '/Output/'

# ------------------------------------------------------------------
# Load DESeq2 results
# ------------------------------------------------------------------

# Use the correct version of DESeq2 output
load(file='../common_robjects/dds_DESeq2_gtex_results_250709.bin')   # current version

tissue_list = names(results_list_gtex)

# ------------------------------------------------------------------
# Run GSEA for each tissue
# ------------------------------------------------------------------

for(t in tissue_list[1:18]){
  print(t)
  
  # Skip tissue
  if(t %in% c("brain_spinal_cord_cervical_c-1")) {
    next
  }
  
  # Extract DESeq2 result table
  data = as.data.frame(results_list_gtex[[t]]$`70-79_vs_20-29`$resall)
  if((dim(data) == 0)[1] == 'TRUE'){
    next
  }
  
  # ----------------------------------------------------------------
  # Prepare ranked gene list for GSEA
  # ----------------------------------------------------------------
  
  data$mlog10QvalxFC <- (-log10(data$padj)) * (data$log2FoldChange)
  toRun = subset(data, select = c(gene_symbol, mlog10QvalxFC))
  colnames(toRun) <- c('Gene', 'mlog10QvalxFC')
  
  # Remove Ensembl version suffix
  toRun$Gene = sub("\\..*", "", toRun$Gene)
  
  # Define named vector for GSEA input
  geneList = toRun[,2]
  names(geneList) = as.character(toRun[,1])
  geneList = sort(geneList, decreasing = TRUE)
  
  head(geneList)
  tail(geneList)
  
  # ----------------------------------------------------------------
  # Run GO Gene Set Enrichment Analysis (GSEA)
  # ----------------------------------------------------------------
  
  go_gsea <- gseGO(
    geneList     = geneList,
    OrgDb        = org.Hs.eg.db,  # Use org.Mm.eg.db for mouse
    keyType      = 'ENSEMBL',
    ont          = c("ALL"),
    pvalueCutoff = 1
  )
  
  go_gsea <- setReadable(go_gsea, 'org.Hs.eg.db', 'ENSEMBL')
  
  # ----------------------------------------------------------------
  # Save results
  # ----------------------------------------------------------------
  
  write.table(go_gsea,
              paste0(outdir,"tables/250710_DESeq2_70svs20s_", t ,"_humanterms_Only_GOGSEA.csv"),
              sep = ",", quote = TRUE, row.names = FALSE)
  
  # Optional: for up/down-regulated separately
  # write.table(go_gsea,
  #             paste0("Results/250709_DESeq2_70svs20s_", t ,"_", set, "regulated_humanterms_Only_GOGSEA.csv"),
  #             sep = ",", quote = TRUE, row.names = FALSE)
}
