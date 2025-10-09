# ------------------------------------------------------------------------------
# Title:          Functional enrichment (GSEA) with clusterProfiler (mouse)
# Author:         
# Date:           Code compiled on 20250808
# Description:    Refactors DESeq2 results to a ranked list (log2FC × -log10 FDR)
#                 and runs GO GSEA (ALL ontologies) per tissue. 
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

## Installation notes (kept as comments from original)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "DOSE"))

# Load packages
library("clusterProfiler") 
library("org.Mm.eg.db") # orgDB for mouse
library('dplyr')
library('tidyr')

# Working directory (uses active RStudio document path)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Input directory
indir = '/Input/'
outdir = '/Output/'

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# DE analysis object (all ages compared to 3 months)
load(paste0(indir,'dds_MACA.bin')) # DE analysis all ages compared to 3mo

# ------------------------------------------------------------------------------
# Refactor DESeq2 results → tidy/wide + ranking metric
# ------------------------------------------------------------------------------

# Keep only the 27_vs_3 comparison
dds.filt <- dds[, grepl("27_vs_3", colnames(dds))] # only select the 27 months vs 3 months condition
dds.filt$ensembl_gene_id <- rownames(dds.filt)

# Convenience alias
df = dds.filt

# Long format (keep Ensembl ID)
df_long <- df %>%
  pivot_longer(
    cols = -c(ensembl_gene_id),
    names_to = "condition",
    values_to = "value"
  )

# Parse tissue, comparison, metric from column names
df_long <- df_long %>%
  extract(condition, into = c("tissue", "comparison", "metric"),
          regex = "(.+?)_(\\d+_vs_\\d+)_(log2FoldChange|padj)") 

# Wide format: side-by-side log2FC and padj
df_wide <- df_long %>%
  pivot_wider(
    names_from = metric,
    values_from = value
  )

# Ranking score (pi-value style): log2FC × -log10(padj)
# Reference: https://pmc.ncbi.nlm.nih.gov/articles/PMC3957066/
df_wide$mlog10Qvalxlog2FC = df_wide$log2FoldChange * (-log10(df_wide$padj))

# ------------------------------------------------------------------------------
# Per-tissue GSEA
# ------------------------------------------------------------------------------

tissue.list = unique(df_wide$tissue)
t = tissue.list[1] # debugging

for(t in tissue.list){
  print(t)
  
  # Ranked list input: Ensembl ID + quantitative score
  # (Must be a 2-column table: ID, ranking metric)
  data = df_wide %>% filter(tissue == t) # *** mouse gene ids
  head(data)
  
  toRun = as.data.frame(data[,c(1,6)])
  
  # Build named numeric vector for GSEA
  geneList = toRun[,2] 
  names(geneList) = as.character(toRun[,1]) # Ensembl IDs as names
  
  # IMPORTANT: sort decreasing for GSEA
  geneList = sort(geneList, decreasing = TRUE)
  
  head(geneList)
  tail(geneList)
  
  # ----------------------------------------------------------------------------
  # GO GSEA (ALL ontologies)
  # Switch OrgDb + keyType as needed (mouse vs human)
  # ----------------------------------------------------------------------------
  go_gsea <- gseGO(geneList     = geneList,
                   OrgDb        = org.Mm.eg.db, # Use org.Hs.eg.db for human; and org.Mm.eg.db for mouse
                   keyType      = 'ENSEMBL',
                   ont          = c("ALL"),
                   pvalueCutoff = 1)
  go_gsea <- setReadable(go_gsea, 'org.Mm.eg.db', 'ENSEMBL') # Use org.Hs.eg.db for human; and org.Mm.eg.db for mouse
  
  # Save results (mouse example filename)
  write.table(go_gsea, paste0(outdir,"tables/mouse_GOGSEA_",t,".csv"),
              sep = ",", quote = TRUE, row.names = FALSE) 
}

# ------------------------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------------------------

rm(list=ls())
