# ------------------------------------------------------------------------------
# Title:          Gene Ontology GSEA using DESeq2 Killifish Data and Human Orthologs
# Author:         Emma Costa
# Date:           Code compiled on 20250808
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description:    Performs Gene Set Enrichment Analysis (GSEA) for GO terms using 
#                 DESeq2 results from killifish tissue comparisons (bin4 vs bin1),
#                 maps killifish genes to human orthologs, and runs GSEA via 
#                 clusterProfiler using human annotations.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------

# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library("DOSE")
library("clusterProfiler") 
library("org.Hs.eg.db")
library("DESeq2")

# I/O directories
indir = '/Input/'
outdir = '/Output/'

# Load DESeq2 results for killifish
load('../common_robjects/dds_DESeq2_TPS_allsamples_agesbinned_240805.bin')

# Load ortholog mapping (NCBI killifish ID to human symbol)
hSymbols = read.table(paste0(indir,"NCBI-Human-orthologs.txt"), head = T, sep = "\t")

# ------------------------------------------------------------------
# Select input tissues and iterate
# ------------------------------------------------------------------

tissue_list = names(results_list_CA1)

for(t in tissue_list){
  print(t)
  data = as.data.frame(results_list_CA1[[t]]$'4_vs_1'$resall)
  
  # ----------------------------------------------------------------
  # Generate GSEA input: signed significance score
  # ----------------------------------------------------------------
  data$mlog10QvalxFC <- (-log10(data$padj))*(data$log2FoldChange)
  toRun = subset(data, select = c(gene_symbol, mlog10QvalxFC))
  colnames(toRun) <- c('Gene', 'mlog10QvalxFC')
  
  # ----------------------------------------------------------------
  # Map killifish genes to human orthologs and get Entrez IDs
  # ----------------------------------------------------------------
  dataH = merge(hSymbols, toRun, by.x = "ncbi", by.y = "Gene") 
  entrezIds = bitr(as.character(dataH[,2]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  dataHE = merge(dataH, entrezIds, by.x = "human", by.y = "SYMBOL")
  
  # Average values across duplicates
  dataHE$mlog10QvalxFC <- as.numeric(dataHE$mlog10QvalxFC)
  unique = aggregate(dataHE[,3], list(dataHE$human), mean)
  dataHEU = merge(unique, entrezIds, by.x = "Group.1", by.y = "SYMBOL")
  colnames(dataHEU) = c("human", "mlog10QvalxFC", "entrez")
  
  geneList = dataHEU[,2]
  names(geneList) = as.character(dataHEU[,1])
  geneList = sort(geneList, decreasing = TRUE)
  
  # ----------------------------------------------------------------
  # Perform GSEA using Gene Ontology annotations
  # ----------------------------------------------------------------
  ego3 <- gseGO(
    geneList     = geneList,
    OrgDb        = org.Hs.eg.db,
    keyType      = 'SYMBOL',
    ont          = c("ALL"),
    pvalueCutoff = 1
  )
  
  # ----------------------------------------------------------------
  # Save results
  # ----------------------------------------------------------------
  write.table(ego3,
              paste0(outdir,"tables/250709_DESeq2_bin4-vs-bin1_", t ,"_humanterms_Only_GOGSEA.csv"),
              sep = ",", quote = TRUE, row.names = FALSE)
  
  # Optional: export killifish–human mapping
  # write.table(dataHE,
  #             paste0("Output/250709_DESeq2_bin4-vs-bin1_", t ,"_Only_GOGSEA_HumanName.csv"),
  #             sep = ",", quote = TRUE, row.names = FALSE)
}

# ------------------------------------------------------------------
# Session Info
# ------------------------------------------------------------------
sessionInfo()
