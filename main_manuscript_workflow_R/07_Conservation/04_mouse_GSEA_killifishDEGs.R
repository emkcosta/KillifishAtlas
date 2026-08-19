# ------------------------------------------------------------------------------
# Title:          Functional enrichment (GSEA) using GO terms (mouse ontology)
# Author:         
# Date:           Code compiled on 20250808
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description:    For each killifish tissue, computes a ranked gene list
#                 (-log10 adj.P × log2FC), maps to mouse symbols/Entrez,
#                 runs gseGO across ALL ontologies, and writes results.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Set up
# ------------------------------------------------------------------------------

# Working directory (uses active RStudio document path)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
library("DOSE")
library("clusterProfiler")
library("org.Mm.eg.db")
library("DESeq2")

# I/O directories
indir  = '/Input/'
outdir = '/Output/'

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# Killifish DESeq2 objects (contains results_list_CA1, etc.)
load('../../common_robjects/dds_DESeq2_TPS_allsamples_agesbinned_240805.bin')

# Orthologs table (killifish NCBI → mouse symbol)
# Note: Using mouse orthologs as the target space for GSEA.
ortho       <- read.csv(paste0(indir,'nfur-ncbi_orthologs_20170302_pps.csv'))
ortho.mus   <- ortho[,c(1,5)]
colnames(ortho.mus) <- c('ncbi', 'mouse')
ortho.mus   <- ortho.mus[-c(which(ortho.mus$mouse == '')),]

mSymbols <- ortho.mus   # alias used below

# ------------------------------------------------------------------------------
# Select input (tissue loop)
# ------------------------------------------------------------------------------

# Full tissue list from loaded DESeq2 results
tissue_list <- names(results_list_CA1)
t <- tissue_list[1] # debugging

for(t in tissue_list){
  print(t)
  
  # DE results (4_vs_1 age-bin comparison) for current tissue
  data <- as.data.frame(results_list_CA1[[t]]$'4_vs_1'$resall)
  
  # ----------------------------------------------------------------------------
  # Build ranked list for GSEA: score = (-log10(padj)) * log2FoldChange
  # ----------------------------------------------------------------------------
  data$mlog10QvalxFC <- (-log10(data$padj)) * (data$log2FoldChange)
  toRun <- subset(data, select = c(gene_symbol, mlog10QvalxFC))
  colnames(toRun) <- c('Gene', 'mlog10QvalxFC')
  
  # ----------------------------------------------------------------------------
  # Map killifish NCBI IDs to mouse symbols, then to Entrez IDs
  # ----------------------------------------------------------------------------
  # Merge Kf NCBI IDs with scores via gene symbol key
  dataH <- merge(mSymbols, toRun, by.x = "ncbi", by.y = "Gene")
  
  # SYMBOL → ENTREZID using mouse org DB
  entrezIds <- bitr(as.character(dataH[,2]), fromType="SYMBOL",
                    toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # Note: some IDs may fail to map (expected)
  
  # Attach the mapping (mouse symbols) to scored genes
  dataHE <- merge(dataH, entrezIds, by.x = "mouse", by.y = "SYMBOL")
  head(dataHE)
  
  # Handle paralogs: average scores by mouse symbol
  dataHE$mlog10QvalxFC <- as.numeric(dataHE$mlog10QvalxFC)
  unique <- aggregate(dataHE[,3], list(dataHE$mouse), mean)
  dataHEU <- merge(unique, entrezIds, by.x = "Group.1", by.y = "SYMBOL")
  colnames(dataHEU) <- c("mouse", "mlog10QvalxFC", "entrez")
  head(dataHEU)
  
  # Gene list for GSEA (named numeric vector)
  geneList <- dataHEU[,2]
  names(geneList) <- as.character(dataHEU[,1])
  
  # IMPORTANT: sort decreasing for GSEA
  geneList <- sort(geneList, decreasing = TRUE)
  
  head(geneList)
  tail(geneList)
  
  # ----------------------------------------------------------------------------
  # GSEA: Gene Ontology (ALL: BP/CC/MP in this call is ALL GO ontologies)
  # ----------------------------------------------------------------------------
  ego3 <- gseGO(geneList     = geneList,
                OrgDb        = org.Mm.eg.db,
                keyType      = 'SYMBOL',
                ont          = c("ALL"),
                pvalueCutoff = 1)
  
  # Write GSEA results
  write.table(ego3,
              paste0(outdir,"tables/250710_DESeq2_bin4-vs-bin1_", t ,"_mouseterms_Only_GOGSEA.csv"),
              sep = ",", quote = TRUE, row.names = FALSE)
  
}

# ------------------------------------------------------------------------------
# Session info
# ------------------------------------------------------------------------------

sessionInfo()
