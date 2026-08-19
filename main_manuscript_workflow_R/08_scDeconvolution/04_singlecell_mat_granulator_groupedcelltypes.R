# ------------------------------------------------------------------------------
# Title:          Generate grouped TPM signature matrices from scKidney Seurat object
# Author:         Emma Costa
# Date:           Code compiled on 20250811
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description:    This script processes single-cell RNA-seq data from kidney tissue
#                 to generate TPM-normalized average expression matrices grouped
#                 by broader cell type categories (male and female separately).
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

# Load libraries
library(Seurat)

# Note: For memory issues, see: https://stackoverflow.com/questions/51248293/error-vector-memory-exhausted-limit-reached-r-3-5-0-macos

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
indir = '/Input/'
outdir = '/Output/'

#load object containing single cell data 
load('../common_robjects/kidneyOnly_Teefy_SeurObj_Fordeconvolution.Rdata')
onlykidney <- SetIdent(onlykidney, value = 'Annotation_v1')
#DimPlot(onlykidney)

# ------------------------------------------------------------------------------
# Load Seurat object and gene length data
# ------------------------------------------------------------------------------

# Load gene lengths
gene.lengths = read.csv(paste0(indir, 'GeneLengths_Atlas_Plate02_Lane03-4_BatchF004_240417.csv'))
gene.lengths = gene.lengths[,1:2]
gene.lengths$gene.name = row.names(gene.lengths)
gene.lengths$A10 = NULL
colnames(gene.lengths)[1] = 'length'

# ------------------------------------------------------------------------------
# Generate male-only signature matrix (TPM)
# ------------------------------------------------------------------------------

male_kidney <- subset(onlykidney, subset = Sex == "M")
gene_set <- rownames(male_kidney)
common.genes = intersect(gene_set, row.names(gene.lengths))
gene.lengths = gene.lengths[common.genes,]
raw_counts <- GetAssayData(male_kidney, assay = "RNA", slot = "counts")
counts_sub <- raw_counts[common.genes, ]
lengths_vector <- gene.lengths$length[match(rownames(counts_sub), gene.lengths$gene.name)]
rate <- sweep(as.matrix(counts_sub), 1, lengths_vector, FUN = "/")
tpm <- sweep(rate, 2, colSums(rate), FUN = "/") * 1e6

# Rename columns
meta <- male_kidney@meta.data
meta$numID <- 1:dim(meta)[1]
meta$cellID <- paste0(gsub("_","-",meta$Annotation_v1), "_",meta$numID)
new_names <- meta[colnames(tpm), "cellID"]
colnames(tpm) <- new_names

# Define broader cell type groups
lymph.prog <- c('B-Cell-Progenitors','Lymphoid-progenitors','NK-T-progenitor-cells')
other.prog <- c('Erythrocyte-Progenitors','HSPCs','Multipotent-progenitors')
mye.prog   <- c('Myeloid-progenitors','Neutrophil-Progenitors')
lymphoid   <- c('B-cells','NK-T-cells')
myeloid    <- c('Macrophages', 'Mast-cells','Neutrophils','Thrombocytes')
non.immune <- c('Endothelial','Fibroblasts','Kidney-distal-tubule','Kidney-prox-tubule')
other      <- c('Erythrocytes')

cell_type_groups <- list(
  lymph_progenitors = lymph.prog,
  other_progenitors = other.prog,
  myeloid_progenitors = mye.prog,
  lymphoid_cells = lymphoid,
  myeloid_cells = myeloid,
  nonimmune_cells = non.immune,
  other_cells = other
)

cell_prefixes <- sub("_.*", "", colnames(tpm))
names(cell_prefixes) <- colnames(tpm)

avg_tpm_grouped <- sapply(cell_type_groups, function(fine_types) {
  cell_cols <- names(cell_prefixes)[cell_prefixes %in% fine_types]
  Matrix::rowMeans(tpm[ , cell_cols, drop = FALSE])
})

avg_tpm_grouped_df <- as.data.frame(avg_tpm_grouped)
write.csv(avg_tpm_grouped_df, file = paste0(outdir,'250624_scKidney_maleOnly_TPM_avg_granulator-grouped.csv'))

# ------------------------------------------------------------------------------
# Generate female-only signature matrix (TPM)
# ------------------------------------------------------------------------------

female_kidney <- subset(onlykidney, subset = Sex == "F")
gene_set <- rownames(female_kidney)
common.genes = intersect(gene_set, row.names(gene.lengths))
gene.lengths = gene.lengths[common.genes,]
raw_counts <- GetAssayData(female_kidney, assay = "RNA", slot = "counts")
counts_sub <- raw_counts[common.genes, ]
lengths_vector <- gene.lengths$length[match(rownames(counts_sub), gene.lengths$gene.name)]
rate <- sweep(as.matrix(counts_sub), 1, lengths_vector, FUN = "/")
tpm <- sweep(rate, 2, colSums(rate), FUN = "/") * 1e6

# Rename columns
meta <- female_kidney@meta.data
meta$numID <- 1:dim(meta)[1]
meta$cellID <- paste0(gsub("_","-",meta$Annotation_v1), "_",meta$numID)
new_names <- meta[colnames(tpm), "cellID"]
colnames(tpm) <- new_names

# Use same cell_type_groups defined above
cell_prefixes <- sub("_.*", "", colnames(tpm))
names(cell_prefixes) <- colnames(tpm)

avg_tpm_grouped <- sapply(cell_type_groups, function(fine_types) {
  cell_cols <- names(cell_prefixes)[cell_prefixes %in% fine_types]
  Matrix::rowMeans(tpm[ , cell_cols, drop = FALSE])
})

avg_tpm_grouped_df <- as.data.frame(avg_tpm_grouped)
write.csv(avg_tpm_grouped_df, file = paste0(outdir,'250624_scKidney_femaleOnly_TPM_avg_granulator-grouped.csv'))
