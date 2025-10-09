# ------------------------------------------------------------------------------
# Title:          Generate TPM signature matrices from single-cell kidney Seurat object
# Author:         Emma Costa
# Date:           Code compiled on 20250811
# Description:    Loads a Seurat object of killifish kidney single-cell data,
#                 calculates TPM values using gene lengths, and generates
#                 average expression profiles by cell type for all, male-only,
#                 and female-only subsets for deconvolution purposes.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

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

# load gene lengths
gene.lengths = read.csv(paste0(indir, 'GeneLengths_Atlas_Plate02_Lane03-4_BatchF004_240417.csv'))
gene.lengths = gene.lengths[,1:2]
gene.lengths$gene.name = row.names(gene.lengths)
gene.lengths$A10 = NULL
colnames(gene.lengths)[1] = 'length'


# ------------------------------------------------------------------------------
# Function to generate TPM matrix and average by cell type
# ------------------------------------------------------------------------------

get_avg_tpm <- function(seurat_obj, suffix) {
  gene_set <- rownames(seurat_obj)
  common.genes = intersect(gene_set, row.names(gene.lengths))
  gene.lengths.sub = gene.lengths[common.genes, ]
  
  raw_counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  counts_sub <- raw_counts[common.genes, ]
  lengths_vector <- gene.lengths.sub$length[match(rownames(counts_sub), gene.lengths.sub$gene.name)]
  
  rate <- sweep(as.matrix(counts_sub), 1, lengths_vector, FUN = "/")
  tpm <- sweep(rate, 2, colSums(rate), FUN = "/") * 1e6
  
  meta <- seurat_obj@meta.data
  meta$numID <- 1:nrow(meta)
  meta$cellID <- paste0(gsub("_","-", meta$Annotation_v1), "_", meta$numID)
  
  new_names <- meta[colnames(tpm), "cellID"]
  colnames(tpm) <- new_names
  
  cell_types <- sub("_.*", "", colnames(tpm))
  cell_type_groups <- split(seq_along(cell_types), cell_types)
  
  avg_tpm <- sapply(cell_type_groups, function(cols) {
    Matrix::rowMeans(tpm[, cols, drop = FALSE])
  })
  
  avg_tpm_df <- as.data.frame(avg_tpm)
  write.csv(avg_tpm_df, file = paste0(outdir, suffix))
}

# ------------------------------------------------------------------------------
# Signature matrix: All cells
# ------------------------------------------------------------------------------

get_avg_tpm(onlykidney, '250612_scKidney_TPM_avg_granulator.csv')

# ------------------------------------------------------------------------------
# Signature matrix: Male only
# ------------------------------------------------------------------------------

male_kidney <- subset(onlykidney, subset = Sex == "M")
get_avg_tpm(male_kidney, '250612_scKidney_maleOnly_TPM_avg_granulator.csv')

# ------------------------------------------------------------------------------
# Signature matrix: Female only
# ------------------------------------------------------------------------------

female_kidney <- subset(onlykidney, subset = Sex == "F")
get_avg_tpm(female_kidney, '250612_scKidney_femaleOnly_TPM_avg_granulator.csv')
