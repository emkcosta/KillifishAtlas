# ------------------------------------------------
# Title: Tissue-wise PCA of Killifish Atlas Samples
# Author: Emma Costa
# Date:           Code compiled on 20250811
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description: Performs PCA separately for each tissue in the Killifish Atlas
#              and visualizes PC1 vs PC2 colored by age, sex, cohort, and RNA batch.
# ------------------------------------------------

# ------------------------------------------------
# Set up
# ------------------------------------------------
rm(list = ls())  # clear environment

# Load libraries
library(DESeq2)
library(PCAtools)
library(dplyr)
library(tidyr)

# Set working directory and paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
indir = "/Input/"
outdir = "/Output/"

dir.create('/Output/Plots/')
dir.create('/Output/Tables/')


# Define colors for plots
colsex <- c('#068ec9','#ba1e2d')  # blue: male, red: female

# Colors by tissue
# Bone, Brain, Eye, Fat, Gut, Heart, Kidney, Liver, Muscle, Ovary, Skin, SpinalCord, Spleen, Testis
colTPS <- c('#d8413f', '#00a550','#b8b8c0', '#eee09b','#010101','#f0932e', '#fcd328',
            '#6cc0ee','#f4c489','#c9bad4' ,'#ab5673', '#f1a8a4','#ef9ac2','#93cca8')

# Colors by age (reversed for plotting order)
colage <- c('#395ea4','#395ea4','#395ea4','#3789c0','#3789c0','#3789c0',
            '#dfeaf0','#dfeaf0','#f3d4ac','#f3d4ac', '#f7b77b','#f7b77b','#f7b77b',
            '#a92b46','#a92b46')
colage.rev <- rev(colage)

# Load DESeq2 object (contains count matrices and sample metadata)
load('../../common_robjects/dds_TPS_allsamples_Gonadcombo_240714.bin')

# ------------------------------------------------
# Run PCA by tissue
# ------------------------------------------------
# Select tissues to include (excluding first entry, 'All')
tissue.list <- names(dds_TPS_list)[2:14]

# Loop through tissues
for (tissue in tissue.list) {
  print(tissue)
  
  # Run variance stabilizing transformation (VST)
  vsd.tissue <- vst(dds_TPS_list[[tissue]])
  
  # Perform PCA using PCAtools
  p <- pca(assay(vsd.tissue), metadata = colData(dds_TPS_list[[tissue]]))
  
  # ------------------------------------------------
  # Plot PCA colored by age
  # ------------------------------------------------
  pdf(file = paste0(outdir,"Plots/Atlas_PCA_",tissue,"_byage-sex_PC1PC2_250606.pdf"),
      width = 4, height = 4)
  print(biplot(p,
               x = 'PC1', y = 'PC2',
               lab = NA,
               colby = 'age_days',
               colkey = colage.rev,
               shape = 'sex',
               legendPosition = 'none'))
  dev.off()
  
  # ------------------------------------------------
  # Plot PCA colored by sex
  # ------------------------------------------------
  pdf(file = paste0(outdir, "Plots/Atlas_PCA_",tissue,"_bysex_PC1PC2_250606.pdf"),
      width = 4, height = 4)
  print(biplot(p,
               x = 'PC1', y = 'PC2',
               lab = NA,
               colby = 'sex',
               colkey = rev(colsex),
               legendPosition = 'none'))
  dev.off()
  
  # ------------------------------------------------
  # Plot PCA colored by cohort
  # ------------------------------------------------
  pdf(file = paste0(outdir,"Plots/Atlas_PCA_",tissue,"_PC1PC2_bycohort_250606.pdf"),
      width = 4, height = 4)
  print(biplot(p,
               x = 'PC1', y = 'PC2',
               lab = NA,
               colby = 'cohort',
               legendPosition = 'top'))
  dev.off()
  
  # ------------------------------------------------
  # Plot PCA colored by RNA batch
  # ------------------------------------------------
  pdf(file = paste0(outdir,"Plots/Atlas_PCA_",tissue,"_RNAbatch_PC1PC2_250606.pdf"),
      width = 4, height = 4)
  print(biplot(p,
               x = 'PC1', y = 'PC2',
               lab = NA,
               colby = 'RNA_batch',
               legendPosition = 'top'))
  dev.off()
}
