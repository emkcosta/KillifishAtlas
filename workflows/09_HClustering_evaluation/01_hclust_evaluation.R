# ------------------------------------------------
# Title: Silhouette analysis of gene expression trajectories in aging tissues
# Author: Emma Costa
# Date:           Code compiled on 20250811
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description: Performs hierarchical clustering and silhouette analysis 
#              on LOESS-smoothed gene expression trajectories across age. 
#              Used to determine optimal cluster number in aging transcriptome data.
# ------------------------------------------------

# ------------------------------------------------
# Load libraries
# ------------------------------------------------
library(tidyverse)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(stats)
library(cluster)

# ------------------------------------------------
# Set working directory and define paths
# ------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
indir = "/Input/"
outdir = "/Output/"

dir.create('./Output/plots/')
dir.create('./Output/tables/')

# ------------------------------------------------
# Load data
# ------------------------------------------------
# Load DESeq2 object and shared gene list
load('../common_robjects/dds_TPS_allsamples_Gonadcombo_240714.bin')
gene.set <- read.csv(file = paste0(indir,'240905_geneexpressedinalltissues.csv'), row.names = 1)

# Select tissue for analysis
tissue <- c('Brain')  # modify as needed

# Extract tissue-specific data
dds <- dds_TPS_list[[tissue]]
sampleTable <- as.data.frame(dds@colData)
countTable <- counts(dds, normalized = TRUE)  # normalized count matrix

# ------------------------------------------------
# Prepare data for clustering
# ------------------------------------------------
# Remove middle timepoint and order samples by age
sampleTable_filtered <- sampleTable %>% filter(!age_days %in% c(102,103))
sampleTable.ordered <- sampleTable_filtered %>% arrange(age_days)
ordered_samples <- sampleTable.ordered$lib

# Subset and order count matrix by selected genes and samples
norm_counts <- countTable[, ordered_samples]
norm_counts <- norm_counts[gene.set$GeneID,]  # keep shared genes only

# Z-score normalize expression per gene
z_scaled_counts <- t(scale(t(norm_counts)))

# ------------------------------------------------
# LOESS fitting of gene expression across age
# ------------------------------------------------
ages <- as.numeric(as.character(unique(sampleTable.ordered$age_days)))
loess_fits <- data.frame(matrix(nrow = nrow(z_scaled_counts), ncol = length(ages)))
rownames(loess_fits) <- rownames(z_scaled_counts)
colnames(loess_fits) <- ages

for (gene in rownames(z_scaled_counts)) {
  gene_data <- data.frame(Age = sampleTable.ordered$age_days, Expression = z_scaled_counts[gene, ])
  gene_data$Age <- as.numeric(as.character(gene_data$Age))
  loess_fit <- loess(Expression ~ Age, data = gene_data)
  loess_fits[gene, ] <- predict(loess_fit, ages)
}

# ------------------------------------------------
# Hierarchical clustering and silhouette analysis
# ------------------------------------------------
# Compute Euclidean distance matrix
dist_matrix <- dist(loess_fits, method = "euclidean")

# Perform hierarchical clustering using complete linkage
hc <- hclust(dist_matrix, method = "complete")

# Loop over selected cluster numbers for silhouette plots
for (num_clusters in c(5,10)) {
  
  # Cut dendrogram into k clusters
  clusters <- cutree(hc, k = num_clusters)
  loess_fits$Cluster <- as.factor(clusters)
  
  # Calculate silhouette widths
  sil <- silhouette(clusters, dist = dist_matrix)
  
  # Save silhouette plot
  pdf(file = paste0(outdir, "plots/hclust_",tissue,"silhouette_k_", num_clusters, "_250611.pdf"), width = 4, height = 6)
  print(plot(sil, border = NA, main = paste("Silhouette Plot for", num_clusters, "Clusters")))
  dev.off()
}

# ------------------------------------------------
# Function: Evaluate silhouette widths over k = 2:20
# ------------------------------------------------
evaluate_silhouette <- function(data, dist_method = "euclidean", hclust_method = "complete", k_range = 2:20) {
  
  # Compute distance matrix
  d <- dist(data, method = dist_method)
  avg_sil_widths <- numeric(length(k_range))
  
  for (i in seq_along(k_range)) {
    k <- k_range[i]
    hc <- hclust(d, method = hclust_method)
    clusters <- cutree(hc, k)
    sil <- silhouette(clusters, dist = d)
    avg_sil_widths[i] <- mean(sil[, 3])  # average silhouette width
  }
  
  # Plot average silhouette width by k
  pdf(file = paste0(outdir, "plots/hclust_",tissue,"silhouette_k_2to20_250611.pdf"), width = 4, height = 4)
  print(plot(k_range, avg_sil_widths, type = "b", pch = 19, frame = FALSE,
             xlab = "Number of clusters (k)", ylab = "Average silhouette width",
             main = "Silhouette Analysis for Hierarchical Clustering"))
  dev.off()
  
  return(data.frame(k = k_range, avg_sil_width = avg_sil_widths))
}

# Run silhouette evaluation over 2–20 clusters
loess_fits$Cluster <- NULL
evaluate_silhouette(loess_fits)
