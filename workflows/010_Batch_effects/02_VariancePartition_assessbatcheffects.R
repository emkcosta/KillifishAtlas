# ------------------------------------------------
# Title: Variance Partitioning Summary and Visualization Across Tissues
# Author: Emma Costa
# Date:           Code compiled on 20250811
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description: Loads precomputed tissue-wise variancePartition results, calculates median
#              variance explained per factor (age, sex, cohort, RNA batch), and generates
#              barplots for cohort and RNA batch contributions to gene expression variance.
# ------------------------------------------------

# ------------------------------------------------
# Set up
# ------------------------------------------------
library(DESeq2)
library(variancePartition) # v1.33.11
library(ggplot2)
library(dplyr)
library(tidyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
indir = "/Input/"
outdir = "/Output/"

# ------------------------------------------------
# Load merged TPM and metadata from prior script
# ------------------------------------------------
# Please run 000_merge_counts_TPM.R to generate this file
tpm <- read.csv(paste0(indir,'TPM_Atlas_allbatches_merged_v3.csv'),
                stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, row.names = 1)

metadata <- read.csv(file = paste0(indir, "ExperimentDesign_allbatches_combined_v7.csv"),
                     stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, row.names = 1)

# Ensure sample order matches between TPM and metadata
col_order <- rownames(metadata)
tpm <- tpm[, col_order]
table(rownames(metadata) == colnames(tpm))  # TRUE = order matches

# Load DESeq2 object for access to tissue metadata
load('../../common_robjects/dds_TPS_allsamples_Gonadcombo_240714.bin')

# Load tissue-specific variancePartition output
varPart.list <- read.csv(paste0(indir,'241121_Alltissue_variancePartition_results_Gonadcombo.csv'),
                         row.names = 1)

# ------------------------------------------------
# Compute median variance explained by factor across tissues
# ------------------------------------------------
data.summary <- data.frame(tissue = unique(varPart.list$tissue))
data.summary$age.median <- NA
data.summary$sex.median <- NA
data.summary$sex.age.median <- NA
data.summary$cohort.median <- NA
data.summary$RNAbatch.median <- NA

# Loop through tissues to compute medians
for(i in 1:nrow(data.summary)) {
  temp <- data.summary[i, ]
  t <- temp$tissue
  df <- varPart.list %>% filter(tissue == t)
  data.summary[i, ]$age.median      <- median(df$age_bin) * 100
  data.summary[i, ]$sex.median      <- median(df$sex) * 100
  data.summary[i, ]$sex.age.median  <- median(df$sex.age_bin) * 100
  data.summary[i, ]$cohort.median   <- median(df$cohort) * 100
  data.summary[i, ]$RNAbatch.median <- median(df$RNA_batch) * 100
}

# ------------------------------------------------
# Barplot: Median % variance explained by COHORT
# ------------------------------------------------

# Tissue colors for cohort plot
plotcols <- c('#b8b8c0','#010101', '#6cc0ee', '#ab5673', '#eee09b', '#d8413f',
              '#fcd328', '#f1a8a4', '#f4c489', '#f0932e', '#ef9ac2', '#7962A3', '#00a550')

# Reorder tissues by cohort contribution
data.summary <- data.summary[order(-data.summary$cohort.median), ]
data.summary$tissue <- factor(data.summary$tissue, levels = data.summary$tissue)

# Plot
pdf(file = paste0(outdir, "Plots/250529_Alltissue_variancePartition_Cohort-Median_barPlot_zeroto20.pdf"),
    width = 2, height = 2)

ggplot(data.summary, aes(x = tissue, y = cohort.median)) +
  geom_bar(stat = 'identity', aes(fill = tissue)) +
  scale_fill_manual(values = plotcols) +
  theme_classic() +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 20)) +
  ggtitle('Cohort - % Variance Explained') +
  ylab('Median % Variance')

dev.off()

# ------------------------------------------------
# Barplot: Median % variance explained by RNA BATCH
# ------------------------------------------------

# Tissue colors for RNA batch plot
plotcols2 <- c('#f4c489', '#d8413f', '#ab5673', '#f0932e', '#6cc0ee', '#7962A3',
               '#010101', '#eee09b', '#fcd328', '#f1a8a4', '#ef9ac2', '#00a550', '#b8b8c0')

# Reorder tissues by RNA batch contribution
data.summary$tissue <- droplevels(data.summary$tissue)
data.summary <- data.summary[order(-data.summary$RNAbatch.median), ]
data.summary$tissue <- factor(data.summary$tissue, levels = data.summary$tissue)

# Plot
pdf(file = paste0(outdir, "Plots/250529_Alltissue_variancePartition_RNAbatch-Median_barPlot_zeroto20.pdf"),
    width = 2, height = 2)

ggplot(data.summary, aes(x = tissue, y = RNAbatch.median)) +
  geom_bar(stat = 'identity', aes(fill = tissue)) +
  scale_fill_manual(values = plotcols2) +
  theme_classic() +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 20)) +
  ggtitle('RNA Batch - % Variance Explained') +
  ylab('Median % Variance')

dev.off()

# ------------------------------------------------
# Save summary statistics table
# ------------------------------------------------
write.csv(data.summary, file = paste0(outdir, "Tables/250606_Alltissue_variancePartition_medians_gonadcombo.csv"))
