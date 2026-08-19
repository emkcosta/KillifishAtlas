# ------------------------------------------------------------------------------
# Title:          Plot distribution of gene-age Spearman correlations across tissues
# Author:         Emma Costa
# Date:           Code finalized on 20250807
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description:    This script loads precomputed age correlation results from bulk
#                 RNA-seq data, aggregates them across tissues, and plots the 
#                 distribution of Spearman correlations per tissue.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Set up
# ------------------------------------------------------------------------------

# Load libraries
library(viridis)
library(dplyr)
library(ggplot2)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load R object with correlation results
load('../common_robjects/dds_BulkSeq_Aging_CorrelationResults_allsamples_Gonadcombo_TPMcutoff_240714.bin') #moved from another directory

# Set output directory
outdir = '/Output/'
dir.create('/Output/Plots/')
dir.create('/Output/Tables/')


# Define tissue-specific color palette
# Tissues: Bone, Brain, Eye, Fat, Gonad, Gut, Heart, Kidney, Liver, Muscle, Skin, SpinalCord, Spleen
colTPS <- c('#d8413f', '#00a550', '#b8b8c0', '#eee09b', '#7961a2', '#010101',
            '#f0932e', '#fcd328', '#6cc0ee', '#f4c489', '#ab5673', '#f1a8a4', '#ef9ac2')

# ------------------------------------------------------------------------------
# Aggregate Spearman correlation results across tissues
# ------------------------------------------------------------------------------

# Extract tissue names
tissue.list <- names(BulkSeq_Agecorrelation_results_list)

# Initialize combined results using the first tissue
all.corres <- BulkSeq_Agecorrelation_results_list[[tissue.list[1]]]$resall
all.corres$tissue <- tissue.list[1]

# Loop through remaining tissues and append
for(t in tissue.list[2:13]){
  tissue.res <- BulkSeq_Agecorrelation_results_list[[t]]$resall 
  tissue.res$tissue <- t
  all.corres <- rbind(tissue.res, all.corres)
}

# Add absolute correlation for downstream filtering
all.corres$abs_spear <- abs(all.corres$cor_spear)

# Define alpha for plotting: more opaque if abs(cor) > 0.5
all.corres$alpha <- ifelse(all.corres$abs_spear > 0.5, 1.0, 0.05)

# ------------------------------------------------------------------------------
# Reorder tissues by number of highly age-correlated genes
# ------------------------------------------------------------------------------

# Count genes with abs(cor_spear) > 0.5
frequency_table <- all.corres %>%
  filter(alpha == 1.0) %>%
  count(tissue, name = "count") %>%
  arrange(count)

# Reorder tissue factor levels based on count
all.corres <- all.corres %>%
  mutate(tissue = factor(tissue, levels = frequency_table$tissue))

# Reordered color palette (match new tissue order)
# Reordered tissues: Bone, Gonad, Liver, Gut, Spleen, Kidney, Heart, Fat, Brain, SpinalCord, Skin, Muscle, Eye
colTPS.reorder <- c('#d8413f', '#7961a2', '#6cc0ee', '#010101', '#ef9ac2', '#fcd328',
                    '#f0932e', '#eee09b', '#00a550', '#f1a8a4', '#ab5673', '#f4c489', '#b8b8c0')

# ------------------------------------------------------------------------------
# Plot: Distribution of Spearman correlations by tissue
# ------------------------------------------------------------------------------

pdf(paste0(outdir, 'Plots/alltissue_sexcombined_cor-spear_scatterplot_bytissue_orderbycorrnum.pdf'),
    width = 2, height = 3)

all.corres %>%
  ggplot(aes(group = tissue, y = tissue, x = cor_spear)) +
  geom_point(position = position_jitter(height = 0.2),
             size = 0.01,
             aes(color = tissue, alpha = alpha)) +
  scale_x_continuous(limits = c(-1, 1)) +
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = 'longdash', col = 'red') +
  scale_color_manual(values = colTPS.reorder) +
  geom_boxplot(outlier.shape = NA, aes(alpha = 0.1)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

dev.off()
