# ------------------------------------------------------------------------------
# Title:          CASA gene set barplots (Killifish vs Mouse, Brain)
# Author:         Emma Costa
# Date:           Code compiled on 20250808
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description:    Loads cross-species unfiltered DE results, intersects with the
#                 CAS gene set, and plots mini barplots (log2FC) for Brain.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Set up
# ------------------------------------------------------------------------------

# Libraries
library('readxl')
library('dplyr')
library('tidyr')
library('ggplot2')

# Working directory (use active RStudio document path)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Output directory
outdir = '/Output/'

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# Cross-species, unfiltered DE results across tissues (from 01_mouse_xspec.R)
combined_unfiltered <- read.csv(paste0(outdir,'tables/250712_maca_allgenes_alltissues.csv'))

# CAS gene set (Ohahn)
cas_genes = read_xlsx(paste0("CAS_geneset_Ohahn.xlsx"))

# ------------------------------------------------------------------------------
# Filter to CAS genes in Brain and prep plotting dataframe
# ------------------------------------------------------------------------------

geneset = cas_genes$`Gene Symbol`

getgenes = combined_unfiltered %>% 
  filter(Mouse.gene %in% geneset & tissue.mus == 'Brain' & padj.kf < 0.05 & padj.mus < 0.05)

# Combine mouse and killifish gene IDs for x-axis labels
getgenes$gene_combo = paste0(getgenes$Mouse.gene, "_", getgenes$gene_symbol)

# Pivot to long format for plotting
plot_data_up <- getgenes %>%
  dplyr::select(gene_combo, log2FoldChange.mus, log2FoldChange.kf, padj.mus, padj.kf) %>%
  pivot_longer(
    cols = starts_with("log2FoldChange"),
    names_to = "species",
    values_to = "log2FC"
  )

# species labels
plot_data_up$species <- recode(plot_data_up$species,
                               log2FoldChange.mus = "Mouse",
                               log2FoldChange.kf  = "Killifish")


write.csv(plot_data_up,paste0(outdir, 'tables/250915_casgenes_overlap.csv'))
# ------------------------------------------------------------------------------
# Plot: mini barplots of log2FC per gene (Mouse vs Killifish)
# ------------------------------------------------------------------------------

p3 <- ggplot(plot_data_up, aes(x = gene_combo, y = log2FC, fill = species)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  theme_minimal() +
  labs(
    x = "Gene",
    y = "log2 Fold Change"
  ) +
  scale_fill_manual(values = c("Mouse" = "#f57a7f", "Killifish" = "#85d0f9")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    legend.position = 'none'
  )



# Save PDF
pdf(paste0(outdir,'plots/250712_CASAgingScore-sigOnly_barplot_KfBrain_vs_MsBrain.pdf'), width = 6, height = 2)
plot(p3)
dev.off()
