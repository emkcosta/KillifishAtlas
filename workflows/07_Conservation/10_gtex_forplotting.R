# ------------------------------------------------------------------------------
# Title:          Cross-species UP/DOWN gene boxplots by physiological system - Human
# Author:         Emma Costa
# Date:           Code compiled on 20250808
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description:    Loads precomputed cross-species DE results (Mouse vs Killifish),
#                 reshapes them, and generates system-specific boxplots of log2FC
#                 for upregulated and downregulated genes across tissues.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Set up
# ------------------------------------------------------------------------------

library('dplyr')
library('readxl')
library('ggplot2')
library('tidyr')
library('biomaRt')

# Working directory (use active RStudio document path)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Output directory
outdir = '/Output/'

# Color palette reference:
# kf, ms, hum = pink, green, blue
# Example: c("#F8766D", "#7CAE00", "#00BFC4")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

combined_ups = read.csv(file = paste0(outdir,'tables/250712_gtex_allUPsiggenes_alltissues.csv'), row.names = 1)
combined_downs = read.csv(file = paste0(outdir,'tables/250712_gtex_allDOWNgenes_alltissues.csv'), row.names = 1)

# ------------------------------------------------------------------------------
# Tissue factor order and filtering
# ------------------------------------------------------------------------------

# Add tissue as factor with levels
my_order <- c("brain_cortex", #nervous system
              "adipose_visceral_omentum",  "ovary", "testis", #endocrine + reproductive
              "esophagus_muscularis", "colon_sigmoid", "colon_transverse", "small_intestine_terminal_ileum", "liver", #digestive
              "muscle_skeletal", "skin_not_sun_exposed_suprapubic", #musculoskeletal + integumentary
              "heart_left_ventricle", "whole_blood", "spleen" #cardiovascular + immune
)

combined_ups$tissue <- factor(combined_ups$tissue, levels = my_order)
combined_downs$tissue <- factor(combined_downs$tissue, levels = my_order)

# ------------------------------------------------------------------------------
# Reshape for plotting (UP genes)
# ------------------------------------------------------------------------------

# Reshape for plotting
plot_data_ups <- combined_ups %>%
  dplyr::select(Human, tissue, log2FoldChange.hum, log2FoldChange.kf, cor_color.kf) %>%
  pivot_longer(cols = starts_with("log2FoldChange"),
               names_to = "species",
               values_to = "log2FC") %>%
  mutate(species = recode(species,
                          log2FoldChange.hum = "Human",
                          log2FoldChange.kf = "Killifish"))

plot_data_ups$dot_alpha = ifelse(plot_data_ups$cor_color.kf == '#f57a7f', 0.6, 0.2)

#recode colors to match the C4B bar plot
#if cor_color.kf =='#f57a7f', change the color for mouse to "#00BFC4" and for killifish to "#F8766D"
plot_data_ups$cor_color.kf[plot_data_ups$cor_color.kf == '#f57a7f' & plot_data_ups$species == 'Human'] <- "#00BFC4"
plot_data_ups$cor_color.kf[plot_data_ups$cor_color.kf == '#f57a7f' & plot_data_ups$species == 'Killifish'] <- "#F8766D"

plot_data_ups$species = factor(plot_data_ups$species, levels = c("Killifish", "Human"))


# ------------------------------------------------------------------------------
# Reshape for plotting (DOWN genes)
# ------------------------------------------------------------------------------

plot_data_downs <- combined_downs %>%
  dplyr::select(Human, tissue, log2FoldChange.hum, log2FoldChange.kf, cor_color.kf) %>%
  pivot_longer(cols = starts_with("log2FoldChange"),
               names_to = "species",
               values_to = "log2FC") %>%
  mutate(species = recode(species,
                          log2FoldChange.hum = "Human",
                          log2FoldChange.kf = "Killifish"))


plot_data_downs$dot_alpha = ifelse(plot_data_downs$cor_color.kf == '#85d0f9', 0.6, 0.2)

#recode colors to match the C4B bar plot
plot_data_downs$cor_color.kf[plot_data_downs$cor_color.kf == '#85d0f9' & plot_data_downs$species == 'Human'] <- "#00BFC4"
plot_data_downs$cor_color.kf[plot_data_downs$cor_color.kf == '#85d0f9' & plot_data_downs$species == 'Killifish'] <- "#F8766D"


plot_data_downs$species = factor(plot_data_downs$species, levels = c("Killifish", "Human"))

# ------------------------------------------------------------------------------
# Define physiological system groupings
# ------------------------------------------------------------------------------

nerv <- c("brain_cortex") #nervous system
endocrine <- c("adipose_visceral_omentum",  "ovary", "testis") #endocrine + reproductive
digestive <- c("esophagus_muscularis", "colon_sigmoid", "colon_transverse", "small_intestine_terminal_ileum", "liver") #digestive
musc.integ <- c("muscle_skeletal", "skin_not_sun_exposed_suprapubic") #musculoskeletal + integumentary
card.imm <- c("heart_left_ventricle", "whole_blood", "spleen") #cardiovascular + immune

# ------------------------------------------------------------------------------
# Plots — UP genes by system
# ------------------------------------------------------------------------------

#nervous system
plot_data_ups.filt = plot_data_ups %>% filter(tissue %in% nerv)
p_1 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Human = "#00BFC4")) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6)) +
  ylim(c(0,5))

pdf(paste0(outdir,'plots/250712_boxplotup_nervous_KfBrain_vs_HumanBrain.pdf'), width = 1.5, height = 2)
plot(p_1)
dev.off()

#endocrine
plot_data_ups.filt = plot_data_ups %>% filter(tissue %in% endocrine)
p_2 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Human = "#00BFC4")) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6),,panel.spacing = unit(0.1, "lines")) +
  ylim(c(0,5))

pdf(paste0(outdir,'plots/250712_boxplotup_endocrine_KfBrain_vs_HumanBrain.pdf'), width = 3, height = 2)
plot(p_2)
dev.off()

#digestive
plot_data_ups.filt = plot_data_ups %>% filter(tissue %in% digestive)
p_3 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Human = "#00BFC4")) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6),panel.spacing = unit(0.1, "lines")) +
  ylim(c(0,5))

pdf(paste0(outdir,'plots/250712_boxplotup_digestive_KfBrain_vs_HumanBrain.pdf'), width = 5, height = 2)
plot(p_3)
dev.off()

#musculoskeletal + integumentary
plot_data_ups.filt = plot_data_ups %>% filter(tissue %in% musc.integ)
p_4 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Human = "#00BFC4")) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6),panel.spacing = unit(0.1, "lines")) +
  ylim(c(0,5))

pdf(paste0(outdir,'plots/250712_boxplotup_muscinteg_KfBrain_vs_HumanBrain.pdf'), width = 2, height = 2)
plot(p_4)
dev.off()

#cardiovascular + immune
plot_data_ups.filt = plot_data_ups %>% filter(tissue %in% card.imm)
p_5 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Human = "#00BFC4")) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6),panel.spacing = unit(0.1, "lines")) +
  ylim(c(0,5))

pdf(paste0(outdir,'plots/250712_boxplotup_cardioimmune_KfBrain_vs_HumanBrain.pdf'), width = 3, height = 2)
plot(p_5)
dev.off()

# ------------------------------------------------------------------------------
# Plots — DOWN genes by system
# ------------------------------------------------------------------------------

#nervous system
plot_data_downs.filt = plot_data_downs %>% filter(tissue %in% nerv)
p_6 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Human = "#00BFC4")) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6)) +
  ylim(c(-5,0))

pdf(paste0(outdir,'plots/250712_boxplotdown_nervous_KfBrain_vs_HumanBrain.pdf'), width = 1.5, height = 2)
plot(p_6)
dev.off()

#endocrine
plot_data_downs.filt = plot_data_downs %>% filter(tissue %in% endocrine)
p_7 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Human = "#00BFC4")) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6),panel.spacing = unit(0.1, "lines")) +
  ylim(c(-5,0))

pdf(paste0(outdir,'plots/250712_boxplotdown_endocrine_KfBrain_vs_HumanBrain.pdf'), width = 3, height = 2)
plot(p_7)
dev.off()

#digestive
plot_data_downs.filt = plot_data_downs %>% filter(tissue %in% digestive)
p_8 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Human = "#00BFC4")) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6),panel.spacing = unit(0.1, "lines")) +
  ylim(c(-5,0))

pdf(paste0(outdir,'plots/250712_boxplotdown_digestive_KfBrain_vs_HumanBrain.pdf'), width = 5, height = 2)
plot(p_8)
dev.off()

#musculoskeletal + integumentary
plot_data_downs.filt = plot_data_downs %>% filter(tissue %in% musc.integ)
p_9 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Human = "#00BFC4")) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6),panel.spacing = unit(0.1, "lines")) +
  ylim(c(-5,0))

pdf(paste0(outdir,'plots/250712_boxplotdown_muscinteg_KfBrain_vs_HumanBrain.pdf'), width = 2, height = 2)
plot(p_9)
dev.off()

#cardiovascular + immune
plot_data_downs.filt = plot_data_downs %>% filter(tissue %in% card.imm)
p_10 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Human = "#00BFC4")) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6),panel.spacing = unit(0.1, "lines")) +
  ylim(c(-5,0))

pdf(paste0(outdir,'plots/250712_boxplotdown_cardioimmune_KfBrain_vs_HumanBrain.pdf'), width = 3, height = 2)
plot(p_10)
dev.off()