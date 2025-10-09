# ------------------------------------------------------------------------------
# Title:          Cross-species UP/DOWN gene boxplots by physiological system - Mouse
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

# Libraries
library(ggplot2)

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

combined_ups = read.csv(file = paste0(outdir,'tables/250712_maca_allUPsiggenes_alltissues.csv'), row.names = 1)
combined_downs = read.csv(file = paste0(outdir,'tables/250712_maca_allDOWNsiggenes_alltissues.csv'), row.names = 1)

# ------------------------------------------------------------------------------
# Tissue factor order and filtering
# ------------------------------------------------------------------------------

# Add tissue as factor with levels
my_order <- c("Brain", #nervous system
              "Brown_Fat","Gonadal_Fat", "Subcutaneous_Fat" ,"Mesenteric_Fat", #endocrine 
              "Small_Intestine", "Liver", #digestive
              "Limb_Muscle", "Bone", "Skin", #musculoskeletal + integumentary
              "Heart", "White_Blood_Cells", "Kidney","Marrow", "Spleen" #cardiovascular + immune
)

combined_ups$tissue <- factor(combined_ups$tissue, levels = my_order)
combined_downs$tissue <- factor(combined_downs$tissue, levels = my_order)

# Exclude selected tissues
combined_ups <- combined_ups %>% subset(!(tissue %in% c('Bone')))
combined_downs <- combined_downs %>% subset(!(tissue %in% c('Bone','Marrow', 'Spleen')))

# ------------------------------------------------------------------------------
# Reshape for plotting (UP genes)
# ------------------------------------------------------------------------------

plot_data_ups <- combined_ups %>%
  dplyr::select(Mouse.gene, tissue.mus, log2FoldChange.mus, log2FoldChange.kf, cor_color.kf) %>%
  pivot_longer(cols = starts_with("log2FoldChange"),
               names_to = "species",
               values_to = "log2FC") %>%
  mutate(species = recode(species,
                          log2FoldChange.mus = "Mouse",
                          log2FoldChange.kf = "Killifish"))

plot_data_ups$dot_alpha = ifelse(plot_data_ups$cor_color.kf == '#f57a7f', 0.6, 0.2)

# Recode colors to match the C4B bar plot
# If cor_color.kf =='#f57a7f', change Mouse to "#7CAE00" and Killifish to "#F8766D"
plot_data_ups$cor_color.kf[plot_data_ups$cor_color.kf == '#f57a7f' & plot_data_ups$species == 'Mouse'] <- "#7CAE00"
plot_data_ups$cor_color.kf[plot_data_ups$cor_color.kf == '#f57a7f' & plot_data_ups$species == 'Killifish'] <- "#F8766D"

# ------------------------------------------------------------------------------
# Reshape for plotting (DOWN genes)
# ------------------------------------------------------------------------------

plot_data_downs <- combined_downs %>%
  dplyr::select(Mouse.gene, tissue.mus, log2FoldChange.mus, log2FoldChange.kf, cor_color.kf) %>%
  pivot_longer(cols = starts_with("log2FoldChange"),
               names_to = "species",
               values_to = "log2FC") %>%
  mutate(species = recode(species,
                          log2FoldChange.mus = "Mouse",
                          log2FoldChange.kf = "Killifish"))

plot_data_downs$dot_alpha = ifelse(plot_data_downs$cor_color.kf == '#85d0f9', 0.6, 0.2)

# Recode colors to match the C4B bar plot
plot_data_downs$cor_color.kf[plot_data_downs$cor_color.kf == '#85d0f9' & plot_data_downs$species == 'Mouse'] <- "#7CAE00"
plot_data_downs$cor_color.kf[plot_data_downs$cor_color.kf == '#85d0f9' & plot_data_downs$species == 'Killifish'] <- "#F8766D"

# ------------------------------------------------------------------------------
# Define physiological system groupings
# ------------------------------------------------------------------------------

nerv <- c("Brain") #nervous system
endocrine <- c("Brown_Fat","Gonadal_Fat", "Subcutaneous_Fat" ,"Mesenteric_Fat") #endocrine 
digestive <- c("Small_Intestine", "Liver") #digestive
musc.integ <- c("Limb_Muscle", "Bone", "Skin") #musculoskeletal + integumentary
card.imm <- c("Heart", "White_Blood_Cells", "Kidney","Marrow", "Spleen") #cardiovascular + immune

# ------------------------------------------------------------------------------
# Plots — UP genes by system
# ------------------------------------------------------------------------------

# nervous system
plot_data_ups.filt = plot_data_ups %>% filter(tissue.mus %in% nerv)
p_1 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Mouse = '#7CAE00')) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6)) + 
  ylim(c(0,5))

pdf(paste0(outdir,'plots/250712_boxplotup_nervous_KfBrain_vs_MsBrain.pdf'), width = 1.5, height = 2)
plot(p_1)
dev.off()

# endocrine
plot_data_ups.filt = plot_data_ups %>% filter(tissue.mus %in% endocrine)
p_2 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Mouse = '#7CAE00')) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6), panel.spacing = unit(0.1, "lines")) + 
  ylim(c(0,5))

pdf(paste0(outdir,'plots/250712_boxplotup_endocrine_KfBrain_vs_MsBrain.pdf'), width = 4, height = 2)
plot(p_2)
dev.off()

# digestive
plot_data_ups.filt = plot_data_ups %>% filter(tissue.mus %in% digestive)
p_3 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Mouse = '#7CAE00')) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6), panel.spacing = unit(0.1, "lines")) + 
  ylim(c(0,5))

pdf(paste0(outdir,'plots/250712_boxplotup_digestive_KfBrain_vs_MsBrain.pdf'), width = 2, height = 2)
plot(p_3)
dev.off()

# musculoskeletal + integumentary
plot_data_ups.filt = plot_data_ups %>% filter(tissue.mus %in% musc.integ)
p_4 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Mouse = '#7CAE00')) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6), panel.spacing = unit(0.1, "lines")) +  
  ylim(c(0,5))

pdf(paste0(outdir,'plots/250712_boxplotup_muscinteg_KfBrain_vs_MsBrain.pdf'), width = 2, height = 2)
plot(p_4)
dev.off()

# cardiovascular + immune
plot_data_ups.filt = plot_data_ups %>% filter(tissue.mus %in% card.imm)
p_5 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Mouse = '#7CAE00')) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6), panel.spacing = unit(0.1, "lines")) +  
  ylim(c(0,5))

pdf(paste0(outdir,'plots/250712_boxplotup_cardio-imm_KfBrain_vs_MsBrain.pdf'), width = 5, height = 2)
plot(p_5)
dev.off()

# ------------------------------------------------------------------------------
# Plots — DOWN genes by system
# ------------------------------------------------------------------------------

# nervous system
plot_data_downs.filt = plot_data_downs %>% filter(tissue.mus %in% nerv)
p_6 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Mouse = '#7CAE00')) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6)) +  
  ylim(c(-5,0))

pdf(paste0(outdir,'plots/250712_boxplotdown_nervous_KfBrain_vs_MsBrain.pdf'), width = 1.5, height = 2)
plot(p_6)
dev.off()

# endocrine
plot_data_downs.filt = plot_data_downs %>% filter(tissue.mus %in% endocrine)
p_7 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Mouse = '#7CAE00')) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6), panel.spacing = unit(0.1, "lines")) +  
  ylim(c(-5,0))

pdf(paste0(outdir,'plots/250712_boxplotdown_endocrine_KfBrain_vs_MsBrain.pdf'), width = 4, height = 2)
plot(p_7)
dev.off()

# digestive
plot_data_downs.filt = plot_data_downs %>% filter(tissue.mus %in% digestive)
p_8 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Mouse = '#7CAE00')) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6), panel.spacing = unit(0.1, "lines")) +  
  ylim(c(-5,0))

pdf(paste0(outdir,'plots/250712_boxplotdown_digestive_KfBrain_vs_MsBrain.pdf'), width = 2, height = 2)
plot(p_8)
dev.off()

# musculoskeletal + integumentary
plot_data_downs.filt = plot_data_downs %>% filter(tissue.mus %in% musc.integ)
p_9 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Mouse = '#7CAE00')) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6), panel.spacing = unit(0.1, "lines")) +  
  ylim(c(-5,0))

pdf(paste0(outdir,'plots/250712_boxplotdown_muscinteg_KfBrain_vs_MsBrain.pdf'), width = 2, height = 2)
plot(p_9)
dev.off()

# cardiovascular + immune
plot_data_downs.filt = plot_data_downs %>% filter(tissue.mus %in% card.imm)
p_10 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  scale_colour_manual(values = c(Killifish = '#F8766D', 
                                 Mouse = '#7CAE00')) +
  theme_classic() +
  theme(legend.position="none", strip.text = element_text(size = 6),panel.spacing = unit(0.1, "lines")) + 
  ylim(c(-5,0))

pdf(paste0(outdir,'plots/250712_boxplotdown_cardioimmun_KfBrain_vs_MsBrain.pdf'), width = 3, height = 2)
plot(p_10)
dev.off()
