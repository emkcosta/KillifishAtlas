# ------------------------------------------------------------------------------
# Title:          Visual summaries of overlapping GO BP terms (killifish vs mouse)
# Author:         
# Date:           Code compiled on 20250808
# Description:    Loads per-tissue overlapping GO BP terms, computes similarity-
#                 based reductions (rrvgo), and generates heatmap, scatter,
#                 treemap, and bubble plots. Also produces summary barplots of
#                 overlap counts across tissues.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

library('dplyr')
library('tidyr')
library('ggplot2')
library('GOSemSim')
library('rrvgo')

# Reference:
# https://yulab-smu.top/biomedical-knowledge-mining-book/useful-utilities.html

# Working directory (uses active RStudio document path)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Output directory
outdir = '/Output/'

# Obtain GO structure (optional; kept commented per original)
# mmGO <- godata('org.Mm.eg.db', ont="BP")

# ------------------------------------------------------------------------------
# Single-tissue selection
# ------------------------------------------------------------------------------

## Select the tissue (uncomment exactly one pair if switching)
# # Brain
# kf.tissue = 'Brain'
# mus.tissue = 'Brain'
# 
# # Liver
kf.tissue = 'Liver'
mus.tissue = 'Liver'
# 
# # Heart
# kf.tissue = 'Heart'
# mus.tissue = 'Heart'

# Fat
# kf.tissue = 'Fat'
# mus.tissue = 'Mesenteric_Fat'

# ------------------------------------------------------------------------------
# Load overlapping BP terms for the selected tissue pair
# ------------------------------------------------------------------------------

commonBP <- read.csv(paste0(outdir, "250710_BPcommonGOIDs_",kf.tissue,"-vs-",mus.tissue,"_atlas-vs-maca.csv"),
                     row.names = 1)

# Scores: -log10(padj) × NES (species-specific) and a summed score
commonBP$mlog10QvalxNES.mus = -log10(commonBP$p.adjust.mus) * commonBP$NES.mus
commonBP$mlog10QvalxNES.kf  = -log10(commonBP$p.adjust.kf)  * commonBP$NES.kf
commonBP$summed_score        = commonBP$mlog10QvalxNES.mus + commonBP$mlog10QvalxNES.kf

# Minimal table, significance filter, and sorting
commonBP.simple <- commonBP[,c(1:3, 6, 8, 16, 18, 26)]
commonBP.simple <- commonBP.simple %>% filter(p.adjust.mus < 0.1 & p.adjust.kf < 0.1)
commonBP.simple <- commonBP.simple %>% arrange(desc(summed_score))

# Notes on NES distributions by tissue (kept as comments from original)
# brain: only +ve NES remain 
# liver: mostly +ve NES
# fat: pos and neg, good distribution
# heart: mostly pos, some negative

# ------------------------------------------------------------------------------
# Similarity reduction of GO BP terms (rrvgo)
# ------------------------------------------------------------------------------

# Similarity matrix among selected GO terms
simMatrix <- calculateSimMatrix(commonBP.simple$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")

# Reduce by semantic similarity (threshold agnostic due to prior padj < 0.1)
scores <- setNames(commonBP.simple$p.adjust.kf, commonBP.simple$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.1,
                                orgdb="org.Mm.eg.db")

# Heatmap of reduced terms
pdf(paste0(outdir, kf.tissue,'_vs_',mus.tissue,'_GSEA_GOBP_overlappingGO-heatmap_250719.pdf'),
    width = 20, height = 8)
print(heatmapPlot(simMatrix,
                  reducedTerms,
                  annotateParent=TRUE,
                  annotationLabel="parentTerm",
                  fontsize=6))
dev.off()

# Scatter of reduced terms
pdf(paste0(outdir, kf.tissue,'_vs_',mus.tissue,'_GSEA_GOBP_overlappingGO-scatter_250719.pdf'),
    width = 20, height = 8)
print(scatterPlot(simMatrix, reducedTerms))
dev.off()

# Treemap of reduced terms
pdf(paste0(outdir, kf.tissue,'_vs_',mus.tissue,'_GSEA_GOBP_overlappingGO-treemap_250719.pdf'),
    width = 8, height = 8)
print(treemapPlot(reducedTerms))
dev.off()

# Save reduced terms and basic frequency summaries
reducedTerms <- reducedTerms %>% arrange(cluster)
write.csv(reducedTerms, paste0(outdir, kf.tissue,'_vs_',mus.tissue,'_GSEA_GOBP_overlappingGO_reducedterms_250719.csv'))

freq.table <- as.data.frame(table(reducedTerms$cluster))
freq.table <- freq.table %>% filter(Freq > 2)
reducedTerms.thresh <- reducedTerms %>% filter(cluster %in% freq.table$Var1) 

length(unique(reducedTerms.thresh$parentTerm)) 
# Fat - 33 parent terms remaining
# Heart - 6 parent terms remaining
# Brain - 1 parent terms remaining -- dont need to filter
# Liver - 3 parent terms remaiing

# ------------------------------------------------------------------------------
# Choose GO terms of interest (by tissue) for bubble plotting
# ------------------------------------------------------------------------------

##### Brain GO terms
# go <- head(commonBP.simple$ID,10) 
go <- c('GO:0002449', 'GO:0002697', 'GO:0031349', 'GO:0002460', 'GO:0002253',
        'GO:0001819', 'GO:0002831', 'GO:0001906', 'GO:0019882', 'GO:0042110')
# GO:0002449 - lymphocyte mediated immunity
# GO:0002697 - regulation of immune effector process
# GO:0031349 - positive regulation of defense response
# GO:0002460 - adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
# GO:0002253 - activation of immune response
# GO:0001819 - positive regulation of cytokine production
# GO:0002831 - regulation of response to biotic stimulus
# GO:0001906 - cell killing
# GO:0019882 - antigen processing and presentation
# GO:0042110 - T cell activation

##### Heart vs Heart GO terms
# go <- c(head(commonBP.simple$ID,5), tail(commonBP.simple$ID,5))
go <- c('GO:0060326','GO:0032944','GO:0050670','GO:0045089','GO:0046631',
        'GO:1903555','GO:000611','GO:0006103','GO:0009060','GO:0032042')
# GO:0060326 - cell chemotaxis
# GO:0032944 - regulation of mononuclear cell proliferation
# GO:0050670 - regulation of lymphocyte proliferation
# GO:0045089 - positive regulation of innate immune response
# GO:0046631 - alpha-beta T cell activation
# GO:1903555 - regulation of TNF superfamily cytokine production
# GO:000611 - oxidative phosphorylation
# GO:0006103 - 2-oxoglutarate metabolic process
# GO:0009060 - aerobic respiration
# GO:0032042 - mitochondrial DNA metabolic process

##### Fat vs Mesenteric Fat GO terms
# go <- c(head(commonBP.simple$ID,5), tail(commonBP.simple$ID,5))
go <- c('GO:0002768','GO:0042113','GO:0002460','GO:0002292','GO:0035418',
        'GO:0030073','GO:0007269','GO:1905952','GO:0098754','GO:0042632')
# GO:0002768 - immune response-regulating cell surface receptor signaling pathway
# GO:0042113 - B cell activation
# GO:0002460 - adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
# GO:0002292 - T cell differentiation involved in immune response
# GO:0035418 - protein localization to synapse
# GO:0030073 - insulin secretion
# GO:0007269 - neurotransmitter secretion
# GO:1905952 - regulation of lipid localization
# GO:0098754 - detoxification
# GO:0042632 - cholesterol homeostasis

##### Liver GO terms
# go <- c(head(commonBP.simple$ID,5), tail(commonBP.simple$ID,5))
go <- c('GO:0030261','GO:1990089','GO:0042110','GO:0045088','GO:0002833',
        'GO:0002757','GO:0043410','GO:0001909','GO:0001819','GO:0043269')
# GO:0030261 - chromosome condensation
# GO:1990089 - response to nerve growth factor
# GO:0042110 - T cell activation
# GO:0045088 - regulation of innate immune response
# GO:0002833 - positive regulation of response to biotic stimulus
# GO:0002757 - immune response-activating signaling pathway
# GO:0043410 - positive regulation of MAPK cascade
# GO:0001909 - leukocyte mediated cytotoxicity
# GO:0001819 - positive regulation of cytokine production
# GO:0043269 - regulation of monoatomic ion transport

# ------------------------------------------------------------------------------
# Prep for bubble plot (common across tissues)
# ------------------------------------------------------------------------------

# Keep only the terms of interest
data_go <- filter(commonBP.simple, commonBP.simple$ID %in% go)

# Rename selected columns for species-specific labels
colnames(data_go)[4:7] = c('NES.mus', 'padj.mus', 'NES.kf', 'padj.kf')

# Long format with .value trick to split NES/padj by species
df_long <- data_go %>%
  pivot_longer(
    cols = c(NES.mus, padj.mus, NES.kf, padj.kf),
    names_to = c(".value", "species"),
    names_sep = "\\."
  )

# Order by summed score; set factor levels
df_long = df_long %>% arrange(desc(df_long$summed_score))
df_long$Description = factor(df_long$Description, levels = unique(df_long$Description))
df_long$species     = factor(df_long$species, levels = c('kf', 'mus'))

# ------------------------------------------------------------------------------
# Bubble plot (NES color, -log10 padj size) and save
# ------------------------------------------------------------------------------

# For brain and fat, width = 15 (per original comment)
pdf(paste0(outdir, kf.tissue,'_vs_',mus.tissue,'_GSEA_GOBP_selectoverlappingGOterms_250715.pdf'),
    width = 15, height = 10)
ggplot(data=df_long, aes(x=species, y=Description)) +
  geom_point(aes(color=NES, size = -log10(padj))) +
  scale_colour_gradient2(low = "blue", mid = "purple", high = "red",
                         midpoint = 0, space = "Lab",
                         na.value = "grey50", guide = "colourbar",
                         aesthetics = "colour")
dev.off()

# ------------------------------------------------------------------------------
# Summative analysis (counts of overlapping BP terms)
# ------------------------------------------------------------------------------

mouse.GO.res.dir = '/Output/'

mouse_files <- list.files(path = mouse.GO.res.dir, full.names = TRUE)
mouse_files <- mouse_files[grepl("250710_BPcommonGOIDs_", mouse_files, ignore.case = TRUE)]
# mouse_files <- mouse_files[-13]  # (optional exclusion per original)

mouse_data_list <- lapply(mouse_files, function(file) {
  df <- read.csv(file, row.names = 1)
  df$source_file <- basename(file)
  return(df)
})

mouse_df <- bind_rows(mouse_data_list)

# Parse tissue name from filename
mouse_df$tissue <- sub(
  pattern = "250710_BPcommonGOIDs_(.*)_atlas-vs-maca.*",
  replacement = "\\1",
  x = mouse_df$source_file
)

# Barplot: all overlapping BP terms per tissue
to_plot1 <- as.data.frame(table(mouse_df$tissue))
to_plot1$Var1 <- factor(to_plot1$Var1, levels = to_plot1$Var1[order(to_plot1$Freq)])

sum1 <- ggplot(to_plot1, aes(x=Var1, y=Freq)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_y_continuous(breaks = seq(0, 6500, by = 250))

pdf(paste0(outdir,'250715_atlas-vs-maca_barplot_allBPcommonGOIDs_bytissue.pdf'))
print(sum1)
dev.off()

# Barplot: padj < 0.1 in both species
mouse_df.filter <- mouse_df %>% filter(p.adjust.mus < 0.1 & p.adjust.kf < 0.1)
to_plot <- as.data.frame(table(mouse_df.filter$tissue))

# Ensure zero bars are shown for tissues with no entries under the filter
zeros  <- setdiff(to_plot1$Var1, to_plot$Var1)
to_add <- data.frame(Var1 = zeros, Freq = rep(0, length(zeros)))
to_plot <- rbind(to_plot, to_add)
to_plot$Var1 <- factor(to_plot$Var1, levels = to_plot$Var1[order(to_plot$Freq)])

sum2 <- ggplot(to_plot, aes(x=Var1, y=Freq)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_y_continuous(breaks = seq(0, 300, by = 50), limits = c(0,300))

pdf(paste0(outdir,'250715_atlas-vs-maca_barplot_pval-less-pt1_BPcommonGOIDs_bytissue.pdf'))
print(sum2)
dev.off()