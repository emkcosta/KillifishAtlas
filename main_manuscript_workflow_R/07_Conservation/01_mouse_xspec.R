# ------------------------------------------------------------------------------
# Title:          Cross-Species Comparison — Analysis of Schaum et al.
# Author:         Emma Costa
# Date:           Code compiled on 20250807
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description:    Reformats MACA mouse DE results, maps orthologs to killifish,
#                 computes combined ranks for concordant DE genes, exports tables,
#                 and generates summary boxplots plus mini barplots per tissue.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Set up
# ------------------------------------------------------------------------------

# Libraries
library('readxl')
library('dplyr')
library('biomaRt')
library('tidyr')
library('ggplot2')

# Working directory (use active RStudio document path)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# I/O directories
indir  = '/Input/'
outdir = '/Output/'

dir.create('/Output/plots/')
dir.create('/Output/tables/')

# ------------------------------------------------------------------------------
# Load ortholog mapping and relevant killifish atlas data
# ------------------------------------------------------------------------------

# Orthologs list
ortho       <- read.csv(paste0(indir,'nfur-ncbi_orthologs_20170302_pps.csv'))
ortho.mouse <- ortho[, c('N..furzeri..NCBI.', 'Mouse')]

# Killifish DESeq2 results and correlation results
cor.res <- read_excel(paste0(indir,'SupFile_4_CorrResults_SexCombined.xlsx'),
                      sheet = 'AllCorrRes_sexcombined')

load('../../common_robjects/dds_DESeq2_TPS_allsamples_agesbinned_240805.bin')

# ------------------------------------------------------------------------------
# Load mouse MACA dataset (DE analysis across ages) - already performed
# ------------------------------------------------------------------------------

load(paste0(indir, 'dds_MACA.bin')) # DE analysis all ages vs 3 months
# Note: MACA dataset used Black6 mice from NIA Aged Rodent Colonies

# ------------------------------------------------------------------------------
# Gene conversion (mouse Ensembl IDs → symbols via biomaRt)
# ------------------------------------------------------------------------------

# Retrieve mappings from Ensembl
mart    <- useDataset('mmusculus_gene_ensembl', useMart('ensembl'))
geneset <- rownames(dds)
G_list  <- getBM(filters= "ensembl_gene_id",
                 attributes= c("ensembl_gene_id","external_gene_name","description"),
                 values=geneset, mart= mart)

# Remove deprecated Ensembl IDs from dds
ids     = which(!(geneset %in% G_list$ensembl_gene_id))
leftout = geneset[ids]
dds     = dds[!rownames(dds)%in%leftout,]

# Filter to comparison of interest (27 vs 3), add gene annotations
dds.filt                   <- dds[, grepl("27_vs_3", colnames(dds))]
dds.filt$ensembl_gene_id   <- rownames(dds.filt)
dds.filt                   <- cbind(dds.filt,G_list)
dds.filt$description       = NULL
dds.filt$ensembl_gene_id   = NULL
dds.filt$ensembl_gene_id   <- rownames(dds.filt)

# ------------------------------------------------------------------------------
# Reformat MACA results to tidy/wide formats
# ------------------------------------------------------------------------------

# Pivot longer: keep Ensembl ID and gene symbol
df = dds.filt
df_long <- df %>%
  pivot_longer(
    cols = -c(ensembl_gene_id, external_gene_name),
    names_to = "condition",
    values_to = "value"
  )

# Parse tissue, comparison, metric from column names
df_long <- df_long %>%
  extract(condition, into = c("tissue", "comparison", "metric"),
          regex = "(.+?)_(\\d+_vs_\\d+)_(log2FoldChange|padj)") 

# Remove rows lacking gene symbols
df_long <- df_long[df_long$external_gene_name != "", ]

# Pivot wider to get log2FC and padj side-by-side
df_wide <- df_long %>%
  pivot_wider(
    names_from = metric,
    values_from = value
  )

# Persist and re-read (keeps a stable on-disk copy)
write.csv(df_wide, file = paste0(outdir,'tables/250708_reformat_maca_deseq2results.csv'), row.names = F)
df_wide <- read.csv(file = paste0(outdir,'tables/250708_reformat_maca_deseq2results.csv'))

# ------------------------------------------------------------------------------
# Prepare killifish→mouse mappings for correlation results
# ------------------------------------------------------------------------------

colnames(ortho.mouse)[1] = "killi_gene"
colnames(cor.res)[1]     = "killi_gene"
cor.res.convert = cor.res %>%
  left_join(ortho.mouse, by = "killi_gene")

# ------------------------------------------------------------------------------
# Harmonize column names and define tissue map
# ------------------------------------------------------------------------------

# Add suffix to mouse-derived columns and standardize gene column
colnames(df_wide) <- paste0(colnames(df_wide), ".mus")
colnames(df_wide)[2]      = 'Mouse.gene'
colnames(cor.res.convert)[10] = 'Mouse.gene'

# Manual tissue mapping: killifish ↔ mouse
tissue_map <- data.frame(
  kf_tissue = c("Bone", "Brain", "Fat", "Fat", "Heart",
                "Kidney", "Muscle", "Liver", "Kidney", "Fat",
                "Skin", "Gut", "Spleen", "Fat", "Kidney"),
  mus_tissue =  c("Bone","Brain","Brown_Fat","Gonadal_Fat","Heart",
                  "Kidney","Limb_Muscle","Liver","Marrow","Mesenteric_Fat",
                  "Skin","Small_Intestine","Spleen","Subcutaneous_Fat","White_Blood_Cells")
)

# Correlation color thresholds
cor_lowerbound = -0.5
cor_upperbound =  0.5

# ------------------------------------------------------------------------------
# Initialize collectors for per-tissue results
# ------------------------------------------------------------------------------

all_ups        <- list()
all_downs      <- list()
all_significant <- list()
all_unfiltered <- list()

# ------------------------------------------------------------------------------
# Main loop: per-tissue cross-species comparison and ranking
# ------------------------------------------------------------------------------

for (row in 1:nrow(tissue_map)) {
  # Tissues to compare
  kf.tissue  <- as.character(tissue_map[row,1])
  mus.tissue <- as.character(tissue_map[row,2])
  
  # Killifish correlation results (per tissue)
  # (Color points by correlation rather than filtering them out)
  cor.res.tissue <- cor.res %>% subset(tissue == kf.tissue)
  colnames(cor.res.tissue)[1] = 'gene_symbol'
  
  # Killifish DEGs (4_vs_1 bin; 50% surv too few samples → use next bin)
  res = as.data.frame(results_list_CA1[[kf.tissue]]$'4_vs_1'$resall)
  colnames(ortho.mouse)[1] = 'gene_symbol'
  res = res %>% left_join(ortho.mouse, by = 'gene_symbol')
  colnames(res)[c(1:6)] <- paste0(colnames(res)[c(1:6)], ".kf") # mark killifish-derived columns
  
  # Mouse DEGs (subset tissue; oldest vs 3mo)
  res.mouse = df_wide %>% filter(tissue.mus == mus.tissue)
  res.mouse = res.mouse %>% filter(Mouse.gene %in% res$Mouse)
  
  # Annotate killifish DEG trend and correlation color
  res <- res %>%
    mutate(trend.kf = case_when(
      log2FoldChange.kf >  0  ~ "up",
      log2FoldChange.kf < 0  ~ "down",
      TRUE              ~ NA_character_
    ))
  
  res$cor_spear.kf <- cor.res.tissue$cor_spear[match(res$gene_symbol, cor.res.tissue$gene_symbol)]
  res$cor_color.kf <- ifelse(
    res$cor_spear.kf < cor_lowerbound, "#85d0f9",
    ifelse(res$cor_spear.kf > cor_upperbound, "#f57a7f", NA)
  )
  res$cor_color.kf[is.na(res$cor_color.kf)] <- "gray50"
  
  res = as.data.frame(res)
  
  # Merge mouse and killifish results by ortholog
  colnames(res)[8] = 'Mouse.gene'
  toplot = res.mouse %>% left_join(res, by ='Mouse.gene')
  
  # Rank scores per species (effect size × significance)
  toplot <- toplot %>%
    mutate(
      rank_score_mouse = log2FoldChange.mus * -log10(padj.mus)
    )
  
  toplot <- toplot %>%
    mutate(
      rank_score_killifish = log2FoldChange.kf * -log10(padj.kf)
    )
  
  # Save unfiltered per-tissue table
  write.csv(toplot,file = paste0(outdir,'tables/250704_unfiltered-combinedrank_',kf.tissue,"_vs_", mus.tissue, ".csv"))
  
  # Prefilter to significant DEGs in both species
  toplot.filtered <- toplot %>% subset(padj.mus < 0.05 & padj.kf < 0.05)
  
  # Concordant up/down sets (+ tissue label)
  toplot.filtered.ups <- toplot.filtered %>% subset(log2FoldChange.mus > 0 & log2FoldChange.kf > 0 ) %>%
    mutate(tissue = mus.tissue)
  toplot.filtered.downs <- toplot.filtered %>% subset(log2FoldChange.mus < 0 & log2FoldChange.kf < 0 ) %>%
    mutate(tissue = mus.tissue)
  
  # Per-species and combined ranks (ups)
  toplot.filtered.ups <- toplot.filtered.ups %>%
    mutate(rank_killifish_top = dense_rank(dplyr::desc(rank_score_killifish)))
  toplot.filtered.ups <- toplot.filtered.ups %>%
    mutate(rank_mouse_top = dense_rank(dplyr::desc(rank_score_mouse)))
  toplot.filtered.ups <- toplot.filtered.ups %>%
    mutate(rank_combo_top = rank_killifish_top + rank_mouse_top)
  
  # Per-species and combined ranks (downs)
  toplot.filtered.downs <- toplot.filtered.downs %>%
    mutate(rank_killifish_bott = dense_rank(rank_score_killifish))
  toplot.filtered.downs <- toplot.filtered.downs %>%
    mutate(rank_mouse_bott = dense_rank(rank_score_mouse))
  toplot.filtered.downs <- toplot.filtered.downs %>%
    mutate(rank_combo_bott = rank_killifish_bott + rank_mouse_bott)
  
  # Save filtered per-tissue tables
  write.csv(toplot.filtered.downs,file = paste0(outdir,'tables/250704_filteredbypval-combinedrank-shareddown_',kf.tissue,"_vs_", mus.tissue, ".csv"))
  write.csv(toplot.filtered.ups,file = paste0(outdir,'tables/250704_filteredbypval-combinedrank-sharedup_',kf.tissue,"_vs_", mus.tissue, ".csv"))
  
  # Top-N (limit 200 for clarity)
  toplot_topn_combo  <- as.data.frame(toplot.filtered.ups)   %>% arrange(rank_combo_top)   %>% slice_head(n = 200)
  toplot_bottn_combo <- as.data.frame(toplot.filtered.downs) %>% arrange(rank_combo_bott)  %>% slice_head(n = 200)
  
  # If too few genes, skip plotting/exports
  if(length(rownames(toplot.filtered)) < 4){
    next
  } 
  
  # --------------------------------------------------------------------
  # Mini barplots (top 2 genes up/down) — saved per tissue to PDF
  # --------------------------------------------------------------------
  
  # UP mini barplot
  plot_data_up <- toplot_topn_combo[1:2,] %>%
    dplyr::select(Mouse.gene, log2FoldChange.mus, log2FoldChange.kf) %>%
    pivot_longer(
      cols = starts_with("log2FoldChange"),
      names_to = "species",
      values_to = "log2FC"
    )
  plot_data_up$species <- recode(plot_data_up$species,
                                 log2FoldChange.mus = "Mouse",
                                 log2FoldChange.kf  = "Killifish")
  p3 <- ggplot(plot_data_up, aes(x = Mouse.gene, y = log2FC, fill = species)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    theme_minimal() +
    labs(x = "Gene", y = "log2 Fold Change") +
    scale_fill_manual(values = c("Mouse" = "#f57a7f", "Killifish" = "#85d0f9")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')
  pdf(paste0(outdir,'plots/250704_barplot_top2-combinedrankup_',kf.tissue,"_vs_", mus.tissue, ".pdf"), width = 2, height = 2)
  print(p3)
  dev.off()
  
  # DOWN mini barplot
  plot_data_down <- toplot_bottn_combo[1:2,][1:2,] %>%
    dplyr::select(Mouse.gene, log2FoldChange.mus, log2FoldChange.kf) %>%
    pivot_longer(
      cols = starts_with("log2FoldChange"),
      names_to = "species",
      values_to = "log2FC"
    )
  plot_data_down$species <- recode(plot_data_down$species,
                                   log2FoldChange.mus = "Mouse",
                                   log2FoldChange.kf  = "Killifish")
  p4 <- ggplot(plot_data_down, aes(x = Mouse.gene, y = log2FC, fill = species)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    theme_minimal() +
    labs(x = "Gene", y = "log2 Fold Change") +
    scale_fill_manual(values = c("Human" = "#f57a7f", "Killifish" = "#85d0f9")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')
  pdf(paste0(outdir,'plots/250704_barplot_bott2-combinedrankdown_',kf.tissue,"_vs_", mus.tissue, ".pdf"), width = 2, height = 2)
  print(p4)
  dev.off()
  
  # Export top-200 tables
  write.csv(toplot_topn_combo, file = paste0(outdir,'tables/250704_top200-combinedrankup_',kf.tissue,"_vs_", mus.tissue, ".csv"))
  write.csv(toplot_bottn_combo,file = paste0(outdir,'tables/250704_top200-combinedrankdown_',kf.tissue,"_vs_", mus.tissue, ".csv"))
  
  # Save to collectors
  all_ups[[mus.tissue]]   <- toplot.filtered.ups
  all_downs[[mus.tissue]] <- toplot.filtered.downs
  
  toplot.filtered$tissue.mus <- mus.tissue
  toplot$tissue.mus          <- mus.tissue
  
  all_significant[[mus.tissue]] <- toplot.filtered
  all_unfiltered[[mus.tissue]]  <- toplot
}


# ------------------------------------------------------------------------------
# Combine, export, and order tissues for final summaries
# ------------------------------------------------------------------------------

# Combine all tissues (significant & unfiltered)
combined_significant <- bind_rows(all_significant)
combined_unfiltered  <- bind_rows(all_unfiltered)
write.csv(combined_significant, file = paste0(outdir,'tables/250712_maca_allsiggenes_alltissues.csv'))
write.csv(combined_unfiltered,  file = paste0(outdir,'tables/250712_maca_allgenes_alltissues.csv'))

# Tissue order for plotting
my_order <- c("Brain", # nervous system
              "Brown_Fat","Gonadal_Fat", "Subcutaneous_Fat" ,"Mesenteric_Fat", # endocrine 
              "Small_Intestine", "Liver", # digestive
              "Limb_Muscle", "Bone", "Skin", # musculoskeletal + integumentary
              "Heart", "White_Blood_Cells", "Kidney","Marrow", "Spleen" # cardiovascular + immune
)

# Combine up/down sets and export
combined_ups   <- bind_rows(all_ups)
combined_downs <- bind_rows(all_downs)
write.csv(combined_ups,   file = paste0(outdir,'tables/250712_maca_allUPsiggenes_alltissues.csv'))
write.csv(combined_downs, file = paste0(outdir,'tables/250712_maca_allDOWNsiggenes_alltissues.csv'))

# Factorize tissue with desired order
combined_ups$tissue   <- factor(combined_ups$tissue,   levels = my_order)
combined_downs$tissue <- factor(combined_downs$tissue, levels = my_order)

# Optional filtering for final plots
combined_ups   <- combined_ups   %>% subset(!(tissue %in% c('Bone')))
combined_downs <- combined_downs %>% subset(!(tissue %in% c('Bone','Marrow', 'Spleen')))

# ------------------------------------------------------------------------------
# Build plotting datasets (ups & downs)
# ------------------------------------------------------------------------------

plot_data_ups <- combined_ups %>%
  dplyr::select(Mouse.gene, tissue.mus, log2FoldChange.mus, log2FoldChange.kf, cor_color.kf) %>%
  pivot_longer(cols = starts_with("log2FoldChange"),
               names_to = "species", values_to = "log2FC") %>%
  mutate(species   = recode(species, log2FoldChange.mus = "Mouse", log2FoldChange.kf = "Killifish"))
plot_data_ups$dot_alpha = ifelse(plot_data_ups$cor_color.kf == '#f57a7f', 0.6, 0.2)

plot_data_downs <- combined_downs %>%
  dplyr::select(Mouse.gene, tissue.mus, log2FoldChange.mus, log2FoldChange.kf, cor_color.kf) %>%
  pivot_longer(cols = starts_with("log2FoldChange"),
               names_to = "species", values_to = "log2FC") %>%
  mutate(species   = recode(species, log2FoldChange.mus = "Mouse", log2FoldChange.kf = "Killifish"))
plot_data_downs$dot_alpha = ifelse(plot_data_downs$cor_color.kf == '#85d0f9', 0.6, 0.2)

# ------------------------------------------------------------------------------
# Subsystem-specific panels (nervous, endocrine, digestive, etc.)
# ------------------------------------------------------------------------------

nerv      <- c("Brain")                                           # nervous system
endocrine <- c("Brown_Fat","Gonadal_Fat", "Subcutaneous_Fat" ,"Mesenteric_Fat") # endocrine 
digestive <- c("Small_Intestine", "Liver")                        # digestive
musc.integ<- c("Limb_Muscle", "Bone", "Skin")                     # musculoskeletal + integumentary
card.imm  <- c("Heart", "White_Blood_Cells", "Kidney","Marrow", "Spleen") # cardiovascular + immune

# UPS panels by subsystem
plot_data_ups.filt = plot_data_ups %>% filter(tissue.mus %in% nerv)
p_1 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  theme_classic() +
  theme(legend.position="none")

plot_data_ups.filt = plot_data_ups %>% filter(tissue.mus %in% endocrine)
p_2 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  theme_classic() +
  theme(legend.position="none")

plot_data_ups.filt = plot_data_ups %>% filter(tissue.mus %in% digestive)
p_3 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  theme_classic() +
  theme(legend.position="none")

plot_data_ups.filt = plot_data_ups %>% filter(tissue.mus %in% musc.integ)
p_4 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  theme_classic() +
  theme(legend.position="none")

plot_data_ups.filt = plot_data_ups %>% filter(tissue.mus %in% card.imm)
p_5 <- ggplot(plot_data_ups.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  theme_classic() +
  theme(legend.position="none")

# DOWNS panels by subsystem
plot_data_downs.filt = plot_data_downs %>% filter(tissue.mus %in% nerv)
p_6 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  theme_classic() +
  theme(legend.position="none")

plot_data_downs.filt = plot_data_downs %>% filter(tissue.mus %in% endocrine)
p_7 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  theme_classic() +
  theme(legend.position="none")

plot_data_downs.filt = plot_data_downs %>% filter(tissue.mus %in% digestive)
p_8 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  theme_classic() +
  theme(legend.position="none")

plot_data_downs.filt = plot_data_downs %>% filter(tissue.mus %in% musc.integ)
p_9 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  theme_classic() +
  theme(legend.position="none")

plot_data_downs.filt = plot_data_downs %>% filter(tissue.mus %in% card.imm)
p_10 <- ggplot(plot_data_downs.filt, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs.filt$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs.filt$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  theme_classic() +
  theme(legend.position="none")

# ------------------------------------------------------------------------------
# Final combined boxplots across all tissues (ups & downs)
# ------------------------------------------------------------------------------

# UPS combined
p_up <- ggplot(plot_data_ups, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  theme_classic() +
  theme(legend.position="none")

pdf(paste0(outdir,'plots/250708_boxplot_top-combinedrankup_alltissue.pdf'), width = 10, height = 2)
print(p_up)
dev.off()

# DOWNS combined
p_down <- ggplot(plot_data_downs, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs$dot_alpha)


p_down <- ggplot(plot_data_downs, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue.mus, nrow = 1) +
  theme_classic() +
  theme(legend.position="none")

pdf(paste0(outdir,'plots/250708_boxplot_top-combinedrankdown_alltissue.pdf'), width = 10, height = 2)
print(p_down)
dev.off()