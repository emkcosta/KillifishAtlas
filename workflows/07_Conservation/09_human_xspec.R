# ------------------------------------------------------------------------------
# Title:          Cross-Species Comparison — Analysis of GTEx data v10 release
# Author:         Emma Costa
# Date:           Code compiled on 20250807
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description:    Reformats GTEx human DE results, maps orthologs to killifish,
#                 computes combined ranks for concordant DE genes, exports tables,
#                 and generates summary boxplots plus mini barplots per tissue.
# ------------------------------------------------------------------------------


#------------------------------------------------
# Set up
#------------------------------------------------
library('dplyr')
library('readxl')
library('ggplot2')
library('tidyr')
library('biomaRt')

# Working directory (active RStudio document path)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# I/O directories
indir  = '/Input/'
outdir = '/Output/'

# ------------------------------------------------------------------------------
# Load GTEx data (DE analysis across ages)
# ------------------------------------------------------------------------------
# load GTEx metadata and DESeq2 results
gtex.meta = read.csv(paste0(indir,'GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.csv'))
load(file='../common_robjects/dds_DESeq2_gtex_results_250709.bin') #deseq results # 50% survival for humans 75-85 - https://www.cdc.gov/nchs/data/nvsr/nvsr69/NVSR69-08-508.pdf

# ------------------------------------------------------------------------------
# Load ortholog mapping and relevant killifish atlas data
# ------------------------------------------------------------------------------
# load orthologs list
ortho = read.csv(paste0(indir,'nfur-ncbi_orthologs_20170302_pps.csv'))
ortho.human = ortho[,c('N..furzeri..NCBI.','Human')] #using "Human" column
rm(ortho)

# load data from biomaRt for ENS ID conversion
mart <- useDataset('hsapiens_gene_ensembl', useMart('ensembl'))

# load atlas deseq results
load('../common_robjects/dds_DESeq2_TPS_allsamples_agesbinned_240805.bin')
load('../common_robjects/dds_TPS_alltissue_femalesOnly_DESeq2results_241106.bin') #this is to get the gonad separated
load('../common_robjects/dds_TPS_alltissue_malesOnly_DESeq2results_241106.bin') #this is to get the gonad separated

# load atlas correlation results
cor.res <- read_excel(paste0(indir,'SupFile_4_CorrResults_SexCombined.xlsx',sheet = 'AllCorrRes_sexcombined'))
cor.res.ovary <- read_excel(paste0(indir,'SupFile_3_CorrResults_SexSplit.xlsx', sheet = "Gonad_Ovary_F"))
cor.res.testis <- read_excel(paste0(indir,'SupFile_3_CorrResults_SexSplit.xlsx', sheet = "Gonad_Testis_M"))
cor.res.ovary$...1 = NULL
rownames(cor.res.ovary) = cor.res.ovary$gene
cor.res.testis$...1 = NULL
rownames(cor.res.testis) = cor.res.testis$gene

# ------------------------------------------------------------------------------
# Harmonize column names and define tissue map
# ------------------------------------------------------------------------------
#brain - neuroanatomical similarities: https://pmc.ncbi.nlm.nih.gov/articles/PMC10464506/ - telencephalon
#gut - lack a stomach https://link.springer.com/article/10.1007/s00441-023-03845-8, esophagus directly to intestine
#kidney medulla doesn't have the 70-79 age group - so I will remove

tissue_map <- data.frame(
  kf_tissue = c('Fat', 'Liver', 'Brain', 'Spleen', 'Skin', 'Heart', 
                'Muscle', 'SpinalCord', 
                'Kidney', 'Kidney',
                'Gut', 'Gut', 'Gut', 'Gut'), #exclude eye and bone bc no equivalent & gonad bc sex combo
  
  human_tissue = c("adipose_visceral_omentum", "liver","brain_cortex", "spleen", "skin_not_sun_exposed_suprapubic","heart_left_ventricle",
                   "muscle_skeletal", "brain_spinal_cord_cervical_c-1", 
                   "kidney_cortex", "whole_blood", 
                   "esophagus_muscularis", "small_intestine_terminal_ileum", "colon_transverse", "colon_sigmoid")  # human equivalent
)

# Correlation color thresholds
cor_lowerbound = -0.5
cor_upperbound = 0.5

# ------------------------------------------------------------------------------
# Initialize collectors for per-tissue results
# ------------------------------------------------------------------------------
# collect all tissues
all_ups <- list()
all_downs <- list()
all_significant <- list()
all_unfiltered <- list()

# ------------------------------------------------------------------------------
# Main loop: per-tissue cross-species comparison and ranking
# ------------------------------------------------------------------------------

for (row in 1:nrow(tissue_map)) {
  #get the tissues you're going to compare
  kf.tissue <- as.character(tissue_map[row,1])
  hum.tissue <- as.character(tissue_map[row,2])
  
  #get the correlation results for that tissue from the killifish dataset
  #instead of filtering out the correlated genes, I can just color them
  cor.res.tissue <- cor.res %>% subset(tissue == kf.tissue)
  colnames(cor.res.tissue)[1] = 'gene_symbol'
  
  # fish DEG res
  res = as.data.frame(results_list_CA1[[kf.tissue]]$'4_vs_1'$resall) #~50% surv for fish has too few samples so went 1 beyond
  colnames(ortho.human)[1] = 'gene_symbol'
  res = res %>% left_join(ortho.human, by = 'gene_symbol')
  colnames(res)[c(1:6)] <- paste0(colnames(res)[c(1:6)], ".kf") #add a suffix to specify killifish derived
  
  # human DEG res
  res.human = as.data.frame(results_list_gtex[[hum.tissue]]$'70-79_vs_20-29'$resall) #~50% survival & oldest group available
  
  #convert ENS to gene IDs
  rownames(res.human) = sub("\\..*", "", rownames(res.human))
  res.human$gene_symbol = rownames(res.human)
  colnames(res.human)[7] = 'ensembl_gene_id'
  geneset <- rownames(res.human)
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                            "external_gene_name", "description"),values=geneset,mart= mart)

  ids = which(!(geneset %in% G_list$ensembl_gene_id)) # these are deprecated
  leftout = geneset[ids] # these are deprecated
  res.human = res.human[!rownames(res.human)%in%leftout,] # remove the deprecated gene names for ease

  #Replace Ensembl IDs with gene names where available
  id_to_name <- setNames(G_list$external_gene_name, G_list$ensembl_gene_id)
  new_rownames <- id_to_name[rownames(res.human)]
  new_rownames[is.na(new_rownames)] <- rownames(res.human)[is.na(new_rownames)]
  res.human$gene_symbol <- new_rownames

  #filter raw counts to only include genes that have KF orthologs
  common_genes = intersect(res.human$gene_symbol, ortho.human$Human)
  res.human <- res.human %>% filter(gene_symbol %in% common_genes)
  
  #ENS IDs without a gene name, remove them
  res.human = res.human %>% filter(!(gene_symbol == ''))
  res.human = res.human %>% filter(gene_symbol %in% res$Human)
  colnames(res.human)[c(1:7)] <- paste0(colnames(res.human)[c(1:7)], ".hum") #add a suffix to specify human derived
  
  ##### plotting results for killifish
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
  
  ##### add results for human 
  colnames(res.human)[8] = 'Human'
  toplot = res.human %>% left_join(res, by ='Human')
  
  ##### calculate scores
  toplot <- toplot %>%
    mutate(
      rank_score_human = log2FoldChange.hum * -log10(padj.hum)
    )
  
  toplot <- toplot %>%
    mutate(
      rank_score_killifish = log2FoldChange.kf * -log10(padj.kf)
    )
  
  write.csv(toplot,file = paste0(outdir,'tables/250710_unfiltered-combinedrank_',kf.tissue,"_vs_", hum.tissue, ".csv"))
  
  toplot.filtered <- toplot %>% subset(padj.hum < 0.05 & padj.kf < 0.05) #prefilter out genes not padj < 0.05
  
  toplot.filtered.ups <- toplot.filtered %>% subset(log2FoldChange.hum > 0 & log2FoldChange.kf > 0 ) %>%
    mutate(tissue = hum.tissue)  #both going up
  toplot.filtered.downs <- toplot.filtered %>% subset(log2FoldChange.hum < 0 & log2FoldChange.kf < 0 ) %>%
    mutate(tissue = hum.tissue) #both going down
  
  
  toplot.filtered.ups <- toplot.filtered.ups %>%
    mutate(rank_killifish_top = dense_rank(dplyr::desc(rank_score_killifish)))
  
  toplot.filtered.downs <- toplot.filtered.downs %>%
    mutate(rank_killifish_bott = dense_rank(rank_score_killifish))
  
  toplot.filtered.ups <- toplot.filtered.ups %>%
    mutate(rank_human_top = dense_rank(dplyr::desc(rank_score_human)))
  
  toplot.filtered.downs <- toplot.filtered.downs %>%
    mutate(rank_human_bott = dense_rank(rank_score_human))
  
  toplot.filtered.ups <- toplot.filtered.ups %>%
    mutate(rank_combo_top = rank_killifish_top + rank_human_top)
  
  toplot.filtered.downs <- toplot.filtered.downs %>%
    mutate(rank_combo_bott = rank_killifish_bott + rank_human_bott)
   
  write.csv(toplot.filtered.downs,file = paste0(outdir,'tables/250710_filteredbypval-combinedrank-shareddown_',kf.tissue,"_vs_", hum.tissue, ".csv"))
  write.csv(toplot.filtered.ups,file = paste0(outdir,'tables/250710_filteredbypval-combinedrank-sharedup_',kf.tissue,"_vs_", hum.tissue, ".csv"))
   
  
  ############
  toplot_topn_combo <- toplot.filtered.ups %>%
    arrange(rank_combo_top) %>%   # sort ascending
    slice_head(n = 200)    # not gonna go higher than 200 bc then the signs get confusing, note some tissues dont even have that many DEGs
  
  toplot_bottn_combo <- toplot.filtered.downs %>%
    arrange(rank_combo_bott) %>%   # sort ascending
    slice_head(n = 200)    # not gonna go higher than 200 bc then the signs get confusing
  
  
  ############
  if(length(rownames(toplot.filtered)) < 4){
    next  # skip to next iteration of for loop
  } 
  
  #plot the mini barplots
  plot_data_up <- toplot_topn_combo[1:2,] %>%
    dplyr::select(Human, log2FoldChange.hum, log2FoldChange.kf) %>%
    pivot_longer(
      cols = starts_with("log2FoldChange"),
      names_to = "species",
      values_to = "log2FC"
    )
  
  plot_data_up$species <- recode(plot_data_up$species,
                                 log2FoldChange.hum = "Human",
                                 log2FoldChange.kf  = "Killifish")
  
  p3 <- ggplot(plot_data_up, aes(x = Human, y = log2FC, fill = species)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    theme_minimal() +
    labs(
      x = "Gene",
      y = "log2 Fold Change"
    ) +
    scale_fill_manual(values = c("Human" = "#f57a7f", "Killifish" = "#85d0f9")) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = 'none'
    )
  
  pdf(paste0(outdir,'plots/250710_barplot_top2-combinedrankup_',kf.tissue,"_vs_", hum.tissue, ".pdf"), width = 2, height = 2)
  print(p3)
  dev.off()
  
  plot_data_down <- toplot_bottn_combo[1:2,][1:2,] %>%
    dplyr::select(Human, log2FoldChange.hum, log2FoldChange.kf) %>%
    pivot_longer(
      cols = starts_with("log2FoldChange"),
      names_to = "species",
      values_to = "log2FC"
    )
  
  plot_data_down$species <- recode(plot_data_down$species,
                                   log2FoldChange.hum = "Human",
                                   log2FoldChange.kf  = "Killifish")
  
  p4 <- ggplot(plot_data_down, aes(x = Human, y = log2FC, fill = species)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    theme_minimal() +
    labs(
      x = "Gene",
      y = "log2 Fold Change"
    ) +
    scale_fill_manual(values = c("Human" = "#f57a7f", "Killifish" = "#85d0f9")) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = 'none'
    )
  
  pdf(paste0(outdir,'plots/250710_barplot_bott2-combinedrankdown_',kf.tissue,"_vs_", hum.tissue, ".pdf"), width = 2, height = 2)
  print(p4)
  dev.off()
  
  
  if(length(rownames(toplot.filtered)) < 400){
     next  # skip to next iteration of for loop
  }

  #export these data
  write.csv(toplot_topn_combo,file = paste0(outdir,'tables/250710_top200-combinedrankup_',kf.tissue,"_vs_", hum.tissue, ".csv"))
  write.csv(toplot_bottn_combo,file = paste0(outdir,'tables/250710_top200-combinedrankdown_',kf.tissue,"_vs_", hum.tissue, ".csv"))

  
  # Save to list
  all_ups[[hum.tissue]] <- toplot.filtered.ups
  all_downs[[hum.tissue]] <- toplot.filtered.downs
  
  toplot.filtered$tissue.hum <- hum.tissue
  toplot$tissue.hum <- hum.tissue
  
  all_significant[[hum.tissue]] <- toplot.filtered
  all_unfiltered[[hum.tissue]] <- toplot

}


# ------------------------------------------------------------------------------
# Main loop: per-tissue cross-species comparison and ranking, Gonads
# ------------------------------------------------------------------------------
tissue_map <- data.frame(
  kf_tissue = c('Ovary', 'Testis'), 
  human_tissue = c("ovary","testis")  # human equivalent
)


for (row in 1:nrow(tissue_map)) {
  #get the tissues you're going to compare
  kf.tissue <- as.character(tissue_map[row,1])
  hum.tissue <- as.character(tissue_map[row,2])
  
  if(kf.tissue == 'Ovary'){
    #get the correlation results for that tissue from the killifish dataset
    #instead of filtering out the correlated genes, I can just color them
    cor.res.tissue <- cor.res.ovary
    colnames(cor.res.tissue)[1] = 'gene_symbol'
    
    # fish DEG res
    res = as.data.frame(results_list_female$Gonad$'4_vs_1'$resall) #~50% surv for fish has too few samples so went 1 beyond
    colnames(ortho.human)[1] = 'gene_symbol'
    res = res %>% left_join(ortho.human, by = 'gene_symbol')
    colnames(res)[c(1:6)] <- paste0(colnames(res)[c(1:6)], ".kf") #add a suffix to specify killifish derived
  } else {
    #get the correlation results for that tissue from the killifish dataset
    #instead of filtering out the correlated genes, I can just color them
    cor.res.tissue <- cor.res.testis
    colnames(cor.res.tissue)[1] = 'gene_symbol'
    
    # fish DEG res
    res = as.data.frame(results_list_male$Gonad$'4_vs_1'$resall) #~50% surv for fish has too few samples so went 1 beyond
    colnames(ortho.human)[1] = 'gene_symbol'
    res = res %>% left_join(ortho.human, by = 'gene_symbol')
    colnames(res)[c(1:6)] <- paste0(colnames(res)[c(1:6)], ".kf") #add a suffix to specify killifish derived
    
  }
  
  # human DEG res
  res.human = as.data.frame(results_list_gtex[[hum.tissue]]$'70-79_vs_20-29'$resall) #~50% survival & oldest group available
  
  #convert ENS to gene IDs
  rownames(res.human) = sub("\\..*", "", rownames(res.human))
  res.human$gene_symbol = rownames(res.human)
  colnames(res.human)[7] = 'ensembl_gene_id'
  geneset <- rownames(res.human)
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                            "external_gene_name", "description"),values=geneset,mart= mart)
  
  ids = which(!(geneset %in% G_list$ensembl_gene_id)) # these are deprecated
  leftout = geneset[ids] # these are deprecated
  res.human = res.human[!rownames(res.human)%in%leftout,] # remove the deprecated gene names for ease
  
  #Replace Ensembl IDs with gene names where available
  id_to_name <- setNames(G_list$external_gene_name, G_list$ensembl_gene_id)
  new_rownames <- id_to_name[rownames(res.human)]
  new_rownames[is.na(new_rownames)] <- rownames(res.human)[is.na(new_rownames)]
  res.human$gene_symbol <- new_rownames
  
  #filter raw counts to only include genes that have KF orthologs
  common_genes = intersect(res.human$gene_symbol, ortho.human$Human)
  res.human <- res.human %>% filter(gene_symbol %in% common_genes)
  
  #ENS IDs without a gene name, remove them
  res.human = res.human %>% filter(!(gene_symbol == ''))
  res.human = res.human %>% filter(gene_symbol %in% res$Human)
  colnames(res.human)[c(1:7)] <- paste0(colnames(res.human)[c(1:7)], ".hum") #add a suffix to specify human derived
  
  
  ##### plotting results for killifish
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
  
  
  ##### add results for human 
  colnames(res.human)[8] = 'Human'
  toplot = res.human %>% left_join(res, by ='Human')
  
  ##### calculate scores
  toplot <- toplot %>%
    mutate(
      rank_score_human = log2FoldChange.hum * -log10(padj.hum)
    )
  
  toplot <- toplot %>%
    mutate(
      rank_score_killifish = log2FoldChange.kf * -log10(padj.kf)
    )
  
  write.csv(toplot,file = paste0(outdir,'tables/250710_unfiltered-combinedrank_',kf.tissue,"_vs_", hum.tissue, ".csv"))
  
  toplot.filtered <- toplot %>% subset(padj.hum < 0.05 & padj.kf < 0.05) #prefilter out genes not padj < 0.05
  
  toplot.filtered.ups <- toplot.filtered %>% subset(log2FoldChange.hum > 0 & log2FoldChange.kf > 0 ) %>%
    mutate(tissue = hum.tissue)  #both going up
  toplot.filtered.downs <- toplot.filtered %>% subset(log2FoldChange.hum < 0 & log2FoldChange.kf < 0 ) %>%
    mutate(tissue = hum.tissue) #both going down
  
  
  toplot.filtered.ups <- toplot.filtered.ups %>%
    mutate(rank_killifish_top = dense_rank(dplyr::desc(rank_score_killifish)))
  
  toplot.filtered.downs <- toplot.filtered.downs %>%
    mutate(rank_killifish_bott = dense_rank(rank_score_killifish))
  
  toplot.filtered.ups <- toplot.filtered.ups %>%
    mutate(rank_human_top = dense_rank(dplyr::desc(rank_score_human)))
  
  toplot.filtered.downs <- toplot.filtered.downs %>%
    mutate(rank_human_bott = dense_rank(rank_score_human))
  
  toplot.filtered.ups <- toplot.filtered.ups %>%
    mutate(rank_combo_top = rank_killifish_top + rank_human_top)
  
  toplot.filtered.downs <- toplot.filtered.downs %>%
    mutate(rank_combo_bott = rank_killifish_bott + rank_human_bott)
  
  write.csv(toplot.filtered.downs,file = paste0(outdir,'tables/250710_filteredbypval-combinedrank-shareddown_',kf.tissue,"_vs_", hum.tissue, ".csv"))
  write.csv(toplot.filtered.ups,file = paste0(outdir,'tables/250710_filteredbypval-combinedrank-sharedup_',kf.tissue,"_vs_", hum.tissue, ".csv"))
  
  
  ############
  toplot_topn_combo <- toplot.filtered.ups %>%
    arrange(rank_combo_top) %>%   # sort ascending
    slice_head(n = 200)    # not gonna go higher than 200 bc then the signs get confusing, note some tissues dont even have that many DEGs
  
  toplot_bottn_combo <- toplot.filtered.downs %>%
    arrange(rank_combo_bott) %>%   # sort ascending
    slice_head(n = 200)    # not gonna go higher than 200 bc then the signs get confusing
  
  
  ############

  #plot the mini barplots
  plot_data_up <- toplot_topn_combo[1:2,] %>%
    dplyr::select(Human, log2FoldChange.hum, log2FoldChange.kf) %>%
    pivot_longer(
      cols = starts_with("log2FoldChange"),
      names_to = "species",
      values_to = "log2FC"
    )
  
  plot_data_up$species <- recode(plot_data_up$species,
                                 log2FoldChange.hum = "Human",
                                 log2FoldChange.kf  = "Killifish")
  
  p3 <- ggplot(plot_data_up, aes(x = Human, y = log2FC, fill = species)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    theme_minimal() +
    labs(
      x = "Gene",
      y = "log2 Fold Change"
    ) +
    scale_fill_manual(values = c("Human" = "#f57a7f", "Killifish" = "#85d0f9")) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = 'none'
    )
  
  pdf(paste0(outdir,'plots/250710_barplot_top2-combinedrankup_',kf.tissue,"_vs_", hum.tissue, ".pdf"), width = 2, height = 2)
  print(p3)
  dev.off()
  
  plot_data_down <- toplot_bottn_combo[1:2,][1:2,] %>%
    dplyr::select(Human, log2FoldChange.hum, log2FoldChange.kf) %>%
    pivot_longer(
      cols = starts_with("log2FoldChange"),
      names_to = "species",
      values_to = "log2FC"
    )
  
  plot_data_down$species <- recode(plot_data_down$species,
                                   log2FoldChange.hum = "Human",
                                   log2FoldChange.kf  = "Killifish")
  
  p4 <- ggplot(plot_data_down, aes(x = Human, y = log2FC, fill = species)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    theme_minimal() +
    labs(
      x = "Gene",
      y = "log2 Fold Change"
    ) +
    scale_fill_manual(values = c("Human" = "#f57a7f", "Killifish" = "#85d0f9")) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = 'none'
    )
  
  pdf(paste0(outdir,'plots/250710_barplot_bott2-combinedrankdown_',kf.tissue,"_vs_", hum.tissue, ".pdf"), width = 2, height = 2)
  print(p4)
  dev.off()

  #export these data
  write.csv(toplot_topn_combo,file = paste0(outdir,'tables/250710_top200-combinedrankup_',kf.tissue,"_vs_", hum.tissue, ".csv"))
  write.csv(toplot_bottn_combo,file = paste0(outdir,'tables/250710_top200-combinedrankdown_',kf.tissue,"_vs_", hum.tissue, ".csv"))
   
  
  # Save to list
  all_ups[[hum.tissue]] <- toplot.filtered.ups
  all_downs[[hum.tissue]] <- toplot.filtered.downs
  
  toplot.filtered$tissue.hum <- hum.tissue
  toplot$tissue.hum <- hum.tissue
  
  all_significant[[hum.tissue]] <- toplot.filtered
  all_unfiltered[[hum.tissue]] <- toplot
  
}


# ------------------------------------------------------------------------------
# Combine, export, and order tissues for final summaries
# ------------------------------------------------------------------------------
# Combine all tissues
combined_significant <- bind_rows(all_significant)
combined_unfiltered <- bind_rows(all_unfiltered)
write.csv(combined_significant, file = paste0(outdir,'tables/250712_gtex_allsiggenes_alltissues.csv'))
write.csv(combined_unfiltered, file = paste0(outdir,'tables/250712_gtex_allgenes_alltissues.csv'))

# Tissue order for plotting
my_order <- c("brain_cortex", #nervous system
              "adipose_visceral_omentum",  "ovary", "testis", #endocrine + reproductive
              "esophagus_muscularis", "colon_sigmoid", "colon_transverse", "small_intestine_terminal_ileum", "liver", #digestive
              "muscle_skeletal", "skin_not_sun_exposed_suprapubic", #musculoskeletal + integumentary
              "heart_left_ventricle", "whole_blood", "spleen" #cardiovascular + immune
              )

# Combine up/down sets and export
combined_ups <- bind_rows(all_ups)
combined_downs <- bind_rows(all_downs)
write.csv(combined_ups, file = paste0(outdir,'tables/250712_gtex_allUPsiggenes_alltissues.csv'))
write.csv(combined_downs, file = paste0(outdir,'tables/250712_gtex_allDOWNgenes_alltissues.csv'))

# Factorize tissue with desired order
combined_ups$tissue <- factor(combined_ups$tissue, levels = my_order)
combined_downs$tissue <- factor(combined_downs$tissue, levels = my_order)

# ------------------------------------------------------------------------------
# Build plotting datasets (ups & downs)
# ------------------------------------------------------------------------------
plot_data_ups <- combined_ups %>%
  dplyr::select(Human, tissue, log2FoldChange.hum, log2FoldChange.kf, cor_color.kf) %>%
  pivot_longer(cols = starts_with("log2FoldChange"),
               names_to = "species",
               values_to = "log2FC") %>%
  mutate(species = recode(species,
                          log2FoldChange.hum = "Human",
                          log2FoldChange.kf = "Killifish"))

plot_data_ups$dot_alpha = ifelse(plot_data_ups$cor_color.kf == '#f57a7f', 0.6, 0.2)


plot_data_downs <- combined_downs %>%
  dplyr::select(Human, tissue, log2FoldChange.hum, log2FoldChange.kf, cor_color.kf) %>%
  pivot_longer(cols = starts_with("log2FoldChange"),
               names_to = "species",
               values_to = "log2FC") %>%
  mutate(species = recode(species,
                          log2FoldChange.hum = "Human",
                          log2FoldChange.kf = "Killifish"))


plot_data_downs$dot_alpha = ifelse(plot_data_downs$cor_color.kf == '#85d0f9', 0.6, 0.2)


# Final wide plot for UPS
p_up <- ggplot(plot_data_ups, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_ups$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_ups$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA, size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue, nrow = 1) +
  theme_classic() +
  theme(legend.position="none")


pdf(paste0(outdir,'plots/boxplots_combined/250710_boxplot_top-combinedrankup_alltissue.pdf'), width = 10, height = 2)
print(p_up)
dev.off()

# Final wide plot for DOWNS
p_down <- ggplot(plot_data_downs, aes(x = species, y = log2FC, color = species)) +
  geom_jitter(color = plot_data_downs$cor_color.kf, width = 0.2, size = 0.5, alpha = plot_data_downs$dot_alpha) +
  geom_boxplot(alpha = 0, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tissue, nrow = 1) +
  theme_classic() +
  theme(legend.position="none")

pdf(paste0(outdir,'plots/boxplots_combined/250708_boxplot_top-combinedrankdown_alltissue.pdf'), width = 10, height = 2)
print(p_down)
dev.off()




