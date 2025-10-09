# ------------------------------------------------------------------------------
# Title:          GO Enrichment Overlap – Killifish vs. Human
# Author:         Emma Costa
# Date:           Code compiled on 20250812
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description:    Loads precomputed GO enrichment results for killifish and human,
#                 parses tissue identity from filenames, and computes overlapping
#                 GO terms by ontology (BP, CC, MF) across matched tissues.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------
# Load packages and set working directory
# ------------------------------------------------------------------

library('dplyr')
library('tidyr')

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
outdir = '/Output/tables/'

human.GO.res.dir = outdir
kf.GO.res.dir = outdir

# ------------------------------------------------------------------
# Load and parse killifish GO enrichment results
# ------------------------------------------------------------------

kf_files <- list.files(path = kf.GO.res.dir, full.names = TRUE)
kf_files_human <- kf_files[grepl("humanterms", kf_files, ignore.case = TRUE)]

kf_data_list <- lapply(kf_files_human, function(file) {
  df <- read.csv(file)
  df$source_file <- basename(file)
  return(df)
})

kf_df <- bind_rows(kf_data_list)

kf_df$tissue <- sub(
  pattern = ".*bin4-vs-bin1_(.*?)_humanterms.*",
  replacement = "\\1",
  x = kf_df$source_file
)

kf_df$source_file = NULL

# ------------------------------------------------------------------
# Load and parse human GO enrichment results
# ------------------------------------------------------------------
human_files <- list.files(path = human.GO.res.dir, full.names = TRUE)
human_files <- grep("250710_DESeq2_70svs20s_.*_humanterms.*", human_files, value = TRUE)


human_data_list <- lapply(human_files, function(file) {
  df <- read.csv(file)
  df$source_file <- basename(file)
  return(df)
})

human_df <- bind_rows(human_data_list)

human_df$tissue <- sub(
  pattern = "250710_DESeq2_70svs20s_(.*)_humanterms.*",
  replacement = "\\1",
  x = human_df$source_file
)

human_df$source_file = NULL

# ------------------------------------------------------------------
# Define tissue mapping between killifish and human
# ------------------------------------------------------------------

tissue_map <- data.frame(
  kf_tissue = c('Fat', 'Liver', 'Brain', 'Spleen', 'Skin', 'Heart', 
                'Muscle', 'SpinalCord', 
                'Kidney', 'Kidney',
                'Gut', 'Gut', 'Gut', 'Gut'),
  
  human_tissue = c("adipose_visceral_omentum", "liver", "brain_cortex", "spleen", "skin_not_sun_exposed_suprapubic", "heart_left_ventricle",
                   "muscle_skeletal", "brain_spinal_cord_cervical_c-1", 
                   "kidney_cortex", "whole_blood", 
                   "esophagus_muscularis", "small_intestine_terminal_ileum", "colon_transverse", "colon_sigmoid")
)

# ------------------------------------------------------------------
# Loop over tissue pairs and compute GO term overlaps
# ------------------------------------------------------------------

for (row in 1:nrow(tissue_map)) {
  kf.tissue <- as.character(tissue_map[row, 1])
  hum.tissue <- as.character(tissue_map[row, 2])
  
  kf.tissue.GO <- kf_df %>% filter(tissue == kf.tissue)
  hum.tissue.GO <- human_df %>% filter(tissue == hum.tissue)
  
  kf.tissue.GO.dir <- kf.tissue.GO
  hum.tissue.GO.dir <- hum.tissue.GO 
  
  colnames(kf.tissue.GO.dir)[c(4:13)] <- paste0(colnames(kf.tissue.GO.dir)[c(4:13)], ".kf")
  colnames(hum.tissue.GO.dir)[c(4:13)] <- paste0(colnames(hum.tissue.GO.dir)[c(4:13)], ".hum")
  
  # Join and output full overlap
  common_GO <- hum.tissue.GO.dir %>%
    left_join(kf.tissue.GO.dir, by = c('ONTOLOGY', 'ID', 'Description'))
  
  write.csv(common_GO, paste0(outdir, "250712_ALLcommonGOIDs_", kf.tissue, "-vs-", hum.tissue, "_atlas-vs-human.csv"))
  
  # Subset by ontology and export
  common_GO_BP <- common_GO %>% filter(ONTOLOGY == 'BP')
  write.csv(common_GO_BP, paste0(outdir, "250712_BPcommonGOIDs_", kf.tissue, "-vs-", hum.tissue, "_atlas-vs-human.csv"))
  
  common_GO_CC <- common_GO %>% filter(ONTOLOGY == 'CC')
  write.csv(common_GO_CC, paste0(outdir, "250712_CCcommonGOIDs_", kf.tissue, "-vs-", hum.tissue, "_atlas-vs-human.csv"))
  
  common_GO_MF <- common_GO %>% filter(ONTOLOGY == 'MF')
  write.csv(common_GO_MF, paste0(outdir, "250712_MFcommonGOIDs_", kf.tissue, "-vs-", hum.tissue, "_atlas-vs-human.csv"))
}
