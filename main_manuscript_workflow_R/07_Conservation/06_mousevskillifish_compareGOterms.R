# ------------------------------------------------------------------------------
# Title:          Overlap of GO enrichment between killifish and mouse
# Author:         
# Date:           Code compiled on 20250808
# Description:    Aggregates per-tissue GO GSEA results for killifish and mouse,
#                 parses tissue labels from filenames, and writes overlap tables
#                 (ALL, BP, CC, MF) for each mapped tissue pair.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

library('dplyr')
library('tidyr')

# Working directory (uses active RStudio document path)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# I/O directories
indir  = '/Input/'
outdir = '/Output/tables/'

# Source result directories
maca.GO.res.dir = outdir  # mouse GO GSEA results
kf.GO.res.dir   = outdir  # killifish GO GSEA results

# ------------------------------------------------------------------------------
# Build mega data frame: killifish GO results
# ------------------------------------------------------------------------------

# List all files; keep those containing "mouseterms"
kf_files       <- list.files(path = kf.GO.res.dir, full.names = TRUE)
kf_files_mouse <- kf_files[grepl("mouseterms", kf_files, ignore.case = TRUE)]

# Read and tag with source filename
kf_data_list <- lapply(kf_files_mouse, function(file) {
  df <- read.csv(file)
  df$source_file <- basename(file)  # keep just filename
  return(df)
})

# Combine all killifish GO results
kf_df <- bind_rows(kf_data_list)

# Extract tissue name from filename pattern
kf_df$tissue <- sub(
  pattern = ".*bin4-vs-bin1_(.*?)_mouseterms.*",
  replacement = "\\1",
  x = kf_df$source_file
)

# Drop helper column
kf_df$source_file = NULL

# ------------------------------------------------------------------------------
# Build mega data frame: mouse GO results
# ------------------------------------------------------------------------------

maca_files <- list.files(path = maca.GO.res.dir, full.names = TRUE)
maca_files_mouse <- maca_files[grepl("mouse_GOGSEA", kf_files, ignore.case = TRUE)]

# Read and tag with source filename
maca_data_list <- lapply(maca_files_mouse, function(file) {
  df <- read.csv(file)
  df$source_file <- basename(file)
  return(df)
})

# Combine all mouse GO results
maca_df <- bind_rows(maca_data_list)

# Extract tissue name from filename pattern
maca_df$tissue <- sub(
  pattern = ".*GOGSEA_(.*)\\.csv",
  replacement = "\\1",
  x = maca_df$source_file
)

# Drop helper column
maca_df$source_file = NULL

# ------------------------------------------------------------------------------
# Overlap: map tissues and join by GO term
# ------------------------------------------------------------------------------

tissue_map <- data.frame(
  kf_tissue  = c("Bone", "Brain", "Fat", "Fat", "Heart",
                 "Kidney", "Muscle", "Liver", "Kidney", "Fat",
                 "Skin", "Gut", "Spleen", "Fat", "Kidney"),
  mus_tissue = c("Bone","Brain","Brown_Fat","Gonadal_Fat","Heart",
                 "Kidney","Limb_Muscle","Liver","Marrow","Mesenteric_Fat",
                 "Skin","Small_Intestine","Spleen","Subcutaneous_Fat","White_Blood_Cells")
)

# directions = c('upregulated', 'downregulated')  # (kept as comment per original)

for (row in 1:nrow(tissue_map)) {
  # Select paired tissues
  kf.tissue  <- as.character(tissue_map[row,1])
  mus.tissue <- as.character(tissue_map[row,2])
  
  # Subset GO results for each species/tissue
  kf.tissue.GO  <- kf_df   %>% filter(tissue == kf.tissue)
  mus.tissue.GO <- maca_df %>% filter(tissue == mus.tissue)
  
  # Duplicate as ".dir" objects (per original)
  kf.tissue.GO.dir  <- kf.tissue.GO
  mus.tissue.GO.dir <- mus.tissue.GO 
  
  # Add species suffixes to shared columns (positions 4:13 as in original)
  colnames(kf.tissue.GO.dir)[c(4:13)]  <- paste0(colnames(kf.tissue.GO.dir)[c(4:13)],  ".kf")
  colnames(mus.tissue.GO.dir)[c(4:13)] <- paste0(colnames(mus.tissue.GO.dir)[c(4:13)], ".mus")
  
  # Left-join on GO identifiers
  common_GO <- mus.tissue.GO.dir %>%
    left_join(kf.tissue.GO.dir, by = c('ONTOLOGY','ID','Description'))
  
  # Write overlap tables (ALL, and by ontology)
  write.csv(common_GO,    paste0(outdir, "250710_ALLcommonGOIDs_", kf.tissue,"-vs-",mus.tissue,"_atlas-vs-maca.csv"))
  
  common_GO_BP <- common_GO %>% filter(ONTOLOGY == 'BP')
  write.csv(common_GO_BP, paste0(outdir, "250710_BPcommonGOIDs_", kf.tissue,"-vs-",mus.tissue,"_atlas-vs-maca.csv"))
  
  common_GO_CC <- common_GO %>% filter(ONTOLOGY == 'CC')
  write.csv(common_GO_CC, paste0(outdir, "250710_CCcommonGOIDs_", kf.tissue,"-vs-",mus.tissue,"_atlas-vs-maca.csv"))
  
  common_GO_MF <- common_GO %>% filter(ONTOLOGY == 'MF')
  write.csv(common_GO_MF, paste0(outdir, "250710_MFcommonGOIDs_", kf.tissue,"-vs-",mus.tissue,"_atlas-vs-maca.csv"))
}
