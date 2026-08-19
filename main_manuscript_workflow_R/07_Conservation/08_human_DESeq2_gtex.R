# ------------------------------------------------------------------------------
# Title:          GTEx processing and DESeq2 analysis (AGE ± SEX) by tissue
# Author:         Emma Costa
# Date:           Code compiled on 20250807
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description:    Reads GTEx raw counts per tissue, builds DESeq2 objects,
#                 saves normalized counts and metadata, then runs DE between
#                 70–79 vs 20–29 (both directions) per tissue.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

library('CePa')
library('DESeq2')
library('dplyr')
library('stringr')
# library('biomaRt')  # Optional, not used
library('PCAtools')

# Working directory (active RStudio document path)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# I/O directories
indir  = '/Input/'
outdir = '/Output/'

# ------------------------------------------------------------------------------
# Load and prepare GTEx metadata
# ------------------------------------------------------------------------------

gtex.meta = read.csv(paste0(indir,'GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.csv'))
gtex.meta <- gtex.meta %>%
  mutate(
    SUBJID = gsub("-",".", SUBJID)
  )

# Sex codes:
# 1 = Male
# 2 = Female
#
# DTHHRDY - Hardy death classification:
# 0) Ventilator case
# 1) Violent & fast (<10 min)
# 2) Fast natural (<1 hr)
# 3) Intermediate (1–24 hrs)
# 4) Slow (>1 day)

# ------------------------------------------------------------------------------
# Discover raw count files and parse tissue names
# ------------------------------------------------------------------------------

files <- list.files(
  path        = paste0(indir, "raw_count_data/"),
  pattern     = '*.gct',
  full.names  = TRUE
)

gtex.file.df <- data.frame(file = files)
gtex.file.df <- gtex.file.df %>%
  mutate(
    tissue = str_extract(file, "v10_.*\\.gct") %>%  # extract 'v10_...' + '.gct'
      str_remove("^v10_") %>%                        # remove leading 'v10_'
      str_remove("\\.gct$")                          # remove trailing '.gct'
  )
colnames(gtex.file.df) = c('file_path', 'tissue')

# ------------------------------------------------------------------------------
# Create DESeq2 object list (one per tissue)
# ------------------------------------------------------------------------------

dds_gtex_list <- list()
tissue.list <- as.character(unique(gtex.file.df$tissue))

#tissue = tissue.list[1] # debugging

for (tissue in tissue.list) {
  print(tissue)
  
  # ----- Get counts -----
  countdata = read.gct(paste0(indir, 'raw_count_data/gene_reads_v10_', tissue, '.gct'))
  
  # Filter out genes with zero counts across all samples
  table(rowSums(countdata) > 0)
  countdata <- countdata[rowSums(countdata) > 0, ]
  
  # Save filtered counts
  write.csv(countdata, paste0(outdir, 'raw_count_filtered/gtex_gene_reads_v10_filtered_lowcountOnly_', tissue, ".csv"))
  
  # ----- Make metadata -----
  metadat = data.frame(sample_ID = colnames(countdata))
  metadat = metadat %>%
    rowwise() %>%
    mutate(
      SUBJID = {
        dots <- str_locate_all(sample_ID, fixed("."))[[1]]
        if (nrow(dots) >= 2) {
          str_sub(sample_ID, 1, dots[2, "start"] - 1)
        } else {
          sample_ID
        }
      }
    ) %>%
    ungroup()
  
  metadat.new = left_join(metadat, gtex.meta, by = 'SUBJID')
  metadat.new$tissue = tissue
  metadat.new$SEX = ifelse(metadat.new$SEX == '1', 'M', 'F')
  
  # Save experiment design
  write.csv(metadat.new, paste0(outdir, 'experimentdesign/ExperimentDesign_gtex_', tissue, "_250709.csv"))
  
  # ----- Build DESeq2 object -----
  sampleTable = metadat.new
  coldata <- DataFrame(sampleTable)
  rownames(coldata) = coldata$sample_ID
  col_order <- rownames(coldata)
  countdata <- countdata[, col_order]
  table(colnames(countdata) == rownames(coldata))
  
  if (tissue %in% c('ovary', 'testis')) {
    dds_tissue <- DESeqDataSetFromMatrix(countData = countdata,
                                         colData   = coldata,
                                         design    = ~ AGE)
  } else {
    dds_tissue <- DESeqDataSetFromMatrix(countData = countdata,
                                         colData   = coldata,
                                         design    = ~ AGE + SEX)
  }
  
  dds_gtex_list[[tissue]] = dds_tissue
  
  # Save normalized counts
  dds_gtex_list[[tissue]] <- estimateSizeFactors(dds_gtex_list[[tissue]])
  ncounts <- counts(dds_gtex_list[[tissue]], normalized = TRUE)
  write.csv(ncounts, file = paste0(outdir, 'normcounts/CountsNormDESeq2_gtex_', tissue, '_250709.csv'))
  
  rm(ncounts, countdata, dds_tissue)
  gc()
}

# Save all tissue DESeq2 objects
save(dds_gtex_list, file = '../common_robjects/dds_gtex_alltissues_250709.bin')

# ------------------------------------------------------------------------------
# Run differential expression (70–79 vs 20–29, both directions)
# ------------------------------------------------------------------------------

load(file = '../common_robjects/dds_gtex_alltissues_250709.bin')

padj_cutoff = 0.05
results_list_gtex <- list()
tissue_list = names(dds_gtex_list)

tissue_list = c("adipose_visceral_omentum","brain_cortex", "brain_frontal_cortex_ba9","brain_spinal_cord_cervical_c-1",
                "colon_sigmoid", "colon_transverse", "esophagus_muscularis", "heart_left_ventricle",
                "kidney_cortex", "kidney_medulla", "liver", "muscle_skeletal", "ovary","skin_not_sun_exposed_suprapubic",
                "small_intestine_terminal_ileum", "spleen", "testis","whole_blood")

for (tissue in tissue_list) {
  print(tissue)
  
  # Adjust design
  if (tissue %in% c('ovary', 'testis')) {
    design(dds_gtex_list[[tissue]]) <- ~ AGE
  } else {
    design(dds_gtex_list[[tissue]]) <- ~ AGE + SEX
  }
  
  # Run DESeq
  dds_gtex_list[[tissue]] <- DESeq(dds_gtex_list[[tissue]], fitType = 'local')
  
  # Extract results
  dds_tissue_temp <- dds_gtex_list[[tissue]]
  results_list_tissue <- list()
  
  start_time <- Sys.time()
  
  # ----- 70–79 vs 20–29 -----
  cond1 <- "70-79"
  cond2 <- "20-29"
  folder <- paste(cond1, cond2, sep = "_vs_")
  print(folder)
  aspect_to_test <- 'AGE'
  
  results_temp <- results(dds_tissue_temp, contrast = c(aspect_to_test, cond1, cond2), cooksCutoff = TRUE)
  results_temp$gene_symbol <- row.names(results_temp)
  resOrdered <- as.data.frame(results_temp[order(results_temp$pvalue),])
  resSig <- subset(resOrdered, padj < padj_cutoff)
  dim(resSig)
  
  results_list_tissue[[folder]]$resall     <- results_temp
  results_list_tissue[[folder]]$resOrdered <- resOrdered
  results_list_tissue[[folder]]$ressig     <- resSig
  
  # ----- 20–29 vs 70–79 -----
  cond1 <- "20-29"
  cond2 <- "70-79"
  folder <- paste(cond1, cond2, sep = "_vs_")
  print(folder)
  aspect_to_test <- 'AGE'
  
  results_temp <- results(dds_tissue_temp, contrast = c(aspect_to_test, cond1, cond2), cooksCutoff = TRUE)
  results_temp$gene_symbol <- row.names(results_temp)
  resOrdered <- as.data.frame(results_temp[order(results_temp$pvalue),])
  resSig <- subset(resOrdered, padj < padj_cutoff)
  dim(resSig)
  
  results_list_tissue[[folder]]$resall     <- results_temp
  results_list_tissue[[folder]]$resOrdered <- resOrdered
  results_list_tissue[[folder]]$ressig     <- resSig
  
  end_time <- Sys.time()
  end_time - start_time
  
  # Save after each tissue
  results_list_gtex[[tissue]]  <- results_list_tissue
  save(results_list_gtex, file='../common_robjects/dds_DESeq2_gtex_results_250701.bin')
  save(dds_gtex_list,  file='../common_robjects/dds_DESeq2_gtex_deseqobjlist_250701.bin')
  
  rm(results_list_tissue, dds_tissue_temp)
  gc()
}

# Final save
save(results_list_gtex,  file='../common_robjects/dds_DESeq2_gtex_results_250709.bin')
