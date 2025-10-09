#Script adapter from Rahul Nagvekar & public vignette
#Vignette: https://bioconductor.org/packages/devel/bioc/vignettes/granulator/inst/doc/granulator.html#ref-Hunt2018


library(granulator)
library(Seurat)
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tibble)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
indir = '/Input/'
outdir = '/Output/'

dir.create('./Output/plots/')
dir.create('./Output/tables/')


#Sex: male, female 
colsex <-c('#068ec9','#ba1e2d')

#------------------------------------------------
# Load in data
#------------------------------------------------

# load object containing bulk data
load('../common_robjects/dds_TPS_allsamples_Gonadcombo_240714.bin')

#load object containing single cell data 
load('../common_robjects/kidneyOnly_Teefy_SeurObj_Fordeconvolution.Rdata')
onlykidney <- SetIdent(onlykidney, value = 'Annotation_v1')
#DimPlot(onlykidney)

#------------------------------------------------
##### Expression matrix - whole kidney - best to be TPM
#------------------------------------------------
# format bk.dat
dds.tissue = dds_TPS_list[['Kidney']]
bk.dat = counts(dds.tissue, normalized = F)
bk.dat_df = as.data.frame(bk.dat)
bk.dat_mat = as.matrix(bk.dat_df)

# load gene lengths
gene.lengths = read.csv(paste0(indir, 'GeneLengths_Atlas_Plate02_Lane03-4_BatchF004_240417.csv'))
row.order = rownames(bk.dat)
gene.lengths = gene.lengths[row.order,]
len = as.vector(gene.lengths[,1])

# make TPM normalized bulk expression matrix
expr_mat = get_TPM(bk.dat_mat, len)

#------------------------------------------------
##### Signature matrix -  should be TPM - choose your input
#------------------------------------------------
avg.expr = read.csv(file = '250612_scKidney_maleOnly_TPM_avg_granulator.csv', row.names = 1) #fine celltypes, male
#avg.expr = read.csv(file = '250612_scKidney_femaleOnly_TPM_avg_granulator.csv', row.names = 1) #fine celltypes, female
#avg.expr = read.csv(file = '250624_scKidney_maleOnly_TPM_avg_granulator-grouped.csv', row.names = 1) #grouped myeloid, lymphoid, male
#avg.expr = read.csv(file = '250624_scKidney_femaleOnly_TPM_avg_granulator-grouped.csv', row.names = 1) #grouped myeloid, lymphoid, female


#------------------------------------------------
# Subset to shared genes
#------------------------------------------------
common_genes <- intersect(rownames(expr_mat), rownames(avg.expr)) #subsetting to shared genes may affect proportion outcomes
expr_mat_common <- subset(expr_mat, rownames(expr_mat) %in% common_genes)
avg.expr_common <- subset(avg.expr, rownames(avg.expr) %in% common_genes)
avg.expr_common.mat <- as.matrix(avg.expr_common)

#------------------------------------------------
# Run granulator
#------------------------------------------------
decon.res <- deconvolute(expr_mat_common, avg.expr_common.mat, use_cores = 4) #takes a long time to run ~ 1-1.5 hours
save(decon.res, file = paste0(outdir,"granulator-res_02_male.RData")) # this is the version scData TPM, male
#save(decon.res, file = paste0(outdir,"granulator-res_02_female.RData")) # this is the version scData TPM, female
#save(decon.res, file = paste0(outdir,"granulator-res_02-grouped_male.RData")) # this is the version scData TPM, male, grouped cell types
#save(decon.res, file = paste0(outdir,"granulator-res_02-grouped_female.RData")) # this is the version scData TPM, female, grouped cell types

#------------------------------------------------
# Male - Plot expected proportions from Teefy et al.
#------------------------------------------------
kidney_male <- onlykidney@meta.data %>% subset(Sex == 'M')
kidney_male_bysample <- table(kidney_male$Batch, kidney_male$Annotation_v1)
kidney_male_bysample <- as.data.frame(kidney_male_bysample)
colnames(kidney_male_bysample) <- c('sample_id', 'cell_type', 'count')

df_pct_male <- kidney_male_bysample %>%
  group_by(sample_id) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ungroup() %>%
  select(sample_id, cell_type, percent) %>%
  pivot_wider(names_from = cell_type, values_from = percent, values_fill = 0)

df_pct_male <- as.data.frame(df_pct_male)

df_pct_long <- df_pct_male %>%
  pivot_longer(
    cols = -sample_id,                  # everything except sample_id
    names_to = "cell_type",
    values_to = "percent"
  )

mean_order <- df_pct_long %>%
  group_by(cell_type) %>%
  summarise(mean_percent = mean(percent, na.rm = TRUE)) %>%
  arrange(-mean_percent) %>%
  pull(cell_type)

df_pct_long$cell_type <- factor(df_pct_long$cell_type, levels = mean_order)

# plot expected proportions
pdf(file = paste0(outdir,'plots/250623_maleOnly_celltype-proportions.pdf'))
ggplot(df_pct_long, aes(x = cell_type, y = percent)) +
  geom_bar(stat = "summary", fun = "mean", fill = "lightgray", width = 0.6) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "#068ec9") +
  theme_classic() +
  labs(
    x = "Cell Type (ordered by mean %)",
    y = "Percent of Total Cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()
dev.off()


#------------------------------------------------
# Female - Plot expected proportions from Teefy et al.
#------------------------------------------------
kidney_female <- onlykidney@meta.data %>% subset(Sex == 'F')
kidney_female_bysample <- table(kidney_female$Batch, kidney_female$Annotation_v1)
kidney_female_bysample <- as.data.frame(kidney_female_bysample)
colnames(kidney_female_bysample) <- c('sample_id', 'cell_type', 'count')


df_pct_female <- kidney_female_bysample %>%
  group_by(sample_id) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ungroup() %>%
  select(sample_id, cell_type, percent) %>%
  pivot_wider(names_from = cell_type, values_from = percent, values_fill = 0)

df_pct_female <- as.data.frame(df_pct_female)

df_pct_long <- df_pct_female %>%
  pivot_longer(
    cols = -sample_id,                  # everything except sample_id
    names_to = "cell_type",
    values_to = "percent"
  )

mean_order <- df_pct_long %>%
  group_by(cell_type) %>%
  summarise(mean_percent = mean(percent, na.rm = TRUE)) %>%
  arrange(-mean_percent) %>%
  pull(cell_type)

df_pct_long$cell_type <- factor(df_pct_long$cell_type, levels = mean_order)

# plot expected proportions
pdf(file = paste0(outdir,'plots/250623_femaleOnly_celltype-proportions.pdf'))
ggplot(df_pct_long, aes(x = cell_type, y = percent)) +
  geom_bar(stat = "summary", fun = "mean", fill = "lightgray", width = 0.6) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "#ba1e2d") +
  theme_classic() +
  labs(
    x = "Cell Type (ordered by mean %)",
    y = "Percent of Total Cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()
dev.off()
