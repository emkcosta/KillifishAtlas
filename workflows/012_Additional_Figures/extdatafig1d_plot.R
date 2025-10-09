# ------------------------------------------------------------------------------
# Title:          Fish Size Distribution
# Author:         Emma Costa
# Date:           Code finalized on 20250807
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description:    Reads cohort metadata, harmonizes animal IDs, merges fish size
#                 measurements into the experiment design, saves an updated metadata
#                 table, and plots size distributions by sex across age bins.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Set up
# ------------------------------------------------------------------------------

# Libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)

# Working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ------------------------------------------------------------------------------
# Load cohort metadata (Cohort 1 & Cohort 2)
# ------------------------------------------------------------------------------

coh1_metadata <- read_excel("/Input/SupFile_1_AtlasExperiment_Metadata_250415.xlsx", sheet = "Cohort1")
coh2_metadata <- read_excel("/Input/SupFile_1_AtlasExperiment_Metadata_250415.xlsx", sheet = "Cohort2")

# ------------------------------------------------------------------------------
# Harmonize ID fields and select relevant columns
# ------------------------------------------------------------------------------

# Standardize column name for Cohort 2 and remove underscores in IDs
colnames(coh2_metadata)[1] = 'animalID'
coh2_metadata <- coh2_metadata %>%
  mutate(animalID = str_replace_all(animalID, "_", ""))

# Keep only ID and size columns in each cohort table
coh1_metadata <- coh1_metadata[c('Animal ID', 'Size (cm)')]
coh2_metadata <- coh2_metadata[c('animalID', 'Size (cm)')]

# Unify column names and combine cohorts
colnames(coh1_metadata)[1] = 'animalID'
coh_comb_metadata <- rbind(coh1_metadata,coh2_metadata)
colnames(coh_comb_metadata) = c('animalID', 'size')

# ------------------------------------------------------------------------------
# Load TPM matrix and experiment design
# ------------------------------------------------------------------------------

tpm <- read.csv('../Input/TPM_Atlas_allbatches_merged_v3.csv', stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, row.names = 1)
metadata <- read.csv(file="../Input/ExperimentDesign_allbatches_combined_v7.csv", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, row.names = 1 )

# ------------------------------------------------------------------------------
# Clean / harmonize animal IDs in experiment design
# ------------------------------------------------------------------------------

# For cohort 2, remove underscores in IDs
metadata <- metadata %>%
  mutate(animalID = if_else(cohort == "2",
                            str_replace_all(animalID, "_", ""),
                            animalID))

# Fix specific ID formatting inconsistencies
metadata$animalID[metadata$animalID == "C1"] <- "C01"
metadata$animalID[metadata$animalID == "J9"] <- "J09"

# ------------------------------------------------------------------------------
# Merge size information into experiment design and save
# ------------------------------------------------------------------------------

metadata$size <- coh_comb_metadata$size[match(metadata$animalID, coh_comb_metadata$animalID)]

write.csv(metadata, file = paste0("../Output/ExperimentDesign_allbatches_combined_v7_withSize.csv"))

# ------------------------------------------------------------------------------
# Plot size distribution by sex across age bins
# ------------------------------------------------------------------------------

# Reload saved metadata (ensures we plot from the persisted file)
metadata <- read.csv(file = paste0("../Output/ExperimentDesign_allbatches_combined_v7_withSize.csv"), row.names = 1)

# Keep relevant columns, drop NAs, and deduplicate records
metadata.subset <- metadata[c('animalID', 'sex', 'age_days','size')]
metadata.subset <- na.omit(metadata.subset)
metadata.subset <- unique(metadata.subset)

# Define discrete age bins from age_days
metadata.subset$age_bin <- ifelse(metadata.subset$age_days %in% c('47', '49', '52'), '1', NA)
metadata.subset$age_bin <- ifelse(metadata.subset$age_days %in% c('75', '77', '78'), '2', metadata.subset$age_bin)
metadata.subset$age_bin <- ifelse(metadata.subset$age_days %in% c('102', '103'), '3', metadata.subset$age_bin)
metadata.subset$age_bin <- ifelse(metadata.subset$age_days %in% c('133', '134'), '4', metadata.subset$age_bin)
metadata.subset$age_bin <- ifelse(metadata.subset$age_days %in% c('133', '134'), '4', metadata.subset$age_bin)
metadata.subset$age_bin <- ifelse(metadata.subset$age_days %in% c('147', '152', '155'), '5', metadata.subset$age_bin)
metadata.subset$age_bin <- ifelse(metadata.subset$age_days %in% c('161', '162'), '6', metadata.subset$age_bin)

# Ensure size is numeric for plotting
metadata.subset$size <- as.numeric(metadata.subset$size)

# Create PDF: size distribution by sex across age bins
pdf(file = paste0("../Output/extdatafig1d_sexdist_bysex_acrossagebin.pdf"), width = 3, height = 4)
metadata.subset %>%
  ggplot(aes(x = sex, y = size, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ age_bin, nrow = 1) +
  theme_classic() +
  labs(title = "Size distribution by sex across age bins",
       x = "Sex",
       y = "Length (cm)") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(breaks = 0:6, limits = c(0, 6))
dev.off()
