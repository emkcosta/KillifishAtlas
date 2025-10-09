rm(list=ls())
# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
library('ComplexHeatmap')
library('ggstatsplot')

# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cwd = dirname(rstudioapi::getActiveDocumentContext()$path)

#RStudio Version 2023.12.1+402 (2023.12.1+402)
#R Version 4.3.3

#Bone, Brain, Eye, Fat, Gut, Heart, Kidney, Liver, Muscle, Ovary, Skin, SpinalCord, Spleen, Testis
colTPS <- c('#d8413f', '#00a550','#b8b8c0', '#eee09b','#010101','#f0932e', '#fcd328', '#6cc0ee','#f4c489','#c9bad4' ,'#ab5673', '#f1a8a4','#ef9ac2','#93cca8' ) 

# ------------------------------------------------------------------
# Load data - edited in 000_merge_counts_TPM.R, checkreads, featurecounts
# ------------------------------------------------------------------
countdata <- read.csv(paste0(cwd,"/Output/Counts_Atlas_allbatches_merged_v3.csv"), row.names = 1)
sampleTable <- read.csv(paste0(cwd,'Output/ExperimentDesign_allbatches_combined_v5.csv'), row.names = 1)

# ------------------------------------------------------------------
# Metadata harmonization and validation
# ------------------------------------------------------------------
table(sampleTable$animalID) #it appears that some of the eye samples have slightly erroneous sampleIDs - let's standardize them
sampleTable$animalID <- gsub("C_03", "C03", as.character(sampleTable$animalID)) 
sampleTable$animalID <- gsub("D_05", "D05", as.character(sampleTable$animalID)) 
sampleTable$animalID <- gsub("I_03", "I03", as.character(sampleTable$animalID)) 
sampleTable$animalID <- gsub("J_04", "J04", as.character(sampleTable$animalID)) 
sampleTable$animalID <- gsub("-", "_", as.character(sampleTable$animalID)) #some of the names had dashes where they should have had underscores
View(table(sampleTable$animalID)) #looks good now

animal_sex <- paste0(sampleTable$animalID, "_", sampleTable$sex)
p <- table(animal_sex) #it appears that some of the eye samples have the sex erroneously labeled
p

p <- as.data.frame(p)
p <- p %>% filter(Freq == 1) %>% filter(!grepl('and',animal_sex)) %>% filter(!grepl('NA',animal_sex))
p$animalID <- sapply(strsplit(as.character(p$animal_sex), "_"), `[`, 1)
p$sex <- sapply(strsplit(as.character(p$animal_sex), "_"), `[`, 2)
p$sex_corr <- ifelse(p$sex == "M", "F", "M")
rownames(p) <- p$animalID

to_repl <- p$animalID
for(i in to_repl){
  idx = which(sampleTable$animalID == i & sampleTable$tissue == 'Eye')
  sampleTable[idx,]$sex <- p[i,]$sex_corr
}

animal_sex <- paste0(sampleTable$animalID, "_", sampleTable$sex)
View(table(animal_sex)) #looks good now



#save the edited Experiment Design file
write.csv(sampleTable, paste0(cwd,"/Output/ExperimentDesign_allbatches_combined_v6.csv"))

View(table(paste0(sampleTable$animalID, "_", sampleTable$age_days))) #looks like there is one sample from animal P_1B_1 erroneously annotated as being of age "7", instead of 77
idx = which(sampleTable$age_days == 7)
sampleTable$age_days[idx] = 77

#save the edited Experiment Design file
write.csv(sampleTable, paste0(cwd,"/Output/ExperimentDesign_allbatches_combined_v7.csv"))







