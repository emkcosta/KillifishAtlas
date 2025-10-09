# Title: Bubble plots for selected GO terms (by Hypergeometric GO analysis) that show opposite directions of changes with age in males and females (sex-divergent terms)
# Author: Jingxun Chen
# Date: code compiled on 20250507
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Associated figures: Extended Data Fig. 3


# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
set.seed(1234)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path = dirname(rstudioapi::getActiveDocumentContext()$path)

# Load packages
library(dplyr)
library(ggplot2)
library(tidyr)

# Create a directory
dir.create(file.path(paste0(path,'/Output/HyperGO/Sex_divergent/')))
dir.create(file.path(paste0(path,'/Output/HyperGO/Sex_divergent/fUpSigMdown/')))
dir.create(file.path(paste0(path,'/Output/HyperGO/Sex_divergent/mUpSigFdown/')))
dir.create(file.path(paste0(path,'/Output/HyperGO/Sex_divergent/fDownSigMup/')))
dir.create(file.path(paste0(path,'/Output/HyperGO/Sex_divergent/mDownSigFup/')))
dir.create(file.path(paste0(path,'/Output/HyperGO/Sex_divergent/mDownSigFup/')))
dir.create(file.path(paste0(path,'/Output/Plots/')))


# ------------------------------------------------------------------
# Combine all the tissues into a single dataframe
# ------------------------------------------------------------------

# List of tissues
tissue = c('Bone','Brain','Fat','Gonad','Gut','Heart','Kidney','Liver','Muscle','Skin','SpinalCord','Spleen','Eye')

# test tissue
# tissue = c('Bone','Brain')
# t = 'Bone'

# -------- Run for upregulation or downregulation data for each sex separately -------

# GO for genes significantly upregulated in females and downregulated in males (fUpSigMdown)
df = data.frame() # Initialize a data frame 
for (t in tissue) {
  t_df <- read.table(paste0(path, "/Output/HyperGO/HyperGO_BP_fUpSigMdown_", t, "_250521.txt"), sep = "\t", header = TRUE, quote = "\"", fill = TRUE, stringsAsFactors = FALSE)# Final result file # read file
  t_df$tissue <- t # mark the tissue
  df <- rbind(df, t_df) # joining the data
}
fUpSigMdown <- df # save the female downregulated into a dataframe
write.csv(fUpSigMdown, './Output/HyperGO/Alltissues_SexDivergent_HyperGO_fUpSigMdown_250521.csv')


# GO for genes significantly upregulated in males and downregulated in females (mUpSigFdown)
df = data.frame() # Initialize a data frame 
for (t in tissue) {
  t_df <- read.table(paste0(path, "/Output/HyperGO/HyperGO_BP_mUpSigFdown_", t, "_250521.txt"), sep = "\t", header = TRUE, quote = "\"", fill = TRUE, stringsAsFactors = FALSE)# Final result file # read file
  t_df$tissue <- t # mark the tissue
  df <- rbind(df, t_df) # joining the data
}
mUpSigFdown <- df # save the female downregulated into a dataframe
write.csv(mUpSigFdown, './Output/HyperGO/Alltissues_SexDivergent_HyperGO_mUpSigFdown_250521.csv')


# GO for genes significantly downregulated in females and upregulated in males (fDownSigMup)
df = data.frame() # Initialize a data frame 
for (t in tissue) {
  t_df <- read.table(paste0(path, "/Output/HyperGO/HyperGO_BP_fDownSigMup_", t, "_250521.txt"), sep = "\t", header = TRUE, quote = "\"", fill = TRUE, stringsAsFactors = FALSE)# Final result file # read file
  t_df$tissue <- t # mark the tissue
  df <- rbind(df, t_df) # joining the data
}
fDownSigMup <- df # save the female downregulated into a dataframe
write.csv(fDownSigMup, './Output/HyperGO/Alltissues_SexDivergent_HyperGO_fDownSigMup_250521.csv')


# GO for genes significantly downregulated in males and upregulated in females (mDownSigFup)
df = data.frame() # Initialize a data frame 
for (t in tissue) {
  t_df <- read.table(paste0(path, "/Output/HyperGO/HyperGO_BP_mDownSigFup_", t, "_250521.txt"), sep = "\t", header = TRUE, quote = "\"", fill = TRUE, stringsAsFactors = FALSE)# Final result file # read file
  t_df$tissue <- t # mark the tissue
  df <- rbind(df, t_df) # joining the data
}
mDownSigFup <- df # save the female downregulated into a dataframe
write.csv(mDownSigFup, './Output/HyperGO/Alltissues_SexDivergent_HyperGO_mDownSigFup_250521.csv')



#------------------------------------------------------------------
# Select GO terms
# ------------------------------------------------------------------

# Reload the data if the script was closed
fUpSigMdown <- read.csv(paste0(path,'Output/HyperGO/Sex_divergent/Alltissues_SexDivergent_HyperGO_fUpSigMdown_250521.csv'), row.names = 1)
mUpSigFdown <- read.csv(paste0(path,'Output/HyperGO/Sex_divergent/Alltissues_SexDivergent_HyperGO_mUpSigFdown_250521.csv'), row.names = 1)
fDownSigMup <- read.csv(paste0(path,'Output/HyperGO/Sex_divergent/Alltissues_SexDivergent_HyperGO_fDownSigMup_250521.csv'), row.names = 1)
mDownSigFup <- read.csv(paste0(path,'Output/HyperGO/Sex_divergent/Alltissues_SexDivergent_HyperGO_mDownSigFup_250521.csv'), row.names = 1)


# Inspect each tissue to select key GO terms to plot
# Selection criteria: A GO term needs to be significant 
 

# ---------------------------- SECTION 1: fUpSigMdown --------------------------------
# These are GO terms for genes significantly upregulated in females and downregulated in males

# After inspection, here are some GO terms of interest:
# - kidney = GO:0007040, GO:0042116, GO:0007033, GO:0016236, GO:0045088, GO:0042119
# - spleen = GO:0030593, GO:0002274, GO:0050729, GO:007S0098, GO:0071621, GO:0070555
# - eye = GO:0048812, GO:0048790, GO:0007019, GO:0099172
# - gonad = GO:0050900, GO:0090382, GO:0006826, GO:0008286
# - gut = GO:0051495
# - heart = GO:0031146, GO:0070498, GO:0010972
# - liver = GO:0070661, GO:0042116, GO:0072676, GO:1990868, GO:0032609, GO:0043062
# - spinal cord = GO:0050807, GO:0006935

# Select GO term ID
IDs.byhallmarks_fUpSigMdown <- c( 'GO:0007040','GO:0007033','GO:0016236','GO:0090382',  # autophagy
                      'GO:0048812', 'GO:0048790','GO:0007019','GO:0099172', 'GO:0050807', # Neuron
                      'GO:0008286', # Insulin
                      'GO:0006826','GO:0051495',# Iron transport & cytoskeleton
                      'GO:0031146', # Proteostasis
                      'GO:0010972', 'GO:0051301', #cell divisions
                      'GO:0043062', 'GO:0022617', #ECM
                      'GO:0042116','GO:0045088','GO:0042119','GO:0030593', 'GO:0002274', # inflammations
                      'GO:0050729', 'GO:0070098', 'GO:0071621', 'GO:0070555' #inflammation
)


# Check if the selected GO terms are in the data
IDs.byhallmarks_fUpSigMdown %in% fUpSigMdown$GOBPID 

# Make a dataframe with the selected GO terms
fUpSigMdown_hallmark <- fUpSigMdown %>% filter(GOBPID %in% IDs.byhallmarks_fUpSigMdown) 

# Set order of the dataframe 
order.df <- data.frame(IDs.byhallmarks_fUpSigMdown) # set order based on hallmark
colnames(order.df) = c('GOBPID') # name the column 
rownames(order.df) = order.df$GOBPID # set row names to be the GOBPID

# Add a column called hallmark that denotes the different types of GO terms
order.df$hallmark = c('autophagy', 'autophagy','autophagy', 'autophagy',
                      'neuron','neuron','neuron','neuron','neuron',
                      'transport',
                      'transport','transport',
                      'proteostasis',
                      'cellCycle','cellCycle',
                      'ECM','ECM',
                      'inflammation', 'inflammation', 'inflammation','inflammation', 'inflammation',
                      'inflammation','inflammation', 'inflammation','inflammation')


# Add a new column in the HyperGO data to keep track of hallmark types
fUpSigMdown_hallmark$hallmark <- NA

# Loop to add the hallmark type to each GO term
for(i in 1:nrow(fUpSigMdown_hallmark)){
  row.temp = fUpSigMdown_hallmark[i,]   # temporary storage for each row
  ID = row.temp$GOBPID   # find the ID of this row
  h = order.df[ID,]$hallmark   # find the hallmark type for this ID
  fUpSigMdown_hallmark[i,]$hallmark = h   # record the hallmark type on the HyperGO data
}

# Order the dataframe based on hallmark
fUpSigMdown_hallmark <- fUpSigMdown_hallmark[order(fUpSigMdown_hallmark$hallmark), ]
fUpSigMdown_hallmark$hallmark <- factor(fUpSigMdown_hallmark$hallmark, levels = sort(unique(fUpSigMdown_hallmark$hallmark )))
fUpSigMdown_hallmark$Term <- factor(fUpSigMdown_hallmark$Term, levels = unique(fUpSigMdown_hallmark$Term))



# ---------------------------- SECTION 2: mDownSigFup --------------------------------
# These are GO terms for genes significantly downregulated in male and upregulated in females
# Similar to SECTION 1 except that the GO terms are significant in males (not in females)

# After inspection, here are some GO terms of interest:
# - Bone: GO:0006937
# - Brain: GO:0030030
# - Gut: GO:0006508, GO:0051051, GO:0062197, GO:0001525
# - Liver: GO:0043062, GO:0030198, GO:0007155, GO:0007596

# Select GO term ID
IDs.byhallmarks_mDownSigFup <- c( 'GO:0006937', # muscle
                                  'GO:0030030',# Neuron
                                  'GO:0051051', # transport 
                                  'GO:0006508', 'GO:0062197', # Proteostasis
                                  'GO:0007596', # blood coagulation
                                  'GO:0007155', 'GO:0043062', 'GO:0030198', #ECM
                                  'GO:0001525' # angiogenesis
)


# Check if the selected GO terms are in the data
IDs.byhallmarks_mDownSigFup %in% mDownSigFup$GOBPID 

# Make a dataframe with the selected GO terms
mDownSigFup_hallmark <- mDownSigFup %>% filter(GOBPID %in% IDs.byhallmarks_mDownSigFup) 

# Set order of the dataframe 
order.df <- data.frame(IDs.byhallmarks_mDownSigFup) # set order based on hallmark
colnames(order.df) = c('GOBPID') # name the column 
rownames(order.df) = order.df$GOBPID # set row names to be the GOBPID

# Add a column called hallmark that denotes the different types of GO terms
order.df$hallmark = c('muscle', 
                      'neuron',
                      'transport',
                      'proteostasis','proteostasis',
                      'coagulation',
                      'ECM','ECM','ECM',
                      'angiogenesis')


# Add a new column in the HyperGO data to keep track of hallmark types
mDownSigFup_hallmark$hallmark <- NA

# Loop to add the hallmark type to each GO term
for(i in 1:nrow(mDownSigFup_hallmark)){
  row.temp = mDownSigFup_hallmark[i,]   # temporary storage for each row
  ID = row.temp$GOBPID   # find the ID of this row
  h = order.df[ID,]$hallmark   # find the hallmark type for this ID
  mDownSigFup_hallmark[i,]$hallmark = h   # record the hallmark type on the HyperGO data
}

# Order the dataframe based on hallmark
mDownSigFup_hallmark <- mDownSigFup_hallmark[order(mDownSigFup_hallmark$hallmark), ]
mDownSigFup_hallmark$hallmark <- factor(mDownSigFup_hallmark$hallmark, levels = sort(unique(mDownSigFup_hallmark$hallmark )))
mDownSigFup_hallmark$Term <- factor(mDownSigFup_hallmark$Term, levels = unique(mDownSigFup_hallmark$Term))




# ---------------------------- SECTION 3: fDownSigMup --------------------------------
# These are GO terms for genes significantly downregulated in male and upregulated in females
# Similar to SECTION 1 except that the GO terms are significant in males (not in females)

# After inspection, here are some GO terms of interest:
# - Bone: GO:0051031, GO:0006406, GO:0000070, GO:0006260
# - Eye: GO:0042254, GO:0045047, GO:0006413, GO:0006458
# - Gonad: GO:0006261, GO:0006298
# - Gut: GO:0044283, GO:0008610
# - Kidney: GO:0048821, GO:0042743, GO:0061515, GO:0072593
# - Liver: GO:0097250
# - Spinal cord: GO:0030516, GO:0030307, GO:1990138
# - Spleen: GO:0060562

# Select GO term ID
IDs.byhallmarks_fDownSigMup <- c( 'GO:0051031','GO:0006406', # tRNA, mRNA transport
                                  'GO:0042254', 'GO:0045047', 'GO:0006413', # Ribosome and translation
                                  'GO:0044283', 'GO:0008610', # metabolism 
                                  'GO:0006458',  # Proteostasis
                                  'GO:0042743','GO:0072593', 'GO:0097250', # oxphos, mito
                                  'GO:0048821', # blood development
                                  'GO:0030516', 'GO:0030307', 'GO:1990138', # neuron
                                  'GO:0000070', 'GO:0006260', 'GO:0006298',  # cell cycle
                                  'GO:0061515',# inflammation
                                  'GO:0060562' # epithelial tube morphogenesis
)


# Check if the selected GO terms are in the data
IDs.byhallmarks_fDownSigMup %in% fDownSigMup$GOBPID 

# Make a dataframe with the selected GO terms
fDownSigMup_hallmark <- fDownSigMup %>% filter(GOBPID %in% IDs.byhallmarks_fDownSigMup) 

# Set order of the dataframe 
order.df <- data.frame(IDs.byhallmarks_fDownSigMup) # set order based on hallmark
colnames(order.df) = c('GOBPID') # name the column 
rownames(order.df) = order.df$GOBPID # set row names to be the GOBPID

# Add a column called hallmark that denotes the different types of GO terms
order.df$hallmark = c('transport', 'transport',
                      'translation','translation','translation',
                      'metabolism','metabolism',
                      'proteostasis',
                      'oxphos','oxphos','oxphos',
                      'blood development',
                      'neuron','neuron','neuron',
                      'cell cycle', 'cell cycle', 'cell cycle', 
                      'inflammation',
                      'epithelial')


# Add a new column in the HyperGO data to keep track of hallmark types
fDownSigMup_hallmark$hallmark <- NA

# Loop to add the hallmark type to each GO term
for(i in 1:nrow(fDownSigMup_hallmark)){
  row.temp = fDownSigMup_hallmark[i,]   # temporary storage for each row
  ID = row.temp$GOBPID   # find the ID of this row
  h = order.df[ID,]$hallmark   # find the hallmark type for this ID
  fDownSigMup_hallmark[i,]$hallmark = h   # record the hallmark type on the HyperGO data
}

# Order the dataframe based on hallmark
fDownSigMup_hallmark <- fDownSigMup_hallmark[order(fDownSigMup_hallmark$hallmark), ]
fDownSigMup_hallmark$hallmark <- factor(fDownSigMup_hallmark$hallmark, levels = sort(unique(fDownSigMup_hallmark$hallmark )))
fDownSigMup_hallmark$Term <- factor(fDownSigMup_hallmark$Term, levels = unique(fDownSigMup_hallmark$Term))


# ---------------------------- SECTION 4: mUpSigFdown --------------------------------
# These are GO terms for genes significantly downregulated in male and upregulated in females
# Similar to SECTION 1 except that the GO terms are significant in males (not in females)

# After inspection, here are some GO terms of interest:
# - Fat: GO:0044782
# - Gut: GO:0032787
# - Skin: GO:0033559, GO:0022900, GO:0072593, GO:0002573

# Select GO term ID
IDs.byhallmarks_mUpSigFdown <- c( 'GO:0044782', # cilium
                                  'GO:0032787', 'GO:0033559', # metabolism 
                                  'GO:0022900',  'GO:0072593',# electron transport chain
                                  'GO:0002573' # inflammation
)


# Check if the selected GO terms are in the data
IDs.byhallmarks_mUpSigFdown %in% mUpSigFdown$GOBPID 

# Make a dataframe with the selected GO terms
mUpSigFdown_hallmark <- mUpSigFdown %>% filter(GOBPID %in% IDs.byhallmarks_mUpSigFdown) 

# Set order of the dataframe 
order.df <- data.frame(IDs.byhallmarks_mUpSigFdown) # set order based on hallmark
colnames(order.df) = c('GOBPID') # name the column 
rownames(order.df) = order.df$GOBPID # set row names to be the GOBPID

# Add a column called hallmark that denotes the different types of GO terms
order.df$hallmark = c('cilium', 
                      'metabolism','metabolism',
                      'ETC','ETC',
                      'inflammation')


# Add a new column in the HyperGO data to keep track of hallmark types
mUpSigFdown_hallmark$hallmark <- NA

# Loop to add the hallmark type to each GO term
for(i in 1:nrow(mUpSigFdown_hallmark)){
  row.temp = mUpSigFdown_hallmark[i,]   # temporary storage for each row
  ID = row.temp$GOBPID   # find the ID of this row
  h = order.df[ID,]$hallmark   # find the hallmark type for this ID
  mUpSigFdown_hallmark[i,]$hallmark = h   # record the hallmark type on the HyperGO data
}

# Order the dataframe based on hallmark
mUpSigFdown_hallmark <- mUpSigFdown_hallmark[order(mUpSigFdown_hallmark$hallmark), ]
mUpSigFdown_hallmark$hallmark <- factor(mUpSigFdown_hallmark$hallmark, levels = sort(unique(mUpSigFdown_hallmark$hallmark )))
mUpSigFdown_hallmark$Term <- factor(mUpSigFdown_hallmark$Term, levels = unique(mUpSigFdown_hallmark$Term))





#------------------------------------------------------------------
# Prepare to plot
# ------------------------------------------------------------------

# # # The code below is used to set the scale of dot size by including the lowest and highest p.adjust values for male and female
# # # *** Minimum padj value = 3.350261e-13 from fUpSigMdown (lower than mDownSigFup)
# # # *** Minimum padj value = 5.632083e-07 from fDownSigMup(lower than mUpSigFdown)
# # # *** Maximum padj value = 1.0000000 from both male and female
# # # *** Maximum enrichment = 50

# Maximum p.adjust value = 1.0000000. Create a row that include this info
max <- as.data.frame(c("BP", 0,0, 0, 0, 0, "None", 0, 1, "Scale", "min_scale"))
max.t <- t(max)
# Set column names
colnames(max.t) <- colnames(fUpSigMdown_hallmark)

# --------------------------- fUpSigMdown and mDownSigFup ---------------------------
# Minimum p.adjust value = 3.350261e-13. Create a row that include this info
min_fUp <- as.data.frame(c("BP", 0,0, 0, 0, 0, "None", 50, 3.350261e-13, "Scale", "max_scale"))
min_fUp.t <- t(min_fUp)
# Set column names
colnames(min_fUp.t) <- colnames(fUpSigMdown_hallmark) # set the column names to be the same as the data

# Create a data frame called scale for the min and max values
scale_fUp <- as.data.frame(rbind(min_fUp.t, max.t))

# Set the enrichment and p.adjust terms to numeric for plotting
scale_fUp$enrichment <- as.numeric(scale_fUp$enrichment)
scale_fUp$padj <- as.numeric(scale_fUp$padj)

# Combine the scale and HyperGO data
fUpSigMdown_hallmark_scale <- rbind(fUpSigMdown_hallmark, scale_fUp)
mDownSigFup_hallmark_scale <- rbind(mDownSigFup_hallmark, scale_fUp)


# --------------------------- fDownSigMup and mUpSigFdown ---------------------------
# Minimum p.adjust value = 5.632083e-07. Create a row that include this info
min_fDown <- as.data.frame(c("BP", 0,0, 0, 0, 0, "None", 50, 5.632083e-07, "Scale", "max_scale"))
min_fDown.t <- t(min_fDown)
# Set column names
colnames(min_fDown.t) <- colnames(fDownSigMup_hallmark) # set the column names to be the same as the data

# Create a data frame called scale for the min and max values
scale_fDown <- as.data.frame(rbind(min_fDown.t, max.t))

# Set the enrichment and p.adjust terms to numeric for plotting
scale_fDown$enrichment <- as.numeric(scale_fDown$enrichment)
scale_fDown$padj <- as.numeric(scale_fDown$padj)

# Combine the scale and HyperGO data
fDownSigMup_hallmark_scale <- rbind(fDownSigMup_hallmark, scale_fDown)
mUpSigFdown_hallmark_scale <- rbind(mUpSigFdown_hallmark, scale_fDown)




# ------------------------------------------------------------------
# Plot the bubble plot and save as PDF
# ------------------------------------------------------------------

# Results for fUpSigMdown
pdf(paste0(path, '/Output/Plots/Alltissue_HyperGO_BP_fUpSigMdown_GOterms_250522.pdf'), width = 15, height = 10)
ggplot(data=fUpSigMdown_hallmark_scale, aes(x=tissue, y=Term)) +
  geom_point(aes(fill=pmin(enrichment, 50), size = -log10(padj)), shape = 21, color = "black", stroke = 0.5) +
  scale_size(range = c(0, 5)) +
  scale_fill_gradient2(high = "#652A0E", low = "white", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar") +
  labs(title = 'HyperGO (BP) for genes significantly upregulated with age in females and downregulated in males', x='Tissues', y = '',size='-log(FDR)', color='Enrichment')
dev.off()

# Results for mDownSigFup
pdf(paste0(path, '/Output/Plots/Alltissue_HyperGO_BP_mDownSigFup_GOterms_250522.pdf'), width = 15, height = 10)
ggplot(data=mDownSigFup_hallmark_scale, aes(x=tissue, y=Term)) +
  geom_point(aes(fill=pmin(enrichment, 50), size = -log10(padj)), shape = 21, color = "black", stroke = 0.5) +
  scale_size(range = c(0, 5)) +
  scale_fill_gradient2(high = "#652A0E", low = "white", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar") +
  labs(title = 'HyperGO (BP) for genes significantly downregulated with age in males and upregulated in females', x='Tissues', y = '',size='-log(FDR)', color='Enrichment')
dev.off()


# Results for fDownSigMup
pdf(paste0(path, '/Output/Plots/Alltissue_HyperGO_BP_fDownSigMup_GOterms_250522.pdf'), width = 15, height = 10)
ggplot(data=fDownSigMup_hallmark_scale, aes(x=tissue, y=Term)) +
  geom_point(aes(fill=pmin(enrichment, 50), size = -log10(padj)), shape = 21, color = "black", stroke = 0.5) +
  scale_size(range = c(0, 5)) +
  scale_fill_gradient2(high = "#55C667FF", low = "white", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar") +
  labs(title = 'HyperGO (BP) for genes significantly downregulated with age in females and upregulated in males', x='Tissues', y = '',size='-log(FDR)', color='Enrichment')
dev.off()

# Results for mUpSigFdown
pdf(paste0(path, '/Output/Plots/Alltissue_HyperGO_BP_mUpSigFdown_GOterms_250522.pdf'), width = 15, height = 10)
ggplot(data=mUpSigFdown_hallmark_scale, aes(x=tissue, y=Term)) +
  geom_point(aes(fill=pmin(enrichment, 50), size = -log10(padj)), shape = 21, color = "black", stroke = 0.5) +
  scale_size(range = c(0, 5)) +
  scale_fill_gradient2(high = "#55C667FF", low = "white", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar") +
  labs(title = 'HyperGO (BP) for genes significantly upregulated with age in males and downregulated in females', x='Tissues', y = '',size='-log(FDR)', color='Enrichment')
dev.off()





# Note: In Extended Data 3, the male column and female column are placed side by side for each tissue. 
# Manual editing in Illustrator (by moving the female columns to proper location) is used to generate this side-by-side configuration. 


rm(list = ls())

# ------------------------------------------------------------------
sessionInfo() 
