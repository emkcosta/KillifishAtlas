# Title: Bubble plots for selected GO terms (by Hypergeometric GO analysis) that show the same changes with age in males and females (sex-shared terms)
# Author: Jingxun Chen
# Date: code compiled on 2025020
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Associated figures: Extended Data Fig. 2


# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
set.seed(1234)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path = dirname(rstudioapi::getActiveDocumentContext()$path)
outpath = "./Output/"

# Load packages
library(dplyr)
library(ggplot2)
library(tidyr)

# Create a directory
dir.create(file.path(paste0(path,'/Output/HyperGO/F_down/')))
dir.create(file.path(paste0(path,'/Output/HyperGO/M_down/')))
dir.create(file.path(paste0(path,'/Output/HyperGO/F_up/')))
dir.create(file.path(paste0(path,'/Output/HyperGO/M_up/')))
dir.create(file.path(paste0(path,'/Output/HyperGO/Sex_shared/')))


# ------------------------------------------------------------------
# Combine all the tissues into a single dataframe
# ------------------------------------------------------------------

# List of tissues
tissue = c('Bone','Brain','Fat','Gonad','Gut','Heart','Kidney','Liver','Muscle','Skin','SpinalCord','Spleen','Eye')

# test tissue
# tissue = c('Bone','Brain')
# t = 'Bone'

# -------- Run for upregulation or downregulation data for each sex separately -------

# female downregulated
sex = 'female'
df = data.frame() # Initialize a data frame 
for (t in tissue) {
  t_df <- read.table(paste0(outpath, "Output/HyperGO/F_down/HyperGO_BP_AgeDown_", t, "_", sex, "_250520.txt"), sep = "\t", header = TRUE, quote = "\"", fill = TRUE, stringsAsFactors = FALSE)# Final result file # read file
  t_df$sex <- sex # mark the sex 
  t_df$tissue <- t # mark the tissue
  df <- rbind(df, t_df) # joining the data
}
df_f_down <- df # save the female downregulated into a dataframe

# male downregulated
sex = 'male'
df = data.frame() # Initialize a data frame 
for (t in tissue) {
  t_df <- read.table(paste0(outpath, "Output/HyperGO/M_down/HyperGO_BP_AgeDown_", t, "_", sex, "_250520.txt"), sep = "\t", header = TRUE, quote = "\"", fill = TRUE, stringsAsFactors = FALSE)# Final result file # read file
  t_df$sex <- sex # mark the sex 
  t_df$tissue <- t # mark the tissue
  df <- rbind(df, t_df) # joining the data
}
df_m_down <- df # save the male downregulated into a dataframe

# female upregulated
sex = 'female'
df = data.frame() # Initialize a data frame 
for (t in tissue) {
  t_df <- read.table(paste0(outpath, "Output/HyperGO/F_up/HyperGO_BP_AgeUp_", t, "_", sex, "_250504.txt"), sep = "\t", header = TRUE, quote = "\"", fill = TRUE, stringsAsFactors = FALSE)# Final result file # read file
  t_df$sex <- sex # mark the sex 
  t_df$tissue <- t # mark the tissue
  df <- rbind(df, t_df) # joining the data
}
df_f_up <- df # save the female downregulated into a dataframe

# male upregulated
sex = 'male'
df = data.frame() # Initialize a data frame 
for (t in tissue) {
  t_df <- read.table(paste0(outpath, "Output/HyperGO/M_up/HyperGO_BP_AgeUp_", t, "_", sex, "_250504.txt"), sep = "\t", header = TRUE, quote = "\"", fill = TRUE, stringsAsFactors = FALSE)# Final result file # read file
  t_df$sex <- sex # mark the sex 
  t_df$tissue <- t # mark the tissue
  df <- rbind(df, t_df) # joining the data
}
df_m_up <- df # save the male downregulated into a dataframe


# Combine all tissues and sex, for downregulated or upregulated GO terms separately
df_down <- rbind(df_f_down, df_m_down) # downregulated, all tissues + sex
df_up <- rbind(df_f_up, df_m_up) # upregulated, all tissues + sex

# Export these data
write.csv(df_down, paste0(outpath, 'Output/HyperGO/Alltissues_sexsplit_corres_HyperGO_AgeDown_250520.csv'))
write.csv(df_up, paste0(outpath, 'Output/HyperGO/Alltissues_sexsplit_corres_HyperGO_AgeUp_250507.csv'))



#------------------------------------------------------------------
# Select GO terms
# ------------------------------------------------------------------

# Reload the data if the script was closed
df_down <- read.csv(paste0(outpath, 'Output/HyperGO/Sex_shared/Alltissues_sexsplit_corres_HyperGO_AgeDown_250520.csv'), row.names = 1)
df_up <- read.csv(paste0(outpath, 'Output/HyperGO/Sex_shared/Alltissues_sexsplit_corres_HyperGO_AgeUp_250507.csv'), row.names = 1)

# ------ For downregulated terms ------
# Filter to keep only significant p.adj
df_down_sig <- df_down[df_down$padj < 0.05, ]

# Select GO term ID
# downregulated
# First check whether the same terms from GSEA was present in the HyperGO list (remove those not found)
# A term needs to be present in at least 2 tissues to be retained in this full list
IDs.byhallmarks_down <- c( 'GO:0009058','GO:0006520',  # metabolism (biosynthetic process, amino acid)
                           'GO:0000723', # telomere maintenance
                           'GO:0045047', # proteostasis - protein targeting to ER
                           'GO:0006396', 'GO:0016071', 'GO:0000398', ## mRNA processing
                           'GO:0042254','GO:0006412', 'GO:0006364', ## Ribosome
                           'GO:0032543', 'GO:0042773', 'GO:0009060', #mitochondrial dysfunction (mito translation, respiration)
                           'GO:0051276','GO:0045787','GO:0007049', 'GO:0006260', ## Cell divisions and DNA replication
                           'GO:0030198', 'GO:0043062') ## ECM (only muscle shows this)


# Check if the selected GO terms are in the data
IDs.byhallmarks_down %in% df_down$GOBPID 

# Make a dataframe with the selected GO terms
df_down_hallmark <- df_down %>% filter(GOBPID %in% IDs.byhallmarks_down) 

# Set order of the dataframe 
order.df_down <- data.frame(IDs.byhallmarks_down) # set order based on hallmark
colnames(order.df_down) = c('GOBPID') # name the column 
rownames(order.df_down) = order.df_down$GOBPID # set row names to be the GOBPID

# Add a column called hallmark that denotes the different types of GO terms
order.df_down$hallmark = c('metabolism', 'metabolism',
                      'telomere',
                      'proteostasis',
                      'mRNAprocess', 'mRNAprocess', 'mRNAprocess',
                      'ribosome','ribosome','ribosome',
                      'mitochondria', 'mitochondria','mitochondria',
                      'cellCycle','cellCycle','cellCycle','cellCycle',
                      'ECM','ECM')

# Add a new column in the HyperGO data to keep track of hallmark types
df_down_hallmark$hallmark <- NA

# Loop to add the hallmark type to each GO term
for(i in 1:nrow(df_down_hallmark)){
  row.temp = df_down_hallmark[i,]   # temporary storage for each row
  ID = row.temp$GOBPID   # find the ID of this row
  h = order.df_down[ID,]$hallmark   # find the hallmark type for this ID
  df_down_hallmark[i,]$hallmark = h   # record the hallmark type on the HyperGO data
}

# Order the dataframe based on hallmark
df_down_hallmark <- df_down_hallmark[order(df_down_hallmark$hallmark), ]
df_down_hallmark$hallmark <- factor(df_down_hallmark$hallmark, levels = sort(unique(df_down_hallmark$hallmark )))
df_down_hallmark$Term <- factor(df_down_hallmark$Term, levels = unique(df_down_hallmark$Term))



# ------ For upregulated terms ------
# Filter to keep only significant p.adj
df_up_sig <- df_up[df_up$padj < 0.05, ]

# Select GO term ID
# Upregulated
# First check whether the same terms from GSEA was present in the HyperGO list (remove those not found)
# A term needs to be present in at least 2 tissues to be retained in this full list

IDs.byhallmarks_up <- c( 'GO:0042254','GO:0042255','GO:0002181', 'GO:0006413', 'GO:0006364',## Ribosome
                           'GO:0006914','GO:0010508','GO:0016236','GO:0007040', 'G0:0032418', #autophagy & lysosome
                           'GO:0032868','GO:0032869','GO:0071396', ## Metabolism
                           'GO:0030198', 'GO:0022617', ## ECM
                           'GO:0032635','GO:0032612','GO:0032637', 'GO:0072593', 'GO:0002250', 'GO:0001816', 'GO:0045087' #inflammation
)

# Check if the selected GO terms are in the data
IDs.byhallmarks_up %in% df_up$GOBPID 

# Make a dataframe with the selected GO terms
df_up_hallmark <- df_up %>% filter(GOBPID %in% IDs.byhallmarks_up) 

# Set order of the dataframe 
order.df_up <- data.frame(IDs.byhallmarks_up) # set order based on hallmark
colnames(order.df_up) = c('GOBPID') # name the column 
rownames(order.df_up) = order.df_up$GOBPID # set row names to be the GOBPID

# Add a column called hallmark that denotes the different types of GO terms
order.df_up$hallmark = c('ribosome', 'ribosome','ribosome','ribosome','ribosome',
                           'autophagy', 'autophagy', 'autophagy','autophagy','autophagy',
                           'metabolism', 'metabolism','metabolism',
                           'ECM','ECM', 
                         'inflammation','inflammation','inflammation','inflammation','inflammation','inflammation','inflammation')

# Add a new column in the HyperGO data to keep track of hallmark types
df_up_hallmark$hallmark <- NA

# Loop to add the hallmark type to each GO term
for(i in 1:nrow(df_up_hallmark)){
  row.temp = df_up_hallmark[i,]   # temporary storage for each row
  ID = row.temp$GOBPID   # find the ID of this row
  h = order.df_up[ID,]$hallmark   # find the hallmark type for this ID
  df_up_hallmark[i,]$hallmark = h   # record the hallmark type on the HyperGO data
}

# Order the dataframe based on hallmark
df_up_hallmark <- df_up_hallmark[order(df_up_hallmark$hallmark), ]
df_up_hallmark$hallmark <- factor(df_up_hallmark$hallmark, levels = sort(unique(df_up_hallmark$hallmark )))
df_up_hallmark$Term <- factor(df_up_hallmark$Term, levels = unique(df_up_hallmark$Term))



#------------------------------------------------------------------
# Prepare to plot
# ------------------------------------------------------------------
df_down_hallmark$Condition = paste(df_down_hallmark$sex, df_down_hallmark$tissue, sep = "_")
df_up_hallmark$Condition = paste(df_up_hallmark$sex, df_up_hallmark$tissue, sep = "_")

# # The code below is used to set the scale of dot size by including the lowest and highest p.adjust values for male and female
# # *** Minimum padj value = 1.496032e-52 from male (downregulated) 
# # *** Minimum padj value = 2.453895e-41 from female (upregulated) 
# # *** Maximum padj value = 1.0000000 from both male and female

# Maximum p.adjust value = 1.0000000. Create a row that include this info
max1 <- as.data.frame(c("BP", 0,0, 0, 0, 0, "None", 50, 1, "male", "Scale", "MinScale", "male_Scale"))
max2 <- as.data.frame(c("BP", 0,0, 0, 0, 0, "None", 50, 1, "female", "Scale", "MinScale", "female_Scale"))
max1.t <- t(max1)
max2.t <- t(max2)
# Set column names
colnames(max1.t) <- colnames(df_down_hallmark)
colnames(max2.t) <- colnames(df_down_hallmark)
max <- rbind(max1.t, max2.t)

# --------------------------- Downregulated --------------------------- 
# Minimum p.adjust value = 1.496032e-52. Create a row that include this info
min1_d <- as.data.frame(c("BP", 0,0, 0, 0, 0, "None", 0, 1.496032e-52, "male", "Scale", "MinScale", "male_Scale"))
min2_d <- as.data.frame(c("BP", 0,0, 0, 0, 0, "None", 0, 1.496032e-52, "female", "Scale", "MinScale", "female_Scale"))
min1_d.t <- t(min1_d) # transpose rows and columns
min2_d.t <- t(min2_d) # transpose rows and columns
# Set column names
colnames(min1_d.t) <- colnames(df_down_hallmark) # set the column names to be the same as the data
colnames(min2_d.t) <- colnames(df_down_hallmark) # set the column names to be the same as the data
min_d <- rbind(min1_d.t, min2_d.t)

# Create a data frame called scale for the min and max values
scale_d <- as.data.frame(rbind(min_d, max))

# Set the enrichment and p.adjust terms to numeric for plotting
scale_d$enrichment <- as.numeric(scale_d$enrichment)
scale_d$padj <- as.numeric(scale_d$padj)

# Combine the scale and HyperGO data
df_down_hallmark_scale <- rbind(df_down_hallmark, scale_d)

# --------------------------- Upregulated --------------------------- 
# Upregulated: Minimum p.adjust value = 2.453895e-41. Create a row that include this info
min1_u <- as.data.frame(c("BP", 0,0, 0, 0, 0, "None", 0, 2.453895e-41, "male", "Scale", "MinScale", "male_Scale"))
min2_u <- as.data.frame(c("BP", 0,0, 0, 0, 0, "None", 0, 2.453895e-41, "female", "Scale", "MinScale", "female_Scale"))
min1_u.t <- t(min1_u) # transpose rows and columns
min2_u.t <- t(min2_u) # transpose rows and columns
# Set column names
colnames(min1_u.t) <- colnames(df_up_hallmark) # set the column names to be the same as the data
colnames(min2_u.t) <- colnames(df_up_hallmark) # set the column names to be the same as the data
min_u <- rbind(min1_u.t, min2_u.t)

# Create a data frame called scale for the min and max values
scale_u <- as.data.frame(rbind(min_u, max))

# Set the enrichment and p.adjust terms to numeric for plotting
scale_u$enrichment <- as.numeric(scale_u$enrichment)
scale_u$padj <- as.numeric(scale_u$padj)

# Combine the scale and HyperGO data
df_up_hallmark_scale <- rbind(df_up_hallmark, scale_u)

# ------------------------------------------------------------------
# Plot the bubble plot and save as PDF
# ------------------------------------------------------------------

# Downregulated
# Results for females
pdf(paste0(outpath, '/Plots/Alltissue_HyperGO_BP_F_sharedDownGOterms_250524.pdf'), width = 15, height = 10)
ggplot(data=df_down_hallmark_scale %>% filter(sex == 'female'), aes(x=tissue, y=Term)) +
  geom_point(aes(fill=pmin(enrichment, 50), size = -log10(padj)), shape = 21, color = "black", stroke = 0.5) +
  scale_size(range = c(0, 5)) +
  scale_fill_gradient2(high = "darkblue", low = "white", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar") +
  labs(title = 'HyperGO by Aging Hallmark shared GO terms for females across tissue', x='Tissues', y = '',size='-log(FDR)', color='Enrichment')
dev.off()

# Results for males
pdf(paste0(outpath, '/Plots/Alltissue_HyperGO_BP_M_sharedDownGOterms_250524.pdf'), width = 15, height = 10)
ggplot(data=df_down_hallmark_scale %>% filter(sex == 'male'), aes(x=tissue, y=Term)) +
  geom_point(aes(fill=pmin(enrichment, 50), size = -log10(padj)), shape = 21, color = "black", stroke = 0.5) +
  scale_size(range = c(0, 5)) +
  scale_fill_gradient2(high = "darkblue", low = "white", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar") +
  labs(title = 'HyperGO by Aging Hallmark shared GO terms for males across tissue', x='Tissues', y = '',size='-log(FDR)', color='Enrichment')
dev.off()



# Upregulated
# Results for females
# Cap the enrichment score be 50 (i.e., anything higher tha 50 will be plotted as red)
pdf(paste0(outpath, '/Plots/Alltissue_HyperGO_BP_F_sharedUpGOterms_250524.pdf'), width = 15, height = 10)
ggplot(data=df_up_hallmark_scale %>% filter(sex == 'female'), aes(x=tissue, y=Term)) +
  geom_point(aes(fill=pmin(enrichment, 50), size = -log10(padj)), shape = 21, color = "black", stroke = 0.5) +
  scale_size(range = c(0, 5)) +
  scale_fill_gradient2(high = "red", low = "white", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar") +
  labs(title = 'HyperGO by Aging Hallmark shared GO terms for females across tissue', x='Tissues', y = '',size='-log(FDR)', color='Enrichment')
dev.off()

# Results for males
# Cap the enrichment score be 50 (i.e., anything higher tha 50 will be plotted as red)
pdf(paste0(outpath, '/Plots/Alltissue_HyperGO_BP_M_sharedUpGOterms_250524.pdf'), width = 15, height = 10)
ggplot(data=df_up_hallmark_scale %>% filter(sex == 'male'), aes(x=tissue, y=Term)) +
  geom_point(aes(fill=pmin(enrichment, 50), size = -log10(padj)), shape = 21, color = "black", stroke = 0.5) +
  scale_size(range = c(0, 5)) +
  scale_fill_gradient2(high = "red", low = "white", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar") +
  labs(title = 'HyperGO by Aging Hallmark shared GO terms for males across tissue', x='Tissues', y = '',size='-log(FDR)', color='Enrichment')
dev.off()




# Note: In Extended Data 3, the male column and female column are placed side by side for each tissue. 
# Manual editing in Illustrator (by moving the female columns to proper location) is used to generate this side-by-side configuration. 










rm(list = ls())

# ------------------------------------------------------------------
sessionInfo() 
