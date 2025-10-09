# ------------------------------------------------------------------------------
# Title:          Visual summaries of overlapping GO BP terms (killifish vs human)
# Author:         Emma Costa
# Date:           Code compiled on 20250812
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description:    Loads per-tissue overlapping GO BP terms, computes semantic
#                 similarity reduction (rrvgo), and generates heatmap, scatter,
#                 treemap, and bubble plots. Also summarizes GO term counts
#                 across tissues.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
library('dplyr')
library('tidyr')
library('ggplot2')
library('GOSemSim')
library('clusterProfiler')
library('rrvgo')

#check out: https://yulab-smu.top/biomedical-knowledge-mining-book/useful-utilities.html

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
outdir = '/Output/'

#obtain the GO structure and GO.db
#hsGO <- godata('org.Hs.eg.db', ont="BP")

# ------------------------------------------------------------------------------
# Select tissue (uncomment only one block)
# ------------------------------------------------------------------------------
# Brain
kf.tissue = 'Brain'
hum.tissue = 'brain_cortex'

# Heart
#kf.tissue = 'Heart'
#hum.tissue = 'heart_left_ventricle'

# Fat
#kf.tissue = 'Fat'
#hum.tissue = 'adipose_visceral_omentum'

# ------------------------------------------------------------------------------
# Load overlapping GO BP terms
# ------------------------------------------------------------------------------
commonBP <- read.csv(paste0(outdir, "250712_BPcommonGOIDs_",kf.tissue,"-vs-",hum.tissue,"_atlas-vs-human.csv"), row.names = 1)

commonBP$mlog10QvalxNES.hum = -log10(commonBP$p.adjust.hum) * commonBP$NES.hum
commonBP$mlog10QvalxNES.kf = -log10(commonBP$p.adjust.kf) * commonBP$NES.kf

commonBP$summed_score = commonBP$mlog10QvalxNES.hum + commonBP$mlog10QvalxNES.kf

commonBP.simple <- commonBP[,c(1:3, 6,8,16, 18, 26)]
commonBP.simple <- commonBP.simple %>% filter(p.adjust.hum < 0.1 & p.adjust.kf < 0.1)
commonBP.simple <- commonBP.simple %>% arrange(desc(summed_score))


# filtering by significance
# brain: +ve and -ve NES remain 
# heart: +ve and -ve NES remain 
# fat: +ve and -ve NES remain 


# ------------------------------------------------------------------------------
# Similarity reduction of GO terms (rrvgo)
# ------------------------------------------------------------------------------
# calculate similarity matrix
simMatrix <- calculateSimMatrix(commonBP.simple$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

# find similarity between all go terms (the thresholding shouldnt matter since padj<0.1 already thresholded)
scores <- setNames(commonBP.simple$p.adjust.kf, commonBP.simple$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.1,
                                orgdb="org.Hs.eg.db")

#heatmap of reduced terms
pdf(paste0(outdir,'plots/', kf.tissue,'_vs_',hum.tissue,'_GSEA_GOBP_overlappingGO-heatmap_250719.pdf'), width = 20, height = 8)
print(heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6))
dev.off()

#dotplot scatter of reduced terms
pdf(paste0(outdir, 'plots/',kf.tissue,'_vs_',hum.tissue,'_GSEA_GOBP_overlappingGO-scatter_250719.pdf'), width = 20, height = 8)
print(scatterPlot(simMatrix, reducedTerms))
dev.off()

#treemap of reduced terms
pdf(paste0(outdir,'plots/', kf.tissue,'_vs_',hum.tissue,'_GSEA_GOBP_overlappingGO-treemap_250719.pdf'), width = 8, height = 8)
print(treemapPlot(reducedTerms))
dev.off()

reducedTerms <- reducedTerms %>% arrange(cluster)
write.csv(reducedTerms, paste0(outdir, kf.tissue,'_vs_',hum.tissue,'_GSEA_GOBP_overlappingGO_reducedterms_250719.csv'))

freq.table = as.data.frame(table(reducedTerms$cluster))
freq.table = freq.table %>% filter(Freq > 2)
reducedTerms.thresh <- reducedTerms %>% filter(cluster %in% freq.table$Var1) 

length(unique(reducedTerms.thresh$parentTerm)) 
# Fat - 9 parent terms remaining
# Heart - 18 parent terms remaining
# Brain - 3 ""

# ------------------------------------------------------------------------------
# Define GO terms of interest (hand-picked by tissue)
# ------------------------------------------------------------------------------
##### Brain GO terms
#go <- c(head(commonBP.simple$ID,5), tail(commonBP.simple$ID,5)) #score based selection
go <- c('GO:1990266', 'GO:0097530', 'GO:0002460', 'GO:0072678', 'GO:0030595', 'GO:0098754', 'GO:0000910', 'GO:0090329', 'GO:0019722', 'GO:0034765') #hand selection

#GO:1990266 - neutrophil migration
#GO:0097530 - granulocyte migration
#GO:0002460 - adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
#GO:0072678 - T cell migration
#GO:0030595 - leukocyte chemotaxis
#GO:0098754 - detoxification
#GO:0000910 - cytokinesis
#GO:0090329 - regulation of DNA-templated DNA replication
#GO:0019722 - calcium-mediated signaling
#GO:0034765 - regulation of monoatomic ion transmembrane transport


###### Heart GO terms
#go <- c(head(commonBP.simple$ID,5), tail(commonBP.simple$ID,5))
#go <- c('GO:0060326', 'GO:0045785', 'GO:0070374', 'GO:0009060', 'GO:0006260', 'GO:0006119', 'GO:0097746', 'GO:0022407', 'GO:0046395', 'GO:0042773')

#GO:0060326 - cell chemotaxis
#GO:0045785 - positive regulation of cell adhesion
#GO:0070374 - positive regulation of ERK1 and ERK2 cascade
#GO:0009060 - aerobic respiration
#GO:0006260 - DNA replication
#GO:0006119 - oxidative phosphorylation
#GO:0097746 - blood vessel diameter maintenance
#GO:0022407 - regulation of cell-cell adhesion
#GO:0046395 - carboxylic acid catabolic process
#GO:0042773 - ATP synthesis coupled electron transport


##### Fat GO terms
#go <- c(head(commonBP.simple$ID,5), tail(commonBP.simple$ID,5))
#go <- c('GO:0002757','GO:0002460', 'GO:0050853','GO:0034440', 'GO:0009060', 'GO:0009145', 'GO:0042254', 'GO:0034470', 'GO:0046395', 'GO:0034248')

#GO:0002757 - immune response-activating signaling pathway
#GO:0002460 - adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
#GO:0050853 - B cell receptor signaling pathway 
#GO:0034440 - lipid oxidation 
#GO:0009060 - aerobic respiration
#GO:0009145 - purine nucleoside triphosphate biosynthetic process
#GO:0042254 - ribosome biogenesis
#GO:0034470 - ncRNA processing
#GO:0046395 - carboxylic acid catabolic process
#GO:0034248 - regulation of amide metabolic process


# ------------------------------------------------------------------------------
# Prepare bubble plot data
# ------------------------------------------------------------------------------
# Keep only the terms of interest
data_go <- filter(commonBP.simple, commonBP.simple$ID %in% go)

colnames(data_go)[4:7] = c('NES.hum', 'padj.hum', 'NES.kf', 'padj.kf')

df_long <- data_go %>%
  pivot_longer(
    cols = c(NES.hum, padj.hum, NES.kf, padj.kf),
    names_to = c(".value", "species"),
    names_sep = "\\."
  )

df_long$Description = factor(df_long$Description, levels = unique(df_long$Description[order(desc(df_long$summed_score))]))
df_long$species = factor(df_long$species, levels = c('kf', 'hum'))

# ------------------------------------------------------------------------------
# Generate and save bubble plot
# ------------------------------------------------------------------------------

#for hand-selected
pdf(paste0(outdir, 'plots/',kf.tissue,'_vs_',hum.tissue,'_GSEA_GOBP_selectoverlappingGOterms-handselect_250719.pdf'), width = 15, height = 10)
ggplot(data=df_long, aes(x=species, y=Description)) +
  geom_point(aes(color=NES, size = -log10(padj))) +
  scale_colour_gradient2(low = "blue", mid = "purple", high = "red", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar", aesthetics = "colour") 
dev.off()

# ------------------------------------------------------------------
# Summative analysis (counts of overlapping BP terms)
# ------------------------------------------------------------------
# how many common GO terms that meet X criteria
human.GO.res.dir = paste0(outdir,'tables/')
human_files <- list.files(path = human.GO.res.dir, full.names = TRUE)
human_files <- human_files[grepl("250712_BPcommonGOIDs_", human_files, ignore.case = TRUE)]
human_files <- human_files[-13]


human_data_list <- lapply(human_files, function(file) {
  df <- read.csv(file, row.names = 1)
  df$source_file <- basename(file)  # just filename, or use file for full path
  return(df)
})

human_df <- bind_rows(human_data_list)

human_df$tissue <- sub(
  pattern = "250712_BPcommonGOIDs_(.*)_atlas-vs-human.*",
  replacement = "\\1",
  x = human_df$source_file
)

to_plot1 <- as.data.frame(table(human_df$tissue))
to_plot1$Var1 <- factor(to_plot1$Var1, levels = to_plot1$Var1[order(to_plot1$Freq)])


sum1 = ggplot(to_plot1, aes(x=Var1, y=Freq)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_y_continuous(breaks = seq(0, 6500, by = 250))

pdf(paste0(outdir,'plots/250715_atlas-vs-gtex_barplot_allBPcommonGOIDs_bytissue.pdf'))
print(sum1)
dev.off()


human_df.filter = human_df %>% filter(p.adjust.hum < 0.1 & p.adjust.kf < 0.1)
to_plot <- as.data.frame(table(human_df.filter$tissue))

zeros = setdiff(to_plot1$Var1, to_plot$Var1)
to_add = data.frame(Var1 = zeros, Freq = rep(0,length(zeros)))
to_plot = rbind(to_plot,to_add)
to_plot$Var1 <- factor(to_plot$Var1, levels = to_plot$Var1[order(to_plot$Freq)])


sum2 = ggplot(to_plot, aes(x=Var1, y=Freq)) + 
  geom_bar(stat = "identity") + 
  scale_y_continuous(breaks = seq(0, 300, by = 50), limits = c(0,300)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


pdf(paste0(outdir,'plots/250715_atlas-vs-gtex_barplot_pval-less-pt1_BPcommonGOIDs_bytissue.pdf'))
print(sum2)
dev.off()


