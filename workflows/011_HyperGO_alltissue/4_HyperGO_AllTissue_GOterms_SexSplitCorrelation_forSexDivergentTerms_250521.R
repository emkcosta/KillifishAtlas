# Title: Perform pathway analysis using Hypergeometric GO enrichment (focusing on sex-divergent terms)
# Author: Jingxun Chen
# Date: code compiled on 20250521
# Related publication: Emma K. Costa, and Jingxun Chen et al., Nature Aging, 2025
# Associated figures: Extended Data Fig. 3, Extended Data Fig. 4

# ------------------------------------------------------------------
# README 
# ------------------------------------------------------------------
# To run this script, you'd need to select one of the four categories of genes to run at a time
# For each run, select the gene category under the section called "Select Gene list to be tested"
# You'd also need to select the corresponding output files under the sections called "Select Output Files"
# The script would then loop through all the tissues in each run

# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path = dirname(rstudioapi::getActiveDocumentContext()$path)

# Create a directory
dir.create(file.path(paste0(path,'/Output/HyperGO/')))

# Load packages
## This is to do GO enrichment analysis based on BH list from zebrafish/human
# *** IMPORTANT: Run it one by one for each of BP, MF and CC in the function #####
library("GOstats")
library('dplyr')
library("GSEABase")

# Select ontology input: BP, MF, CC (only ran for BP)
ontolg = "BP"

# List of tissues
tissue = c('Bone','Brain','Fat','Gonad','Gut','Heart','Kidney','Liver','Muscle','Skin','SpinalCord','Spleen','Eye')

# test tissue
# tissue = c('Bone','Brain')
# t = 'Bone'

for (t in tissue) {
  # ------------------------------------------------------------------
  # Load data
  # ------------------------------------------------------------------
  # The age-correlated genes were imported for each tissue and sex selected
  df_m <- read.csv(paste0(path, "/Input/List_by_tissue/241208_allCorrresults", t ,"_maleOnly.csv"))
  df_f <- read.csv(paste0(path, "/Input/List_by_tissue/241208_allCorrresults", t ,"_femaleOnly.csv"))
  
  # Clean up data frames to keep only the spearman's rank correlation and the adjusted p-values
  data_m <- subset(df_m, select = c(gene, cor_spear, padj_spear))
  data_f <- subset(df_f, select = c(gene, cor_spear, padj_spear))
  
  # Rename the columns to keep track of male vs. female data
  colnames(data_m) <- c('gene', "cor_spear.m", "padj_spear.m")
  colnames(data_f) <- c('gene', "cor_spear.f", "padj_spear.f")
  
  # --------------------- Process the data ------------------------
  # Find the age-correlated genes, then categorize them based on up-regulation or down-regulation
  # male: 
  m_up_sig <- data_m[data_m$cor_spear.m > 0.5, ]
  m_up <- data_m[data_m$cor_spear.m > 0, ]
  m_down_sig <- data_m[data_m$cor_spear.m < -0.5, ]
  m_down <- data_m[data_m$cor_spear.m < 0, ]
  # female:
  f_up_sig <- data_f[data_f$cor_spear.f > 0.5, ]
  f_up <- data_f[data_f$cor_spear.f > 0, ]
  f_down_sig <- data_f[data_f$cor_spear.f < -0.5, ]
  f_down <- data_f[data_f$cor_spear.f < 0, ]
  
  # Check the number of age-correlated genes (absolute value of spearman's correlation > 0.5) but are opposite signs in males and females:
  mf_sig_1 <- inner_join(m_up_sig, f_down_sig, by = 'gene')
  mf_sig_2 <- inner_join(m_down_sig, f_up_sig, by = 'gene')
  print(paste0(t, ': # genes upregulated with age in males and down in females = ', nrow(mf_sig_1)))
  print(paste0(t, ': # genes downregulated with age in males and up in females = ', nrow(mf_sig_2)))
  
  # Generate four categories of gene sets:
  # 1) male = spearman's correlation > 0.5, female < 0. This condition is (m_up_sig + f_down)
  m_up_sig_f_down <- inner_join(m_up_sig, f_down, by = 'gene')
  # 2) male = spearman's correlation < -0.5, female > 0. This condition is (m_down_sig + f_up)
  m_down_sig_f_up <- inner_join(m_down_sig, f_up, by = 'gene')
  # 3) female = spearman's correlation > 0.5, male < 0. This condition is (f_up_sig + m_down)
  f_up_sig_m_down <- inner_join(f_up_sig, m_down, by = 'gene')
  # 4) female = spearman's correlation < -0.5, male > 0. This condition is (f_down_sig + m_up)
  f_down_sig_m_up <- inner_join(f_down_sig, m_up, by = 'gene')
  

  # ---------------------- Select Output files ----------------------
  # Output file name given the specific data set, matching each category of genes
  outfilename = paste0(path, "/Output/HyperGO/HyperGO_",ontolg,"_mUpSigFdown_", t,"_250521.txt")
  # outfilename = paste0(path, "/Output/HyperGO/HyperGO_",ontolg,"_mDownSigFup_", t,"_250521.txt")
  # outfilename = paste0(path, "/Output/HyperGO/HyperGO_",ontolg,"_fUpSigMdown_", t,"_250521.txt")
  # outfilename = paste0(path, "/Output/HyperGO/HyperGO_",ontolg,"_fDownSigMup_", t,"_250521.txt")
  
  # Output file name to record the genes associated with each GO term.
  gotermlist = paste0(path, "/Output/HyperGO/GeneHyperGO_",ontolg,"_mUpSigFdown_", t,"_250521.txt")
  # gotermlist = paste0(path, "/Output/HyperGO/GeneHyperGO_",ontolg,"_mDownSigFup_", t,"_250521.txt")
  # gotermlist = paste0(path, "/Output/HyperGO/GeneHyperGO_",ontolg,"_fUpSigMdown_", t,"_250521.txt")
  # gotermlist = paste0(path, "/Output/HyperGO/GeneHyperGO_",ontolg,"_fDownSigMup_", t,"_250521.txt")
  
  
  # --------------------------------------------- 
  # Define background and genes of interest
  # ---------------------------------------------
  # List of universe genes (removed padj = NA). Background.
  # Note: The male and female datasets have different number of genes. I used the intersection of the two.
  data <- inner_join(data_f, data_m, by = 'gene') 
  universe <- subset(data, !is.na(data$padj_spear.m)) # drop any gene with NA in male padj
  universe <- subset(universe, !is.na(data$padj_spear.f)) # drop any gene with NA in female padj
  universe <- subset(universe, select = 'gene')
  colnames(universe) <- 'id'
  
  # ---------------------- Select Gene list to be tested -----------------------
  # Select one of the 4 options, given the specific data set and the condition of interest (up-regulation or down-regulation)  
  
  # Pick one of these options:
  genes <- subset(m_up_sig_f_down, select = 'gene') # male significantly up, female down
  # genes <- subset(m_down_sig_f_up, select = 'gene') # male significantly down, female up
  # genes <- subset(f_up_sig_m_down, select = 'gene') # female significantly up, male down
  # genes <- subset(f_down_sig_m_up, select = 'gene') # female significantly down, male up
  
  # clean up dataframe
  colnames(genes) <- 'id'
  
  
  # --------------------------------------------- 
  # Hypergeometric test for GO analysis
  # ---------------------------------------------
  
  # --------------------------------------------- INPUT FILES --------------------------------------------------------
  # Go terms list. Use either zebrafish or human. Human works a bit better due to better annotations.
  frame = read.table(file ="./Input/GO_terms/GO_killifish-human_best-hits.txt", header = T, colClasses=c(rep("factor",3)))
  # ontology MF, BP, CC was selected at the beginning of the script. Need to run each type separately.
  
  # Minimum number of genes for a term to filter
  mingenes = 5 # Bare minimum is 2. More will get more general terms, less more specific. 5-10 is a good number.
  # Relative enrichment filter
  relenrich = 0 # I generally use 0 to get all terms
  
  # ------------------------------------------------------------------------------------------------------------------
  
  # This is just to get the 3 column. I already have these so I don't need it, still.
  goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)
  
  # put your data into a GOFrame object
  goFrame=GOFrame(goframeData,organism="Human")
  head(goframeData)
  
  # cast this object to a GOAllFrame object will tap into the GO.db package and populate this object with the implicated GO2All mappings for you
  goAllFrame=GOAllFrame(goFrame)
  
  # generate geneSetCollection objects
  gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
  
  # Process the universe list
  universe = universe$id
  universe = lapply(universe, as.character)
  universe = unlist(universe)
  head(universe)
  
  # Process the gene list of interest
  genes = genes$id
  genes = lapply(genes, as.character)
  genes = unlist(genes)
  head(genes)
  
  params <- GSEAGOHyperGParams(name="My Custom GSEA based annotation parameters", 
                               geneSetCollection=gsc, 
                               geneIds = genes, 
                               universeGeneIds = universe, 
                               ontology = ontolg,
                               pvalueCutoff = 1, # 1 will get all terms, and then we can filter later
                               conditional = F, # To consider GO DAG structure or not. Doesn't affect much, but we can always try
                               testDirection = "over") # Default hyper geometric test gives both over and under enriched terms. I am specifying the direction by "over". Fo depleted terms use under.
  
  # call hyperGTest to do the test
  Over <- hyperGTest(params)
  head(summary(Over))
  
  # calculate enrichment and add it to data frame.
  # Relative enrichment factor (E-value) for a GO term = (count/size)/(size/universe)
  enrichment = (summary(Over)[5]$Count / summary(Over)[6]$Size) / (summary(Over)[6]$Size / length(universe))
  
  # create a new frame
  SummaryOver = data.frame(summary(Over), enrichment)
  head(SummaryOver)
  
  # Filter the Over variable on parameters other than P-value
  # Filter the summary of OVER with size of the term, at least 2 genes for a go term
  FilteredSummaryOver = SummaryOver[(SummaryOver$Count >= mingenes & SummaryOver$enrichment >= relenrich),]
  head(FilteredSummaryOver)
  
  # adjust p value for multiple correction
  padj = p.adjust(FilteredSummaryOver$Pvalue, "BH")
  
  # Add padj to the data frame
  FinalSummaryOver = data.frame(FilteredSummaryOver, padj)
  
  # write to a file
  write.table(FinalSummaryOver, outfilename, quote = F, row.names = F, sep = "\t") # Final result file
  
  
  # --------------------- To get the genes for each GO terms ---------------------------------------------
  
  # isolate indexes for the go terms in final results
  ind.GO <- is.element(names(Over@goDag@nodeData@data), eval(parse(text=paste("FinalSummaryOver$", "GO",ontolg,"ID", sep=''))))
  selected.GO <- Over@goDag@nodeData@data[which(ind.GO)]
  
  # get a go terms and genes in a new variable for all the terms in the results of enrichment
  goTerms <- lapply(selected.GO, 
                    function(x) x$geneIds)
  names(goTerms) <- names(Over@goDag@nodeData@data)[ind.GO]
  
  # This will create a new file that will have GO terms and gene names in each GO terms.
  # Number of GO terms or lines should be equal to the enriched go terms as in the other file generated by this script.
  # This needs to be further processed to generate the desired files.
  
  for (i in 1:length(goTerms)){
    
    test = as.data.frame(do.call(rbind, goTerms[i]))
    write.table(test, file = gotermlist, quote = F, col.names = F, append = T) # append each line - so make sure the file is empty for each run, or renamed after each run
    rm(test)
  }
  
  # ------------------------------------------------------------------
  # Clear list to run the script again  
  # ------------------------------------------------------------------
  
  # rm(list=ls()) 
}

# ------------- Below are the printed results from the script ---------- 
# "Bone: # genes upregulated with age in males and down in females = 1"
# "Bone: # genes downregulated with age in males and up in females = 0"
# "Brain: # genes upregulated with age in males and down in females = 0"
# "Brain: # genes downregulated with age in males and up in females = 0"
# "Fat: # genes upregulated with age in males and down in females = 3"
# "Fat: # genes downregulated with age in males and up in females = 0"
# "Gonad: # genes upregulated with age in males and down in females = 27"
# "Gonad: # genes downregulated with age in males and up in females = 14"
# "Gut: # genes upregulated with age in males and down in females = 0"
# "Gut: # genes downregulated with age in males and up in females = 0"
# "Heart: # genes upregulated with age in males and down in females = 0"
# "Heart: # genes downregulated with age in males and up in females = 0"
# "Kidney: # genes upregulated with age in males and down in females = 2"
# "Kidney: # genes downregulated with age in males and up in females = 3"
# "Liver: # genes upregulated with age in males and down in females = 1"
# "Liver: # genes downregulated with age in males and up in females = 3"
# "Muscle: # genes upregulated with age in males and down in females = 1"
# "Muscle: # genes downregulated with age in males and up in females = 1"
# "Skin: # genes upregulated with age in males and down in females = 9"
# "Skin: # genes downregulated with age in males and up in females = 1"
# "SpinalCord: # genes upregulated with age in males and down in females = 0"
# "SpinalCord: # genes downregulated with age in males and up in females = 0"
# "Spleen: # genes upregulated with age in males and down in females = 0"
# "Spleen: # genes downregulated with age in males and up in females = 0"
# "Eye: # genes upregulated with age in males and down in females = 6"
# "Eye: # genes downregulated with age in males and up in females = 8"

sessionInfo() 