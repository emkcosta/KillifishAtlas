# Title: Perform pathway analysis using Hypergeometric GO enrichment (focusing on sex-shared terms)
# Author: Jingxun Chen
# Date: code compiled on 20250520
# Related publication: Emma K. Costa, and Jingxun Chen et al., Nature Aging, 2025
# Associated figures: Extended Data Fig. 3, Extended Data Fig. 4

# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------

# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path = dirname(rstudioapi::getActiveDocumentContext()$path)

# Load packages
## This is to do GO enrichment analysis based on BH list from zebrafish/human
# *** IMPORTANT: Run it one by one for each of BP, MF and CC in the function #####
library("GOstats")
library('dplyr')
library("GSEABase")

# Select ontology input: BP, MF, CC (only ran for BP)
ontolg = "BP"

# Select one sex at a time to run
sex = 'female'
# sex = 'male'

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
  data <- read.csv(paste0(path, "/Input/List_by_tissue/241208_allCorrresults", t ,"_", sex,"Only.csv"))
  
  # --------------------- Process the data ------------------------
  # Find the age-correlated genes, then categorize them based on up-regulation or down-regulation 
  
  # Select either the upregulated or downregulated genes to run:
  # data_up <- subset(data, cor_spear > 0.5) # upregulated with age
  data_down <- subset(data, cor_spear < -0.5) # downregulated with age
  
  # ---------------------- Select Output files ----------------------
  # Output file name given the specific data set, matching whether this is for upregulated or downregulated
  # outfilename = paste0(path, "/Output/HyperGO/HyperGO_",ontolg,"_AgeUp_", t,"_",sex,"_250504.txt") 
  outfilename = paste0(path, "/Output/HyperGO/HyperGO_",ontolg,"_AgeDown_", t,"_",sex,"_250520.txt")
  
  # Output file name to record the genes associated with each GO term.
  # gotermlist = paste0(path, "/Output/HyperGO/GeneHyperGO",ontolg,"_AgeUp_", t,"_",sex,"_250504.txt") 
  gotermlist = paste0(path, "/Output/HyperGO/GeneHyperGO_",ontolg,"_AgeDown_", t,"_",sex,"_250520.txt")
  
  
  
  # --------------------------------------------- 
  # Define background and genes of interest
  # ---------------------------------------------
  # List of universe genes (removed padj = NA). Background.
  universe <- subset(data, !is.na(data$padj_spear))
  universe <- subset(universe, select = 'gene')
  colnames(universe) <- 'id'
  
  # ---------------------- Select Gene list to be tested -----------------------
  # Select one of the 4 options, given the specific data set and the condition of interest (up-regulation or down-regulation)  
  
  # Pick one of these options:
  # genes <- subset(data_up, select = 'gene') # upregulated
  genes <- subset(data_down, select = 'gene') # downregulated
  
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

sessionInfo() 