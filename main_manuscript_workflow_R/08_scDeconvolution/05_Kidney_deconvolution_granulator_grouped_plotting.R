# ------------------------------------------------
# Title: Plot age-associated changes in grouped cell type proportions (granulator output)
# Author: Emma Costa
# Date:           Code compiled on 20250811
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description: Deconvolutes bulk kidney RNA-seq data using granulator output (grouped cell types) 
#              and visualizes age-associated trends by sex. 

# ------------------------------------------------


#------------------------------------------------
# Set up
#------------------------------------------------
library(granulator)

# Deconvolute Kidney bulk data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
indir = '/Input/'
outdir = '/Output/'

# load object containing bulk data
load('../common_robjects/dds_TPS_allsamples_Gonadcombo_240714.bin')

# load the granulator data results
sexes <- c("male", "female")
shortnames <- c("M", "F")
colsex <- c("M" = "#068ec9", "F" = "#ba1e2d")

#------------------------------------------------
# Plotting
#------------------------------------------------
for (i in seq_along(sexes)) {
  sex <- sexes[i]
  shortname <- shortnames[i]
  
  # Load deconvolution results
  load(paste0(outdir, "granulator-res_02-grouped_", sex, ".RData"))
  
  # Get all model names
  model_names <- names(decon.res$proportions)
  
  for (model in model_names) {
    # Extract model results by name
    model_output <- decon.res$proportions[[model]]
    
    # Combine with bulk metadata
    coldata <- as.data.frame(dds_TPS_list[['Kidney']]@colData)
    df <- cbind(coldata, model_output)
    
    # Truncate negative values
    df <- df %>%
      mutate(across(all_of(colnames(df)[21:27]), ~ ifelse(. < 0, 0, .)))
    
    # Normalize rows so proportions sum to 100%
    total <- rowSums(df[, 21:27], na.rm = TRUE)
    df[, 21:27] <- sweep(df[, 21:27], 1, total, FUN = "/") * 100
    
    # Pivot to long format
    df_long <- df %>%
      pivot_longer(
        cols = colnames(df)[21:27],  # update if needed
        names_to = "cell_type",
        values_to = "proportion"
      )
    
    # Add age bins
    df_long$age_bin <- NA
    df_long$age_bin <- ifelse(df_long$age_days %in% c(47,49,52), 1, df_long$age_bin)
    df_long$age_bin <- ifelse(df_long$age_days %in% c(75,77,78), 2, df_long$age_bin)
    df_long$age_bin <- ifelse(df_long$age_days %in% c(102,103), 3, df_long$age_bin)
    df_long$age_bin <- ifelse(df_long$age_days %in% c(133,134), 4, df_long$age_bin)
    df_long$age_bin <- ifelse(df_long$age_days %in% c(147,152,155), 5, df_long$age_bin)
    df_long$age_bin <- ifelse(df_long$age_days %in% c(161,162), 6, df_long$age_bin)
    
    # Subset to specific groups and sex
    subset_1 <- c('lymph.progenitors', 'myeloid.progenitors', 'lymphoid.cells', 'myeloid.cells')
    df_long_subset <- df_long %>%
      subset(cell_type %in% subset_1 & sex == shortname)
    
    # Plot
    p <- ggplot(df_long_subset, aes(x = age_bin, y = proportion)) +
      geom_point(alpha = 0.6, aes(colour = sex), size = 0.5) +
      scale_color_manual(values = colsex) +
      facet_grid(cell_type ~ sex) +
      theme_classic() +
      labs(
        x = "Age",
        y = "Estimated Proportion"
      ) +
      geom_smooth(method = "loess", se = T, color = "black", size = 0.2)+
      ylim(-5, 100) +
      theme(legend.position = "none")
    
    pdf(file = paste0(outdir, 'plots/250624_scdecon-cellpropbyage_grouped_', model,'_',sex, '.pdf'), width = 1.5, height = 2.5)
    print(p)
    dev.off()
  }
}
