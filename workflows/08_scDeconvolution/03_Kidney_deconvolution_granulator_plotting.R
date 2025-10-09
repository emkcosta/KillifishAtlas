# ------------------------------------------------------------------------------
# Title:          Deconvoluted Cell Proportions by Age and Sex in Kidney
# Author:         Emma Costa
# Date:           Code compiled on 20250808
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Description:    Loads granulator-derived deconvolution results for kidney tissue,
#                 summarizes proportions by cell type groups, performs pairwise
#                 Wilcoxon tests across age bins, and plots trajectories.
# ------------------------------------------------------------------------------

#------------------------------------------------
# Set up
#------------------------------------------------
set.seed(123)

library(granulator)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

indir = "/Input/"
outdir = "/Output/"

# load object containing bulk data
load('../common_robjects/dds_TPS_allsamples_Gonadcombo_240714.bin')

# load the granulator data results
sexes <- c("male", "female")
shortnames <- c("M", "F")
colsex <- c("M" = "#068ec9", "F" = "#ba1e2d")
old_bin_touse <- c(6,4) #age_bin6 for male ,age_bin4 for female

#------------------------------------------------
# Plotting
#------------------------------------------------
for (i in seq_along(sexes)) {
  sex <- sexes[i]
  shortname <- shortnames[i]
  
  # Load deconvolution results
  load(paste0(outdir, "granulator-res_02_", sex, ".RData"))

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
      mutate(across(all_of(colnames(df)[21:39]), ~ ifelse(. < 0, 0, .)))
    
    # Normalize rows so proportions sum to 100%
    total <- rowSums(df[, 21:39], na.rm = TRUE)
    df[, 21:39] <- sweep(df[, 21:39], 1, total, FUN = "/") * 100
    
    # Define cell type groups
    lymph.prog <- c('B.Cell.Progenitors','Lymphoid.progenitors','NK.T.progenitor.cells')
    other.prog <- c('Erythrocyte.Progenitors','HSPCs','Multipotent.progenitors')
    mye.prog <- c('Myeloid.progenitors','Neutrophil.Progenitors')
    other <- c('Erythrocytes')
    lymphoid <- c('B.cells','NK.T.cells')
    myeloid <- c('Macrophages', 'Mast.cells','Neutrophils','Thrombocytes')
    non.immune <- c('Endothelial','Fibroblasts','Kidney.distal.tubule','Kidney.prox.tubule')
    
    # Add group sums
    df$lymph_progenitors <- rowSums(df[, lymph.prog], na.rm = TRUE)
    df$other_progenitors <- rowSums(df[, other.prog], na.rm = TRUE)
    df$myeloid_progenitors <- rowSums(df[, mye.prog], na.rm = TRUE)
    df$lymphoid_cells <- rowSums(df[, lymphoid], na.rm = TRUE)
    df$myeloid_cells <- rowSums(df[, myeloid], na.rm = TRUE)
    df$nonimmune_cells <- rowSums(df[, non.immune], na.rm = TRUE)
    
    # Add age bins
    df$age_bin <- NA
    df$age_bin <- ifelse(df$age_days %in% c(47,49,52), 1, df$age_bin)
    df$age_bin <- ifelse(df$age_days %in% c(75,77,78), 2, df$age_bin)
    df$age_bin <- ifelse(df$age_days %in% c(102,103), 3, df$age_bin)
    df$age_bin <- ifelse(df$age_days %in% c(133,134), 4, df$age_bin)
    df$age_bin <- ifelse(df$age_days %in% c(147,152,155), 5, df$age_bin)
    df$age_bin <- ifelse(df$age_days %in% c(161,162), 6, df$age_bin)
    
    df <- df %>% subset(sex == shortname & age_bin %in% c(1,old_bin_touse[i])) %>% 
      mutate(
        age_bin = factor(age_bin, levels = c(1, old_bin_touse[i]))  # ensure Young then Old
      )
    
    ##### First granular cell types plot #####
    df_longer <- df %>%
      pivot_longer(
        cols = c('B.cells','NK.T.cells','Macrophages', 'Mast.cells','Neutrophils','Thrombocytes'),  # update if needed
        names_to = "cell_type",
        values_to = "proportion"
       ) %>%
    mutate(cell_type = factor(cell_type,
                            levels = c('NK.T.cells',
                                       'B.cells',
                                       'Thrombocytes',
                                       'Neutrophils',
                                       'Mast.cells',
                                       'Macrophages')))
    
    ###### Plot ##### 
    p <- ggplot(df_longer, aes(x = as.character(age_bin), y = proportion)) +
      stat_summary(fun = mean, geom = "crossbar", 
                   width = 0.3, colour = "black", fatten = 1) +
      geom_jitter(alpha = 0.6, aes(colour = sex), size = 0.5, width = 0.3, height = 0) +
      facet_grid(cell_type ~ sex) +
      scale_color_manual(values = colsex) +
      theme_classic() +
      labs(
        x = "Age",
        y = "Estimated Proportion"
      ) +
      ylim(-3, 30) +
      scale_x_discrete(limits = c("1",as.character(old_bin_touse[i]))) +
      theme(legend.position = "none")
    
    pdf(file = paste0(outdir, 'plots/250820_scdecon-cell-granular-propbyage_', model,'_',sex, '.pdf'), width = 1.5, height = 4)
    print(p)
    dev.off()
    
    ##### Stats #####
    # Two-sample t-test per facet
    t_res <- df_longer %>%
      group_by(cell_type, sex) %>%
      do(broom::tidy(t.test(proportion ~ age_bin, data = .))) %>%
      ungroup()
    
    # Wilcoxon rank-sum per facet
    w_res <- df_longer %>%
      group_by(cell_type, sex) %>%
      do(broom::tidy(wilcox.test(proportion ~ age_bin, data = ., exact = FALSE))) %>%
      ungroup()
    
    # output the table results
    write.csv(t_res, paste0(outdir, 'tables/250820_scdecon_ttest-results_granularcelltypes_', model,'_',sex, '.csv'))
    write.csv(w_res, paste0(outdir, 'tables/250820_scdecon_wilcox-results__granularcelltypes_', model,'_',sex, '.csv'))
    
    ##### Second cells grouped plot #####
    # Pivot to long format
    df$age_bin <- as.numeric(as.character(df$age_bin))
    df_long <- df %>%
      pivot_longer(
        cols = colnames(df)[40:45],  # update if needed
        names_to = "cell_type",
        values_to = "proportion"
      )

    # Subset to specific groups and sex
    subset_1 <- c('lymphoid_cells', 'myeloid_cells')
    df_long_subset <- df_long %>%
      subset(cell_type %in% subset_1 & sex == shortname & age_bin %in% c(1,old_bin_touse[i])) %>% 
      mutate(
        age_bin = factor(age_bin, levels = c(1, old_bin_touse[i]))  # ensure Young then Old
      )
    
    ##### Stats #####
    # Two-sample t-test per facet
    t_res <- df_long_subset %>%
      group_by(cell_type, sex) %>%
      do(broom::tidy(t.test(proportion ~ age_bin, data = .))) %>%
      ungroup()
    
    # Wilcoxon rank-sum per facet
    w_res <- df_long_subset %>%
      group_by(cell_type, sex) %>%
      do(broom::tidy(wilcox.test(proportion ~ age_bin, data = ., exact = FALSE))) %>%
      ungroup()
    
    # output the table results
    write.csv(t_res, paste0(outdir, 'tables/250820_scdecon_ttest-results_', model,'_',sex, '.csv'))
    write.csv(w_res, paste0(outdir, 'tables/250820_scdecon_wilcox-results_', model,'_',sex, '.csv'))
    
    ###### Plot #####
    p <- ggplot(df_long_subset, aes(x = as.character(age_bin), y = proportion)) +
      stat_summary(fun = mean, geom = "crossbar", 
                   width = 0.3, colour = "black", fatten = 1) +
      geom_jitter(alpha = 0.6, aes(colour = sex), size = 0.5, width = 0.3, height = 0) +
      scale_color_manual(values = colsex) +
      facet_grid(cell_type ~ sex) +
      theme_classic() +
      labs(
        x = "Age",
        y = "Estimated Proportion"
      ) +
      ylim(-5, 50) +
      scale_x_discrete(limits = c("1",as.character(old_bin_touse[i]))) +
      theme(legend.position = "none")

    pdf(file = paste0(outdir, 'plots/250820_scdecon-cellpropbyage_', model,'_',sex, '.pdf'), width = 1.5, height = 2.5)
    print(p)
    dev.off()
    
    
  }
}
