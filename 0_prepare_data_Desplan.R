#### 

# using scRNA-seq data from Desplan lab (Neset Özel)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142787
# 
# Two datasets: 
#     original data (GSE142787_Log_normalized_average_expression.xlsx)
#     binarized expression (GSE142787_Mixture_modeling.xlsx)
#       according to https://elifesciences.org/articles/50901 binarized expression 
#       is more reliable descriptors of protein abundance in a cluster
# 
# corresponding publication from Desplan lab
# https://www.nature.com/articles/s41586-020-2879-3
#
#




library(tidyverse)
library(readxl)




prepare_data <- function() {  
  
  # create a list to hold all data
  full_list <- list()
  
  ## data preparation #### 
  # preprocessing cluster names
  clusters <- read_xlsx("41586_2020_2879_MOESM4_ESM.xlsx") %>% 
    mutate("Cluster" = as.character(`Cluster number`)) %>% 
    select(Cluster, Annotation)
  
  
  # preprocessing modelled expression values
  # load all sheets of the excel file
  modeled_expression <- "GSE142787_Mixture_modeling.xlsx" %>% 
    excel_sheets() %>% 
    set_names() %>% 
    map(read_excel, path = "GSE142787_Mixture_modeling.xlsx") %>% 
    bind_rows(.id = "Timepoint" ) %>% 
    mutate(Gene = `...1`) %>% 
    select(!c("bimodal.on_level", "bimodal.off_level", "bimodal.sd", "bimodal.pi", 
              "unimodal.level", "unimodal.sd", "diff.elpd", "diff.elpd.se", 
              "exprType", "freq_cell_ON", "freq_cell_OFF", "...1"))
  
  
  # Merge in Cluster-names
  modeled_expression <- modeled_expression %>% 
    pivot_longer(cols = !c("Gene", "Timepoint"), names_to = "Cluster") %>% 
    left_join(clusters) %>% 
    # some clusters are specific to timepoints, the following search requires 
    # them to be present in all timepoints -> get complete cases and replace missing with 0 expression
    complete(Timepoint, Annotation, Gene) %>% 
    mutate(value = if_else(is.na(value), 0, value),
           Timepoint = factor(Timepoint, levels = c("P15_MM_final", "P30_MM_final", "P40_MM_final","P50_MM_final", "P70_MM_final", "Adult_MM_final"),
                              labels = c("P15", "P30", "P40", "P50", "P70", "Adult")))
  
  
  
  # preprocessing raw expression values
  log_expression <- "GSE142787_Log_normalized_average_expression.xlsx" %>% 
    excel_sheets() %>% 
    set_names() %>% 
    map(read_excel, path = "GSE142787_Log_normalized_average_expression.xlsx") %>% 
    bind_rows(.id = "Timepoint" ) %>% 
    mutate(Gene = `...1`) %>% select(!c("...1"))
  
  
  # Merge in Cluster-names
  log_expression <- log_expression %>% 
    pivot_longer(cols = !c("Gene", "Timepoint"), names_to = "Cluster") %>% 
    left_join(clusters) %>% 
    # some clusters are specific to timepoints, the following search requires 
    # them to be present in all timepoints -> get complete cases and replace missing with 0 expression
    complete(Timepoint, Annotation, Gene) %>% 
    mutate(value = if_else(is.na(value), 0, value),
           Timepoint = factor(Timepoint, levels = c("P15_average_expression", "P30_average_expression", "P40_average_expression","P50_average_expression", "P70_average_expression", "Adult_average_expression"),
                              labels = c("P15", "P30", "P40", "P50", "P70", "Adult"))) 
  
  
  
  
  # CG_name - Gene - Chromosome - Site table from flymine.org (accessed: 2022/05/29)
  full_list$Chromosome_lookup <- read_tsv("Chromosome_lookup.tsv")
  
  # create a lookup for both gene-name and CG-name
  lookup1 <- full_list$Chromosome_lookup %>% 
    mutate(Gene2 = Gene) %>%
    select(!CG_name)
  
  lookup2 <- full_list$Chromosome_lookup %>% 
    mutate(Gene2 = Gene,
           Gene = CG_name) %>%
    select(!CG_name) 
  
  lookup <- rbind(lookup1, lookup2)
  
  
  

  # merging datasets
  
  full_list$full_data <- full_join(modeled_expression, log_expression, 
                        by = c("Timepoint", "Annotation", "Gene", "Cluster" ), 
                        suffix = c("_modeled", "_raw") ) %>% 
    mutate(value_modeled = replace_na(value_modeled, 0), # replace empty values with 0 for nicer plots
           value_sum = value_raw + value_modeled,
           value_prod = value_raw * value_modeled) %>%
    filter(Gene != "CG3809") %>% # this gene is Adk1 which is already present in the data, also barely any expression
    
    # merge in Chromosomal location (and replace some of the CG-names with meaningful names)
    left_join(lookup, by = c( "Gene" = "Gene" )) %>%
    mutate(Gene = if_else(is.na(Gene2), Gene, Gene2)) %>%
    select(!Gene2) %>%
    distinct() 
  
  
  
  return(full_list)
  
}

data_list <- prepare_data()
gc()

write_csv2(data_list$full_data, "Desplan_log_mean_expression_all.csv")
