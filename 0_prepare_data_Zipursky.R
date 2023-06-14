#### 

# using scRNA-seq data from Zipursky lab
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156455
# 
# Two datasets: 
#     early (P0 - P25)
#     late (P25-P100)
# 
# corresponding publication
# https://pubmed.ncbi.nlm.nih.gov/33125872/ 
#
#


# preparation of raw data of each individual cell 
# averaging per cluster and log10-transformed (log10(mean(rawcount+1)))
# this only needs to be run once to create the averaged dataset

library(tidyverse)


  # create a list to hold all data
  full_list <- list()
  
  ## data preparation #### 
  # preprocessing cluster names
  metadata_early <- read_tsv("GSE156455_tsne_early.tsv") 
  metadata_main <- read_tsv("GSE156455_tsne_main.tsv") 
  
  features <- read_tsv("GSE156455_features_main.tsv")
  
  # load expression matrices
  expression_early <- Matrix::readMM("GSE156455_matrix_early.mtx")
  expression_main <- Matrix::readMM("GSE156455_matrix_main.mtx")
  
  
  # merge data
  dimnames(expression_early) <- list(Gene = features$Gene,
                                    metadata_early$barcode)
  dimnames(expression_main) <- list(Gene = features$Gene,
                                     metadata_main$barcode)
  
  expression_early_table <- as_tibble(as.matrix(expression_early), rownames = "Gene")
  expression_main_table <- as_tibble(as.matrix(expression_main), rownames = "Gene")
  
  # this is a tiny litte bit memory expensive..
  mean_log_expression_early <- expression_early_table %>%
    select(!Gene) %>% 
    pivot_longer(cols = everything(), names_to = "barcode") %>% 
    left_join(metadata_early) %>% 
    group_by(Gene, type, time, class) %>% 
    summarise(value = log10(mean(value + 1))) # to make it strictly positive
  
  
  write_csv(mean_log_expression_early, "Log_mean_expression_early.csv")
  
  
  
  # main expression is too large to process in one batch -> split by genes
  pb = txtProgressBar(min = 0, max = 12999, initial = 1, style = 3) 
  for (i in 1:12999) { 
    expression_main_table[i,]  %>% 
      pivot_longer(cols = !Gene, names_to = "barcode") %>% 
      left_join(metadata_main, by = "barcode") %>% 
      group_by(Gene, type, time, class) %>% 
      summarise(value = log10(mean(value + 1)), .groups = "keep") %>% # to make it strictly positive
      write_csv("Log_mean_expression_main.csv", append = T)    
    setTxtProgressBar(pb,i)
    Sys.sleep(time = 1)
    close(pb)   
    print(i)
  }
 
  # batch2
  pb = txtProgressBar(min = 13000, max = nrow(expression_main_table), initial = 13000, style = 3) 
  for (i in 13000:nrow(expression_main_table)) { 
    expression_main_table[i,]  %>% 
      pivot_longer(cols = !Gene, names_to = "barcode") %>% 
      left_join(metadata_main, by = "barcode") %>% 
      group_by(Gene, type, time, class) %>% 
      summarise(value = log10(mean(value + 1)), .groups = "keep") %>% # to make it strictly positive
      write_csv("Log_mean_expression_main2.csv", append = T)    
    setTxtProgressBar(pb,i)
    Sys.sleep(time = 1)
    close(pb)   
    print(i)
  }
  
  
  # and now combine the three files and filter and save again
  early <- read_csv("Log_mean_expression_early.csv", skip = 1, col_names = c("Gene", "Cluster", "Timepoint", "Class", "Expression"))
  part1 <- read_csv("Log_mean_expression_main.csv", col_names = c("Gene", "Cluster", "Timepoint", "Class", "Expression"))
  part2 <- read_csv("Log_mean_expression_main2.csv", col_names = c("Gene", "Cluster", "Timepoint", "Class", "Expression"))
  
  expression_all_long <- bind_rows(part1, part2, early)
  
  
  expression_all_long %>% distinct() %>% write_csv("Log_mean_expression_all.csv")
  