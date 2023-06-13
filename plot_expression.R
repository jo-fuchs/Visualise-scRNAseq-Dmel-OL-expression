### Visualizing expression of single genes



library(tidyverse)

# read data
data_list <- list()
data_list$full_data <- read_csv2("Desplan_log_mean_expression_all.csv")

### define plotting functions ####

## plotting modeled expression of an individual gene  #####
plot_expression_modeled <- function(gene1) {
  
  data_list$full_data %>% 
    filter(Gene == gene1) %>%
    filter(str_detect(Annotation, pattern = "^[A-Za-z].*") ) %>% 
    # filter out clusters without any single expression
    
    
    pivot_wider(values_from = value_modeled,
                id_cols = c("Annotation"), 
                names_from = c("Timepoint", "Gene")) %>% 
    mutate(Expressed = rowMeans(select(., !c(Annotation))),
           Annotation = fct_reorder(Annotation, Expressed, .desc = F)) %>% 
     
    filter(Expressed > 0) %>% 
    select(!Expressed) %>% 
    pivot_longer(cols = !c("Annotation"), names_to = c("Timepoint", "Gene"), names_pattern = "(.*)_(.*)") %>% 
    
    
    # correct order of genes & Timepoints
    mutate(Timepoint = factor(Timepoint, levels = c("P15", "P30", "P40", "P50", "P70", "Adult"))) %>% 
    
    # plot
    ggplot(aes(y = Annotation, fill = value, x = Timepoint)) +
    geom_tile() + 
    scale_fill_viridis_c(limits= c(0,1), option = "cividis", breaks = c(0, 1), labels = c("off", "on")) + 
    theme_minimal() + 
    theme(legend.position = "right") + 
    labs(title = paste("Modeled expression values of ",  gene1, sep = ""),
         x = "",
         y = "",
         fill = "Expression")
  
}


## plotting raw (log10(expression+1)) expression of an individual gene  #####
plot_expression_raw <- function(gene1, cutoff) {
  
  plotting <- data_list$full_data %>% 
    filter(Gene == gene1) %>% 
    filter(str_detect(Annotation, pattern = "^[A-Za-z].*") ) %>%
    filter(!str_detect(Annotation, pattern = "^Im|^GMC|Apop|LQ")) %>% 
    # filter out clusters without any single expression
    pivot_wider(values_from = value_raw,
                id_cols = c("Annotation"), 
                names_from = c("Timepoint", "Gene")) %>%
    mutate(Expressed = rowMeans(select(., !c(Annotation))),
           #Annotation = fct_reorder(Annotation, Expressed, .desc = F)
           ) %>% 
    # arrange(Expressed) %>% head(20) %>%  
     filter(Expressed > cutoff) %>%
     select(!Expressed) %>% 
    pivot_longer(cols = !c("Annotation"),
                 names_to = c("Timepoint", "Gene"), 
                 names_pattern = "(.*)_(.*)") %>% 
    
    
    # correct order of genes & Timepoints
    mutate(Timepoint = factor(Timepoint, levels = c("P15", "P30", "P40", "P50", "P70", "Adult"))) %>% 
    
    # plot
    ggplot(aes(y = Annotation, fill = value, x = Timepoint)) +
    geom_tile() + 
    scale_fill_viridis_c(limits= c(0,NA))+#, na.value = "#440154FF") + 
    theme_minimal() + 
    theme(legend.position = "right") + 
    labs(title =  gene1,
         x = "",
         y = "",
         fill = "Raw\nexpression")
  return(plotting)
 # ggsave(file.path("Expression_plots", paste(gene1, "_Expression.png", sep = "")), plotting, device = "png", width = 13, height = 35, units = "cm", bg = "white" )
}



#### plot data ####


plot_expression_modeled("Inx2")
(plot_expression_raw("chp", cutoff = 0))


# to check if and how a gene is called in the dataset use this
data_list$full_data |> select(Gene) |> filter(str_detect(Gene, "Dna"))

