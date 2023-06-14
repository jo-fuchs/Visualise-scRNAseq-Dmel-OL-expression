
library(tidyverse)


# load log-transformed expression data per cluster
mean_log_expression_all <- read_csv("Log_mean_expression_all.csv") |> 
  mutate(Cluster = if_else(Cluster == "N136", "LC14?", Cluster)) # cluster N136 likely corresponds to LC14 cluster in Desplan dataset


# find the name of the gene of interest
mean_log_expression_all|> select(Gene) |>  distinct() |>  filter(str_detect(Gene, "Act"))



plot_expression <- function(gene) {
  # plotting the expression for a certain gene
  
  plotting<-  mean_log_expression_all |> 
    filter(Gene == gene) |> 
    filter(!str_detect(Cluster, pattern = "^N.*|^x.*|^C..+|R184")) |>  
    mutate(Cluster = fct_rev(Cluster)) |> 
    #filter(str_detect(Cluster, pattern = "R7.8|R1.6|L1|L2|L3|L4|L5")) %>% 
    ggplot(aes(x = Timepoint, y = Cluster, fill = Expression)) +
    geom_tile() +
    #facet_wrap(~time) +
    theme_minimal() + 
    #ggprism::theme_prism()+
    scale_fill_viridis_c(limits = c(0,NA), 
                         na.value = "yellow") +
    labs(y = "",
         title = gene)
  return(plotting) 
  # ggsave(paste(gene, "_expression_Zipursky.png", sep = ""), plotting, device = "png", width = 13, height = 20, units = "cm", bg = "white" )
}

(plot_expression("chp"))

