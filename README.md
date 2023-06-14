# Visualising scRNAseq gene expression data of developmental Dropophila optic lobe 

This repo contains functions to visualise scRNAseq expression data of developmental Drosophila melanogaster optic lobe from Desplan-lab (corresponding publication <https://www.nature.com/articles/s41586-020-2879-3>).

Original data accessed at: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142787>. The authors published two datasets:

-   log normalized raw expression data per cluster (*GSE142787_Log_normalized_average_expression.xlsx*)
-   binarized expression (*GSE142787_Mixture_modeling.xlsx*)\
    according to <https://elifesciences.org/articles/50901> binarized expression is more reliable descriptors of protein abundance in a cluster

The metadata (annotation of the clusters) is described in: *41586_2020_2879_MOESM4_ESM.xlsx*

## Using the plotting functions

These raw data can be merged and annotated using *0_prepare_data.R* from this repository using tidyverse functions. For plotting just use the functions in *plot_expression.R*:

*plot_expression_raw()* will produce a heatmap of expression values over time for a gene of interest in all annotated clusters of the dataset.

*plot_expression_modeled()* will produce a heatmap of binarized expression values for a gene of interest for only those clusters, where at least one timepoint is predicted to express the gene.

The plots below show raw and modeled expression values for chaoptin (chp), a known photoreceptor-specific gene.

+------------------------------------------------------------------------------+------------------------------------------------------------------+
| ![Raw expression values of an example gene](chp_Expression.png){width="400"} | ![modeled expression values of an example gene](chp_modeled.png) |
+------------------------------------------------------------------------------+------------------------------------------------------------------+
