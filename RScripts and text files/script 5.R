install.packages("gplots")
#install.packages("GSEABase")
BiocManager::install('GSEABase')
#install.packages("Biobase")
BiocManager::install('Biobase')
BiocManager::install("GSVA")
BiocManager::install("gprofiler2")
#install.packages("gprofiler2")
BiocManager::install("clusterProfiler")
#install.packages("clusterProfiler")
#install.packages("msigdbr")
BiocManager::install("msigdbr")
BiocManager::install("enrichplot")
#install.packages("enrichplot")
BiocManager::install("qusage", force=TRUE)
install.packages("patchwork")
install.packages("query")

library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
msigdbr
library(enrichplot) # great for making the standard GSEA enrichment plots
library(qusage) # Quantitative Set Analysis for Gene Expression
library(heatmaply)

myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=50, sort.by="logFC")
nrow(myTopHits)
myTopHits
#also used website
gost.res_up <- gost(rownames(myTopHits), organism = "slycopersicum", correction_method = "fdr")

gostplot(gost.res_up, interactive = T, capped = T)
gost.res_down <- gost(rownames(myTopHits), organism = "slycopersicum", correction_method = "fdr")
gostplot(gost.res_down, interactive = T, capped = T)

publish_gosttable(
   gost.res,
   highlight_terms = NULL,
   use_colors = TRUE,
   show_columns = c("source", "term_name", "term_size", "intersection_size"),
   filename = NULL,
   ggplot=TRUE)

sl_gsea_c2 <- msigdbr(species = "slycopersicum", # change depending on species your data came from
                      category = "C2") %>% # choose your msigdb collection of interest
dplyr::select(gs_name, gene_symbol)
#don't have one of the available species, so doesn't work


