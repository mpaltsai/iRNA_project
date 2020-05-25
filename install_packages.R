#install required packages

relatedPackages = c("stats", "Seurat", "Matrix", "psych", "gprofiler2", "optparse")
for(p in relatedPackages){
  if(!require(p,character.only = TRUE, quietly = TRUE)) 
    suppressMessages(install.packages(p))
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  suppressMessages(install.packages("BiocManager", quietly = TRUE))

if(!require("SingleCellExperiment", character.only = TRUE, quietly = TRUE))
  suppressMessages(BiocManager::install("SingleCellExperiment", quietly = TRUE))

requiredQcAnalysis = c("ggplot2", "scales")
for(p in requiredQcAnalysis){
  if(!require(p,character.only = TRUE, quietly = TRUE)) 
    suppressMessages(install.packages(p))
}

BiocManager::install("tidyverse")
BiocManager::install("AnnotationHub")
BiocManager::install("ensembldb")
BiocManager::install("scater")
BiocManager::install("BioQC")
BiocManager::install("qvalue")
