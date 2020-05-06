#!/usr/bin/env Rscript
#initialize packrat
#packrat::init()

#install related packages
#install.packages(c("stats", "Seurat", "Matrix", "psych", "gprofiler2", "plotly"))
#install.packages("optparse")
#take a snapshot 
#packrat::snapshot()

#clean working environment
#rm(list=ls())
suppressPackageStartupMessages(library(optparse))
#arguments list
option_list = list(
  make_option(c("-d", "--dataset"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-R", "--mitoRatio"), type = "numeric", default = 0.03,
              help = "declare the mitochondrial threshold used to filter the input dataset [default= %default]", metavar = "numeric"),
  make_option(c("-m", "--method"), type = "character", default = "pearson",
              help = "correlation method [default= %default]\n 
              other options: kendall, spearman", metavar = "character"),
  make_option(c("-p", "--pAdjustMethod"), type = "character", default = "fdr",
              help = "p value adjustment [default = %default]", metavar = "character"),
  make_option(c("-g", "--geneSparsity"), type = "numeric", default = 0.02,
              help = "define the percentage of cells a gene is expressed [default = %default]", metavar = "numeric"),
  make_option(c("-t", "--targetGenes"), type = "character", default = 'Erf',
              help = "define target gene(s) [default = %default]\n
              In case of more than one genes, provide a comma separated list of gene names\n
              e.g. Erf, Sp7, Acaca", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

message('loading libraries, please wait')
#load packages
#suppressPackageStartupMessages(library(shiny))
#suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(stats))
#suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
#suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(psych))
suppressPackageStartupMessages(library(gprofiler2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(optparse))
message("libraries loaded")






start.time = Sys.time()