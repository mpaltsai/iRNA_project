#!/usr/bin/env Rscript
#initialize packrat
#packrat::init()

#install related packages
#install.packages(c("stats", "Seurat", "Matrix", "psych", "gprofiler2", "plotly"))
#install.packages("optparse")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("SingleCellExperiment")
#take a snapshot 
#packrat::snapshot()

#clean working environment
#rm(list=ls())
suppressPackageStartupMessages(library(optparse))
#arguments list
option_list = list(
  #data file (the count matrix) as a SingleCellExperiment object
  make_option(c("-d", "--data"), type="character", default=NULL,
              help="data input", metavar="character"),
  #dataset name
  make_option(c("-n", "--dataset.name"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  #output file
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
  #which mitochondrial threshold was used to filter the count matrix
  make_option(c("-R", "--mitoRatio"), type = "numeric", default = 0.03,
              help = "declare the mitochondrial threshold used to filter the input dataset [default= %default]", metavar = "numeric"),
  #a method to run correlation
  make_option(c("-m", "--method"), type = "character", default = "pearson",
              help = "correlation method [default= %default]\n 
              other options: kendall, spearman", metavar = "character"),
  #level of significance
  make_option(c("-a", "--alpha"), type = "numeric", default = 0.05,
              help = "significance level [default= %default]", metavar = "character"),
  #a correction method for p values returned from correlation
  make_option(c("-p", "--pAdjustMethod"), type = "character", default = "fdr",
              help = "p value adjustment [default = %default]", metavar = "character"),
  #Include features detected in at least this many cells. It will be used as input to the argument min.cells of
  #the function CreateSeuratObject (Seurat package)
  make_option(c("-g", "--geneSparsity"), type = "numeric", default = 0.02,
              help = "define the percentage of cells a gene is expressed [default = %default]", metavar = "numeric"),
  #gene(s) to run the analysis over
  make_option(c("-t", "--targetGenes"), type = "character", default = 'Erf',
              help = "define target gene(s) [default = %default]\n
              In case of more than one genes, provide a comma separated list of gene names, without spaces\n
              e.g. Erf,Sp7,Acaca", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

message('\nloading libraries, please wait')
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
suppressPackageStartupMessages(library(SingleCellExperiment))
#suppressPackageStartupMessages(library(optparse))
message("libraries loaded\n")


#set parameters
dataset.name        = opt$dataset.name
mitoRatio           = opt$mitoRatio
paramCorrMethod     = opt$method
paramCorrPAdjust    = opt$pAdjustMethod   
paramSeuMinFeatures = 200
paramGeneSparcityThr= opt$geneSparsity
alpha.level         = opt$alpha


#read data
sce <- readRDS(opt$data)
paramSeuMinCells    = paramGeneSparcityThr*dim(sce)[2]

#creaye a seurat object to work with
seurat                 = CreateSeuratObject(counts = sce@assays@data@listData$counts, project = "craniosynostosis", min.cells = paramSeuMinCells)

#discard genes that are not detected in any cell
#Ismini note: I think is already implemented with paramSeuMinCells
seurat.sub.nonzero = seurat[rowSums(seurat)!=0]

#discard cells that none gene is detected ### NO NEED, it has been taken care of when the seurat object was created, with min. features argument
#or with thw subsetting when checking the assay for QC
seurat.sub.nonzero = seurat.sub.nonzero[,which(colSums(seurat.sub.nonzero)!=0)]

###insert normalized data in the seurat object
#after revoving zero-rows and columns
seurat.sub.nonzero <-  NormalizeData(seurat.sub.nonzero)

#initialize a list to keep correlation analysis results
cor.res <- list()

#capitalize genes names in the seurat object
#rownames(x=seurat.sub.nonzero$RNA) <- toupper(rownames(x=seurat.sub.nonzero$RNA))
seurat.sub.nonzero$RNA@counts@Dimnames[1][[1]] <- toupper(dimnames(seurat.sub.nonzero)[1][[1]])
seurat.sub.nonzero$RNA@data@Dimnames[1][[1]] <- toupper(dimnames(seurat.sub.nonzero)[1][[1]])

#set target genes
genes.incl = strsplit(opt$targetGenes, ",")
genes.incl = genes.incl[[1]]
genes.incl = toupper(genes.incl)

message("\nAnalysis will run with ", length(genes.incl), " gene(s):\n")

cat(genes.incl, sep = " ", "\n")


#initialize a list to keep genes excluded from analysis
genes.xcl = c()
#setdiff(genes.incl, intersect(names(cor.res), genes.incl))

#create a directory to store results of the form ~/directory-running-iRNA/results/timestamp
#check if the directory exists and if not create it
res.dir = ifelse(!dir.exists(file.path(getwd(), "results")), dir.create(file.path(getwd(), "results")), FALSE)
#create a sub directory, name it with a timestamp
st=format(Sys.time(), "%Y-%m-%d_%H:%M")
dir.create(file.path(getwd(), "results", st))
the.path = file.path(getwd(), "results", st)
message("\nCreating a directory to store results\n")
message(" ",the.path )

start.cor.time = Sys.time()

for( x in 1:length(genes.incl)){
  
  selectedGene         = genes.incl[x]
  message("\n ------Working on ", selectedGene, "------")

  refRowIndex = which(row.names(seurat.sub.nonzero)==selectedGene)
  
  #in case a target gene is not included in the assay (or after the subsetting of the matrix)
  if(length(refRowIndex)==0){ 
  message("\n ERROR: the selected gene ", selectedGene, " is not included in the RNA-seq assay\n Proceeding to the next gene")
  next
  }
  
  ######subsetting data##############
  
  #the row of the selected gene is in table div
  div = seurat.sub.nonzero[refRowIndex,]
  #div = div[["RNA"]]@data DON"T DO THAT!!!!! it skips rownames and colnames. Stick with the seurat object!!!!
  
  #Ea table holds only cells where the selected gene is NOT expressed
  
  #we must check if the selected gene is expressed in ALL cells, otherwise Ea gives an error
  
  #this variable is needed to go to the next iteration if becomes TRUE
  skip_to_next <- FALSE
  
  check.it <- tryCatch( div[,which(div[["RNA"]]@data[1,]==0)], error = function(e) { 
    
    e
    genes.xcl <<- c(genes.xcl, selectedGene)
    skip_to_next <<- TRUE 
    
    e$message <- message(e, selectedGene, " is expressed in ALL cells")
  })
  
  if(skip_to_next) { next }
  
  #if everything is fine continue with Ea
  Ea = div[,which(div[["RNA"]]@data[1,]==0)]
  
  #Eb holds the cells where the selected gene is expressed
  Eb = div[, which(div[["RNA"]]@data[1,]>0)]
  
  #table F: a subset of the original data  where the selected gene is not zero
  F = seurat.sub.nonzero[,colnames(Eb)]
  #exclude the row of the selected gene
  F = F[-refRowIndex,]
  
  #compute the percentage of cells (in the whole dataset) where the selected gene is expressed
  sel.gene.percent = ncol(Eb)/ncol(div)
  #Check for correlation between the selected gene and all other genes in all cells in F
  #use the transpose matrices
  Eb.t = as.matrix(t(Eb[["RNA"]]@data))
  
  F.t = as.matrix(t(F[["RNA"]]@data))
  
  #make a function to supress warnings
  cor.fun <- function(){
    #initialize a matrix to keep the results of the correlation test and the wilcoxon test
    corr.matrix = matrix(nrow = nrow(F), ncol = 5)
    rownames(corr.matrix) <- rownames(F)
    colnames(corr.matrix) <- c("r", "p_r", "p_w", "percentage in Ea", "percentage in Eb")
    
    for(j in 1:nrow(F)){
      
      #run the correlation test
      corr.output =  corr.test(as.numeric(Eb.t), as.numeric(F.t[,j]), method = paramCorrMethod, 
                               adjust = paramCorrPAdjust, alpha = alpha.level, ci = FALSE )
      #append results to corr.matrix
      #corr.matrix[j,1] = colnames(F.t)[j]
      corr.matrix[j,1] = corr.output$r
      corr.matrix[j,2] = corr.output$p
      
    }
    corr.matrix
  }
  
  message("\n Running correlation...")
  corr.matrix = suppressWarnings(cor.fun())
  message("\n Correlation ended")
  #keep only those genes where p values are less than 0.05
  corr.sign <- corr.matrix[which(corr.matrix[,2]<0.05),]
  
  ####check whether the distibution of expression of each correlated gene in cells 
  ####where the selected gene is expressed or not is the same 
  
  ###run a Wilcoxon test###
  message("\n Running a Wilcoxon test...")
  for(i in 1:nrow(corr.sign)){
    
    #define the subsets of the data to include expression values for the correlated genes 
    
    #in cells where the selected gene is not expressed
    x.set = as.numeric(seurat.sub.nonzero[["RNA"]]@data[rownames(corr.sign)[i], colnames(Ea)])
    
    #the percentage of cells in Ea the selected gene is expressed
    percent.Ea = length(x.set[which(x.set>0)])/length(x.set)
    #append it to corr.sign
    corr.sign[i,4] = percent.Ea
    
    #where the selected gene is expressed
    y.set = as.numeric(seurat.sub.nonzero[["RNA"]]@data[rownames(corr.sign)[i], colnames(Eb)])
    
    #the percentage of cells in Eb the selected gene is expressed
    percent.Eb= length(y.set[which(y.set>0)])/length(y.set)
    #append it to corr.sign
    corr.sign[i,5] = percent.Eb
    
    #run the test
    p.value.wtest = wilcox.test(x.set,y.set)$p.value
    
    #append the p value to the corr.sign matrix
    corr.sign[i,3] = p.value.wtest
    
  }
  message("\n Wilcoxon test run completed")
  #again keep only those genes where p values are less than 0.05
  corr.sign.genes <- corr.sign[which(corr.sign[,3]<0.05),]
  
  #order the list by the r coefficient
  corr.ordered <- corr.sign.genes[order(abs(as.numeric(corr.sign.genes[,1])),decreasing = TRUE),]
  
  #keep genes with r>0.1
  corr.ordered <- corr.ordered[which(abs(as.numeric(corr.ordered[,1]))>0.1),]

  message("\n Writing results to files")
  
  #Writing results to files
  outputPath          = paste0(the.path,'/', selectedGene, '/')
 
  #check if the directory exists and if not create it
  ifelse(!dir.exists(outputPath), dir.create(outputPath), FALSE)
  
  f.path = paste0(outputPath, paste0(selectedGene,"_", mitoRatio,"_", dataset.name,'.csv'))
  
  #append parameters of analysis to file
  cat(paste0("paramCorrMethod: ",paramCorrMethod,"\n"), 
      paste0("level of significance: ", alpha.level, "\n"),
      paste0("paramCorrPAdjust: ",paramCorrPAdjust,"\n"), 
      paste0("paramGeneSparcityThr: ", paramGeneSparcityThr, "\n"), 
      paste0("dataset: ", dataset.name, "\n"), 
      paste0("percentage of cells where the selected gene is expressed: ", sel.gene.percent, "\n"),
      paste0("mitoRatio: ", mitoRatio), file=f.path)
  
  #write the ordered matrix of correlated genes with r coefficients and p values to the designated directory
 suppressWarnings(write.table(corr.ordered, f.path, sep=",", append=TRUE, col.names=NA))
  #write.csv(corr.ordered, f.path)
  
  #append to cor.res
  cor.res[selectedGene] <- list(corr.ordered)
  cor.res[paste0(selectedGene, "_pos.cor")] <-list(corr.ordered[which(corr.ordered[,"r"]>0), , drop =FALSE])
  cor.res[paste0(selectedGene, "_neg.cor")] <-list(corr.ordered[which(corr.ordered[,"r"]<0), , drop =FALSE])
  #cor.res[selectedGene]
  
  f.path = paste0(outputPath, paste0(selectedGene, "_pos.cor"),"_", mitoRatio,"_", dataset.name,'.csv')
  #append parameters of analysis to file
  cat(paste0("paramCorrMethod: ",paramCorrMethod,"\n"), 
      paste0("level of significance: ", alpha.level, "\n"),
      paste0("paramCorrPAdjust: ",paramCorrPAdjust,"\n"), 
      paste0("paramGeneSparcityThr: ", paramGeneSparcityThr, "\n"), 
      paste0("dataset: ", dataset.name, "\n"), 
      paste0("percentage of cells where the selected gene is expressed: ", sel.gene.percent, "\n"),
      paste0("mitoRatio: ", mitoRatio), file=f.path)
  #write the ordered matrix of correlated genes with r coefficients and p values to the designated directory
  suppressWarnings(write.table( cor.res[paste0(selectedGene, "_pos.cor")], f.path, sep=",", append=TRUE, col.names=NA))
  
  f.path = paste0(outputPath, paste0(selectedGene, "_neg.cor"),"_", mitoRatio,"_", dataset.name,'.csv')
  #append parameters of analysis to file
  cat(paste0("paramCorrMethod: ",paramCorrMethod,"\n"), 
      paste0("level of significance: ", alpha.level, "\n"),
      paste0("paramCorrPAdjust: ",paramCorrPAdjust,"\n"), 
      paste0("paramGeneSparcityThr: ", paramGeneSparcityThr, "\n"), 
      paste0("dataset: ", dataset.name, "\n"), 
      paste0("percentage of cells where the selected gene is expressed: ", sel.gene.percent, "\n"),
      paste0("mitoRatio: ", mitoRatio), file=f.path)
  #write the ordered matrix of correlated genes with r coefficients and p values to the designated directory
  suppressWarnings(write.table( cor.res[paste0(selectedGene, "_neg.cor")], f.path, sep=",", append=TRUE, col.names=NA))
  
  
}
end.cor.time = Sys.time()

message("\n End of correlation analysis \n ")
end.cor.time - start.cor.time

