#!/usr/bin/env Rscript
#initialize packrat
#packrat::init()
#packrat::snapshot()

#install related packages if not installed yet

if(FALSE){
relatedPackages = c("stats", "Seurat", "Matrix", "psych", "gprofiler2", "optparse")
for(p in relatedPackages){
  if(!require(p,character.only = TRUE, quietly = TRUE)) 
    suppressMessages(install.packages(p))
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  suppressMessages(install.packages("BiocManager", quietly = TRUE))

if(!require("SingleCellExperiment", character.only = TRUE, quietly = TRUE))
 suppressMessages(BiocManager::install("SingleCellExperiment", quietly = TRUE))
}

suppressPackageStartupMessages(library(optparse))
#arguments list
option_list = list(
  #data file (the count matrix) as a SingleCellExperiment object
  make_option(c("-d", "--data"), type="character", default=NULL,
              help="data input", metavar="character"),
  #dataset name
  make_option(c("-n", "--dataset.name"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  #which mitochondrial threshold was used to filter the count matrix
  make_option(c("-R", "--mitoRatio"), type = "numeric", default = 0.03,
              help = "declare the mitochondrial threshold used to filter the input dataset [default= %default]", metavar = "numeric"),
  #a method to run correlation
  make_option(c("-m", "--method"), type = "character", default = "pearson",
              help = "correlation method [default= %default] Other options: kendall, spearman", metavar = "character"),
  make_option(c("-r", "--correlation_threshold"), type = "numeric", default = 0.1,
               help = "define a minimum correlation level [default = %default]", metavar = "numeric"),
  #level of significance
  make_option(c("-a", "--alpha"), type = "numeric", default = 0.05,
              help = "significance level [default= %default]", metavar = "numeric"),
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
              In case of more than one genes, provide a comma separated list of gene names, without spaces e.g. Erf,Sp7,Acaca", metavar = "character"),
  #enrichment analysis
  make_option(c("-e", "--enrichment"), type = "character", default = 'TRUE',
              help = "run enrichment analysis [default = %default]\n", metavar = "boolean"),
  #enrichment only run
  make_option(c("-o", "--only_enrichment"), type ="character", default = NULL,
              help = "an enrichment only run, if a previous enrichment analysis has failed 
              and a cor_res.rds file is available [default = %default]\n
              In this case, please provide the full path for cor_res.rds 
              e.g ~/<directory iRNA is running>/results/<timestamp>/cor_res.rds Note that you still have to provide additional information on parameters to be written on .csv file
              (-n, -R, -m etc.). If skipped, no such information will be included", metavar = "character"),
  #significant gene sets
  make_option(c("-s", "--signGeneSets"), type = "character", default = 'FALSE',
              help = "search for significant gene sets [default = %default]", metavar = "boolean")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#load packages

message('\nloading libraries, please wait')
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(psych))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(SingleCellExperiment))

message("libraries loaded\n")


#set parameters
dataset.name        = opt$dataset.name 
mitoRatio           = opt$mitoRatio 
paramCorrMethod     = opt$method 
paramCorrPAdjust    = opt$pAdjustMethod   
paramSeuMinFeatures = 200
paramGeneSparcityThr= opt$geneSparsity 
alpha.level         = opt$alpha 
enrich              = opt$enrichment
gen.set.sign        = opt$signGeneSets
enrich.only         = opt$only_enrichment
c.r                 = opt$correlation_threshold 

#an enrichment only run
if(!is.null(enrich.only)){

  message("\n Reading data...\n")
  cor.res = readRDS(enrich.only)
  output.path = gsub("(.*)[/].*", "\\1", enrich.only)

  
  library(gprofiler2)
  
  message("\n Running enrichment analysis\n")
  
  #Read cor.res by 3 to follow the routine of enrichment analysis using the selectedGene variable
  for (i in seq(1, length(cor.res), by=3)){
      
    selectedGene = names(cor.res[i])
    message("\n------Working on ", selectedGene, "------\n")

    g.sets = lapply(c(cor.res[selectedGene],
                      cor.res[paste0(selectedGene, "_pos.cor")],
                      cor.res[paste0(selectedGene, "_neg.cor")]), rownames)

    #gost_S198 = readRDS("gost_iRNA_filtered_S198.rds")
    
    gost_S198 = list()
    
    for (x in 1:length(g.sets)){
      
      if(length(g.sets[x][1][[1]])==0){
        message("\n ", names(g.sets[x]), " is empty\n")
        next
      }
      message("\n Awaiting a response from server...\n")
      #a helper function to handle errors when calling gost
      gost_ja <- function(x){
        tryCatch(
          expr = {
            gost.res <<-gost(query= g.sets[x], organism = "mmusculus", domain_scope = "annotated", significant = T, evcodes = TRUE,
                             sources = c("GO", "KEGG", "REAC", "WP", "MIRNA", "HPA", "CORUM", "HP"))
            message("Successfully executed the gost call.")
            
          },
          error = function(e){
            message('Error:')
            gost.response <<- TRUE
            
            message(e$message)
          },
          warning = function(w){
            message('Warning:')
            message(w)
          },
          finally = {
            message('\n')
          }
          
        ) 
        
      }
      gost.response <- FALSE
      gost_ja()
      b=0
      while(gost.response==TRUE && b<5) { 
        
        gost_ja()
        
        
        
        Sys.sleep(2)
        b= b+1
        if(b==5){
          message("\nEnrichment analysis failed. Writing results of correlation analysis to an .rds file,
              to use it at another enrichment only run\n")
          saveRDS(cor.res, file = file.path(the.path,"cor_res.rds"))
          stop("Exiting...",call. = FALSE)
          
        }
      }
      
      if(is.null(gost.res)){
        message("\n gene set for ", names(g.sets[x]), " didn't return results from enrichment analysis\n")
        next
      }
      
      message("\n enrichment analysis for ", names(g.sets[x])," completed successfully\n")
      gost_S198[names(g.sets[x])] <-list(as.matrix(gost.res$result[,c("p_value","term_id", "source", "term_name", "intersection")]))
      
      outputPath          = paste0(output.path,'/', selectedGene, '/')
      f.path = paste0(outputPath, paste0(names(g.sets[x]),"_", mitoRatio, "_", dataset.name, '_GO.csv'))
      #append parameters of analysis to file
      cat(paste0("paramCorrMethod: ",paramCorrMethod,"\n"), 
          paste0("correlation level: ", c.r, "\n"),
          paste0("level of significance: ", alpha.level, "\n"),
          paste0("paramCorrPAdjust: ",paramCorrPAdjust,"\n"), 
          paste0("paramGeneSparcityThr: ", paramGeneSparcityThr, "\n"), 
          paste0("dataset: ", dataset.name, "\n"), 
          #paste0("percentage of cells where the selected gene is expressed: ", sel.gene.percent, "\n"),
          paste0("mitoRatio: ", mitoRatio), file=f.path)
      #write enrichment analysis results
      #write.csv(gost_S198[genes.incl[x]], 
      #paste0(file.path(getwd(), outputPath), paste0(genes.incl[x],'_GO.csv')), row.names=FALSE)
      message("\n Writing results to file\n")
      suppressWarnings(write.table(gost_S198[names(g.sets[x])], f.path, sep=",", append=TRUE, col.names=NA))
      #wait a second before calling again gost
      Sys.sleep(1)
    } 
     
      
  } 
  stop("\r Enrichment analysis completed")
}

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
  message("\n Correlation completed")
  #keep only those genes where p values are less than alpha.level = 0.05 c.r = 0.1
  corr.sign <- corr.matrix[which(corr.matrix[,2]<alpha.level & abs(as.numeric(corr.matrix[,1]))>c.r), ,drop = FALSE ]
  if(length(corr.sign)==0)
    stop("None correlated gene returned with a < ", alpha.level, " and r > ", c.r)
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
    p.value.wtest = wilcox.test(x.set,y.set, conf.int = FALSE)$p.value

    #append the p value to the corr.sign matrix
    corr.sign[i,3] = p.value.wtest
    
  }

  
  message("\n Wilcoxon test run completed")
  
  #again keep only those genes where p values are less than alpha.level
  corr.sign.genes <- corr.sign[which(corr.sign[,3]<alpha.level), , drop = FALSE]
  if(length(corr.sign.genes)==0)
    stop("None correlated gene returned with a < ", alpha.level)
  #order the list by the r coefficient
  corr.ordered <- corr.sign.genes[order(abs(as.numeric(corr.sign.genes[,1])),decreasing = TRUE),]
  
  #keep genes with r> c.r
  corr.ordered <- corr.ordered[which(abs(as.numeric(corr.ordered[,1]))>c.r), , drop = FALSE]

  message("\n Writing results to files")
  
  #Writing results to files
  outputPath          = paste0(the.path,'/', selectedGene, '/')
 
  #check if the directory exists and if not create it
  ifelse(!dir.exists(outputPath), dir.create(outputPath), FALSE)
  
  f.path = paste0(outputPath, paste0(selectedGene,"_", mitoRatio,"_", dataset.name,'.csv'))
  
  #append parameters of analysis to file
  cat(paste0("paramCorrMethod: ",paramCorrMethod,"\n"), 
      paste0("correlation level: ", c.r, "\n"),
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

  
  f.path = paste0(outputPath, paste0(selectedGene, "_pos.cor"),"_", mitoRatio,"_", dataset.name,'.csv')
  #append parameters of analysis to file
  cat(paste0("paramCorrMethod: ",paramCorrMethod,"\n"), 
      paste0("correlation level: ", c.r, "\n"),
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
      paste0("correlation level: ", c.r, "\n"),
      paste0("level of significance: ", alpha.level, "\n"),
      paste0("paramCorrPAdjust: ",paramCorrPAdjust,"\n"), 
      paste0("paramGeneSparcityThr: ", paramGeneSparcityThr, "\n"), 
      paste0("dataset: ", dataset.name, "\n"), 
      paste0("percentage of cells where the selected gene is expressed: ", sel.gene.percent, "\n"),
      paste0("mitoRatio: ", mitoRatio), file=f.path)
  #write the ordered matrix of correlated genes with r coefficients and p values to the designated directory
  suppressWarnings(write.table( cor.res[paste0(selectedGene, "_neg.cor")], f.path, sep=",", append=TRUE, col.names=NA))
  
  #################Enrichment analysis#################
  
  if(enrich){
    
    message("\n Running enrichment analysis\n")
    g.sets = lapply(c(cor.res[selectedGene],
                      cor.res[paste0(selectedGene, "_pos.cor")],
                      cor.res[paste0(selectedGene, "_neg.cor")]), rownames)
    
    library(gprofiler2)
    #library(plotly)
    
    
    #gost_S198 = readRDS("gost_iRNA_filtered_S198.rds")
    
    gost_S198 = list()
    
    for (x in 1:length(g.sets)){
      
      if(length(g.sets[x][1][[1]])==0){
        message("\n ", names(g.sets[x]), " is empty\n")
        next
      }
      #a helper function to handle errors when calling gost
      gost_ja <- function(x){
        tryCatch(
          expr = {
            gost.res <<-gost(query= "ERF", organism = "mmusculus", domain_scope = "annotated", significant = T, evcodes = TRUE,
                             sources = c("GO", "KEGG", "REAC", "WP", "MIRNA", "HPA", "CORUM", "HP"))
            message("Successfully executed the gost call.")
            
          },
          error = function(e){
            message('Error:')
            gost.response <<- TRUE
            
            message(e$message)
          },
          warning = function(w){
            message('Warning:')
            message(w)
          },
          finally = {
            message('\n')
          }
          
        ) 
        
      }
      gost.response <- FALSE
      gost_ja()
      b=0
      while(gost.response==TRUE && b<5) { 
        
        gost_ja()
        
        
        
        Sys.sleep(2)
        b= b+1
        if(b==5){
          message("\nEnrichment analysis failed. Writing results of correlation analysis to an .rds file,
          to use it at another enrichment only run\n")
          saveRDS(cor.res, file = file.path(the.path,"cor_res.rds"))
          stop("Exiting...",call. = FALSE)
          
        }
      }
      
      if(is.null(gost.res)){
        message("\n gene set for ", names(g.sets[x]), " didn't return results from enrichment analysis\n")
        next
      }
      
      message("\n enrichment analysis for ", names(g.sets[x])," completed successfully\n")
      gost_S198[names(g.sets[x])] <-list(as.matrix(gost.res$result[,c("p_value","term_id", "source", "term_name", "intersection")]))

      f.path = paste0(outputPath, paste0(names(g.sets[x]),"_", mitoRatio, "_", dataset.name, '_GO.csv'))
      #append parameters of analysis to file
      cat(paste0("paramCorrMethod: ",paramCorrMethod,"\n"), 
          paste0("correlation level: ", c.r, "\n"),
          paste0("level of significance: ", alpha.level, "\n"),
          paste0("paramCorrPAdjust: ",paramCorrPAdjust,"\n"), 
          paste0("paramGeneSparcityThr: ", paramGeneSparcityThr, "\n"), 
          paste0("dataset: ", dataset.name, "\n"), 
          paste0("percentage of cells where the selected gene is expressed: ", sel.gene.percent, "\n"),
          paste0("mitoRatio: ", mitoRatio), file=f.path)
      #write enrichment analysis results
      #write.csv(gost_S198[genes.incl[x]], 
      #paste0(file.path(getwd(), outputPath), paste0(genes.incl[x],'_GO.csv')), row.names=FALSE)
      message("\n Writing results to file\n")
      suppressWarnings(write.table(gost_S198[names(g.sets[x])], f.path, sep=",", append=TRUE, col.names=NA))
      #wait a second before calling again gost
      Sys.sleep(1)
    }
  }
  #################Search for significant gene sets#################
  if(gen.set.sign){
    message("\n Loading gene lists to compare with\n")
    
    #####getting an error:
    #Error in header + nrows : non-numeric argument to binary operator
   # Calls: sapply -> read.csv -> read.table
    #Execution halted
    
    #read frog file and call it genes_OldProt
    genes_OldProt <- read.csv("./input_data/Frog_Clusters.csv", sep=",", header = T)
    dim_genes_OldProt <- dim(genes_OldProt)
    genes_OldProt <- apply(genes_OldProt, 2, toupper)
    genes_OldProt <- as.data.frame(genes_OldProt)
    
    #New_Targets.csv: read the 2 first rows as header and call it genes_batch
    header2 <- sapply(read.csv("./input_data/New_Targets.csv", header=F, sep=",", 
                               check.names=FALSE, nrow=2) , paste, collapse="_")
    
    genes_batch <- read.csv("./input_data/New_Targets.csv", sep=",", header = F,
                            check.names=FALSE, skip=2, col.names=header2)
    dim_genes_batch <- dim(genes_batch)
    genes_batch <- apply(genes_batch, 2, toupper)
    genes_batch <- as.data.frame(genes_batch)
    
    #ChIPseq_genes.csv: read the 2 first rows as header and call it genes_oldVSall
    header3 <- sapply(read.csv("./input_data/ChIPseq_genes.csv", header=F, sep=",", 
                               check.names=FALSE, nrow=2) , paste, collapse="_")
    
    genes_oldVSall <- read.csv("./input_data/ChIPseq_genes.csv", sep=",", header = F,
                               check.names=FALSE, skip=2, col.names=header3)
    dim_genes_OldVSall <- dim(genes_oldVSall)
    print( dim(genes_oldVSall))
    genes_oldVSall <- apply(genes_oldVSall, 2, toupper)
    genes_oldVSall <- as.data.frame(genes_oldVSall)
    
  }
  
}
end.cor.time = Sys.time()

message("\n End of analysis \n ")
end.cor.time - start.cor.time

