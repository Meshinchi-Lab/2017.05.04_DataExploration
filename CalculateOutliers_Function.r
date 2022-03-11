#Jenny Smith

#June 13, 2017 

#purpose: Given a expression matrix with group information column, detect outliers (75% plus 1.5*IQR) or greater than normal BM 
setwd(file.path(TARGET,"RNA/Nanostring/2017.05.04_DataExploration"))

detHighExpressors <- function(expnMatrix, phenoVector, geneList){
  library(magrittr)
  #order and subset if necessary
  
  expnMatrix <- expnMatrix[, intersect(names(phenoVector), colnames(expnMatrix))]
  
  if (ncol(expnMatrix)  < length(phenoVector)){
    phenoVector <- intersect(colnames(expnMatrix), names(phenoVector))
  }
  
  tmp <- data.frame(t(expnMatrix), 
                    Group=phenoVector)
  
  aboveNorm <- function(tmp, gene){
    # print(gene)
    Norm <- tapply(tmp[,gene], INDEX = list(tmp$Group), max) 
    Norm <- Norm[grepl("BM", names(Norm))]
    
    outlierExpn <- tmp[which(tmp[,gene] > Norm & tmp$Group == "EOI"),]
    outlierIDs <- rownames(outlierExpn) %>% gsub("\\.1", "", .)
    
    
    SimilarToBM <- tmp[which(tmp$Group == "DX"), ]
    SimilarToBM <- setdiff(rownames(SimilarToBM), outlierIDs)
    
    list <- list(outlierIDs, SimilarToBM)
    names(list) <- c("AboveMaxBM", "similarToBM")
    
    return(list)
  }
  
  groups <- lapply(geneList, aboveNorm, tmp=tmp)
  names(groups) <- geneList
  
  phenoVectors_MultipleGroups <- function(listofgoupsIDs){
    library(magrittr)
    #listofgoupsIDs contains a list with each item containing the IDs for each factor level.
    #See GroupIDs function - this produced the input called  "listofgroupIDs"
    group <- names(listofgoupsIDs)
    
    vector <- NULL
    names <- NULL
    for (i in 1:length(listofgoupsIDs)){
      g <- group[i]
      vector <- c(vector, rep(g, length(listofgoupsIDs[[i]])))
      names <- c(names, listofgoupsIDs[[i]])
    }
    
    names(vector) <- names
    return(vector)
  }
  
  phenoVectors <- lapply(groups, phenoVectors_MultipleGroups)
  anno_df <- do.call(cbind, phenoVectors)
  
  return(anno_df)
  
}


