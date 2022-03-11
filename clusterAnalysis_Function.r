#Jenny Smith

#May 8, 2017 

#Purpose: Create a PCoA, MDS, nMDS ordination plots given expression data and genes of interest. 
#Updated 6/26/17 to be identical to differential expression pipeline
setwd(file.path(TARGET,"RNA/Nanostring/2017.05.04_DataExploration"))

plotPCoA <- function(expnData,clinData, geneList, factor){
  #df is the dataframe with genes as colnames, patient IDs as rownames
  #df also has the factor data, such as mutation status pos,neg. 
  #cols specifies the numeric columns with expn values. 
  #factor is the name of the factor column 
  
  library(vegan)
  df <- merge_CDE_Expn(clinData,expnData, geneList) #merge causes loss of rownames
  geneColnames <- gsub("\\-", "\\.", geneList)  #remove hyphens to match colnames
  geneColnames <- gsub("(^[0-9].+)", "X\\1", geneColnames) #add an X to genes names that begin with numbers
  
  PCoA <- capscale(df[,geneColnames] ~ 1, distance = "bray", add=TRUE)
  
  scores <- as.data.frame(scores(PCoA, display="sites"))
  
  p <- ggplot(scores, aes(x=MDS1, MDS2)) + 
    geom_point(aes(color=df[,factor])) + 
    theme_bw() + 
    labs(title="")
  
  aov <- aov(scores$MDS1 ~ df[,factor]) #NOTE: This is only valid for balanced experimental designs! Equal # of obs in each factor level. 
  
  list <- list(df, PCoA,scores, p,aov)
  names(list) <- c("MDS_df","PCoA","scores","plot","anova")
  
  return(list)
}


#Updated on 6/9/17 to use variance stabilized transformed data as input (not center scaled log2, like in princomp)
PCA <- function(expnData,phenovector,round=TRUE){
  library(DESeq2)
  #expnData has patient IDs as colnames and genes as rownames. 
  
  countData <- expnData[,match(names(phenovector), colnames(expnData))]
  countData <- round(expnData, digits = 0)
  colData <- as.data.frame(phenovector)
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ phenovector)
  
  dds <- dds[ rowSums(counts(dds)) > 5, ]
  
  varianceStab <- vst(dds, blind = TRUE)
  plot <- plotPCA(varianceStab, intgroup = "phenovector") + theme_numX
  pca.dat <- plotPCA(varianceStab, intgroup = "phenovector", returnData=TRUE)
  
  list <- list(dds, varianceStab, plot, pca.dat)
  names(list) <- c("dds", "vst", "pca_plot","pca_data")
  
  return(list)
}


# plotPCA <- function(expnData,clinData, factor,log2=TRUE){
#   #ExpnData is the expression data with patients as columns, genes as rows. 
#   #clindata  has the factor data, such as mutation status pos,neg. Patient USI as rownames
#   #cols specifies the numeric columns with expn values. 
#   #factor is the name of the factor column 
#   #log2 indicates wether the input expn data matrix is already log2 or not. 
#   
#   source("/Volumes/jlsmith3/scripts/RNAseq_Analysis/DifferentialExpn_PathwayAnalysis/merge_clinData_ExpnData_Functions.r")
#   library(genefilter)
#   library(magrittr)
#   
#   Ngenes <- length(which(colnames(expnData) %in% rownames(clinData))) - 1 #ensure that number of genes is less than #samples
#   
#   topVarGenes <- rownames(expnData[order(rowVars(expnData),decreasing=TRUE), ] %>% .[1:Ngenes, ])
#   pca_df <- merge_CDE_Expn(clinData,expnData, topVarGenes) #merge causes loss of rownames. 
#   
#   topVarGenes <- gsub("\\-", "\\.", topVarGenes) #remove hyphens to match colnames
#   topVarGenes <- gsub("(^[0-9].+)", "X\\1", topVarGenes) #add an X to genes names that begin with numbers
#   
#   
#   if (log2 == FALSE){
#     expn <- log2(pca_df[,intersect(colnames(pca_df),topVarGenes)] + 0.01)
#   }else if (log2 == TRUE){
#     expn = pca_df[,intersect(colnames(pca_df),topVarGenes)]
#   }
#   
#   pca.scaled <- scale(expn)
#   
#   pca <- princomp(pca.scaled,cor = T, scores = T)
#   scores <- as.data.frame(unclass(pca$scores))
#   percVar <- round((pca$sdev^2)/sum(pca$sdev^2)*100, digits = 2)
#   
#   pca_plot <- ggplot(scores, aes(x=scores$Comp.1, y=scores$Comp.2))
#   
#   pca_plot <- pca_plot + geom_point(aes(color=factor(pca_df[,factor]))) +
#     theme_bw() +
#     labs(title = "PCA of TARGET AML: Most Varied Genes",
#          x = paste("principal Component 1:", percVar[1], "% Variance", sep=" "),
#          y = paste("principal Component 2:  ", percVar[2], "% Variance", sep=" "))
#   
#   # scree_plot <-  plot(percVar[1:10], xlab = "Principal Component",
#   #                       ylab = "Percentage of Variance Explained",
#   #                       type = "b", pch = 16,
#   #                       main = "Percent Variance Exlpained in PCA analysis",
#   #                       col = "dark blue")
#   
#   list <- list(pca_df,pca, scores, pca_plot, topVarGenes)
#   names(list) <- c("pca_df", "pca", "scores", "pca_plot", "topVarGenes")
#   
#   return(list)
# }








