rm(list = ls())

library(pheatmap)
library(ggplot2)
library(org.Hs.eg.db)
library(edgeR)
library(singscore)
library(Rtsne)
library(factoextra)
library(reshape2)
library(ggpubr)
library(pROC)
library(caret)
library(e1071)
library(glmnet)
library(cluster)
library("randomForest")
library("FactoMineR")


Calculate_Subtype = function(CIBERSORT_MAT, centroids){
  Subtype = c()
  Subtype_1_dist = c()
  Subtype_2_dist = c()
  Subtype_1_centroid = centroids$Subtype1
  names(Subtype_1_centroid) = rownames(centroids)
  Subtype_2_centroid = centroids$Subtype2
  names(Subtype_2_centroid) = rownames(centroids)
  CIBERSORT_MAT = t(CIBERSORT_MAT)
  for (i in 1:ncol(CIBERSORT_MAT)) {
    a = CIBERSORT_MAT[,i]
    dist_mat_1 = rbind(a, Subtype_1_centroid)
    Dist_subtype1 = dist(dist_mat_1)
    Subtype_1_dist[i] = Dist_subtype1
    dist_mat_2 = rbind(a, Subtype_2_centroid)
    Dist_subtype2 = dist(dist_mat_2)
    Subtype_2_dist[i] = Dist_subtype2
    if(Dist_subtype1 < Dist_subtype2){
      Subtype[i] = "1"
    }
    else if(Dist_subtype1 > Dist_subtype2){
      Subtype[i] = "2"
    }
  }
  names(Subtype) = colnames(CIBERSORT_MAT)
  Subtype = as.data.frame(Subtype)
  Subtype$Subtype_1_dist = Subtype_1_dist
  Subtype$Subtype_2_dist = Subtype_2_dist
  
  return(Subtype)
}


# Check centroids
if(F){
  setwd("C:/Users/sjzou2/OneDrive - City University of Hong Kong/Macrophage and fibrosis Bioinformatics/Redo")
  load("Heart_clustered_update.Rdata")
  load("DiscoverySet_KMeans.Rdata")
  Discovery_CIBERSORT = Heart_clustered_update[,-which(colnames(Heart_clustered_update) %in% c("kmeans_Heart_silhouette","Fibrosis_score","High_low","Further_divide"))]
  Discovery_CIBERSORT = t(Discovery_CIBERSORT)
  rownames(Discovery_CIBERSORT) == rownames(centroids)
  Discovery_CIBERSORT = cbind(Discovery_CIBERSORT, centroids)
  Discovery_CIBERSORT = as.data.frame(t(Discovery_CIBERSORT))
  set.seed(2)
  Heart_tsne_out <- Rtsne(Discovery_CIBERSORT,pca=FALSE,perplexity=30,theta=0.0)
  
  str(Heart_tsne_out)
  Heart_tsnes=Heart_tsne_out$Y
  colnames(Heart_tsnes) <- c("tSNE1", "tSNE2")
  Heart_tsnes=as.data.frame(Heart_tsnes)
  Type = rownames(Discovery_CIBERSORT)
  for (i in 1:length(Type)) {
    if(Type[i] %in% rownames(Heart_clustered_update)){
      Type[i] = Heart_clustered_update[Type[i],"kmeans_Heart_silhouette"]
    }
    else{Type[i] = Type[i]}
  }
  
  Subtype_1_centroid = as.data.frame(Heart_tsnes[which(Heart_tsnes$cluster == "Subtype1"),])
  Subtype_2_centroid = as.data.frame(Heart_tsnes[which(Heart_tsnes$cluster == "Subtype2"),])
  
  Heart_tsnes$cluster = Type
  ggplot(Heart_tsnes, aes(x = tSNE1, y = tSNE2, colour = cluster))+ geom_point(data = Subtype_1_centroid, aes(x = tSNE1, y = tSNE2), color = "blue", size = 5) + geom_point() + 
    geom_point(data = Subtype_2_centroid, aes(x = tSNE1, y = tSNE2), color = "Purple", size = 5)+ theme_bw()
}


# Identify cluster for validating set by distance to centroids
if(F){
  setwd("C:/Users/sjzou2/OneDrive - City University of Hong Kong/Macrophage and fibrosis Bioinformatics/Redo")
  load("DiscoverySet_KMeans.Rdata")
  
  centroids = km_Heart_silhouette$centers
  
  rownames(centroids) = c("Subtype2", "Subtype1")
  centroids = as.data.frame(t(centroids))
  
  
  # Clustering GSE57338
  if(F){
    setwd("C:/Users/sjzou2/OneDrive - City University of Hong Kong/Macrophage and fibrosis Bioinformatics/Redo/Cross dataset validation")
    GSE57338_CIBERSORT = read.csv("GSE57338 CIBERSORT.csv", header = TRUE)
    rownames(GSE57338_CIBERSORT) = GSE57338_CIBERSORT$Mixture
    GSE57338_CIBERSORT = GSE57338_CIBERSORT[,-which(colnames(GSE57338_CIBERSORT) %in% c("Mixture", "P.value", "Correlation", "RMSE"))]
    
    Subtype_GSE57338 = Calculate_Subtype(GSE57338_CIBERSORT, centroids)
    
    save(Subtype_GSE57338, file = "Subtype_GSE57338.Rdata")
  }
  

}


