rm(list = ls())

library(stringr)
library(singscore)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(GEOquery)
library(Rtsne)
library(ggfortify)
library(monocle)
library(factoextra)
library(org.Hs.eg.db)
library(EDASeq)
library(reshape2)
library(WGCNA)

# setting working disk
#setwd("")

# load CIBERSORTx calculated matrix
load("CIBERSOTx_DiscoverySet.Rdata")

# Cluster by CIBERSORTX
if(F){
  fviz_nbclust(CIBERSOTx_DiscoverySet, kmeans, method = "wss") + geom_vline(xintercept = 5, linetype = 2)
  fviz_nbclust(CIBERSOTx_DiscoverySet, kmeans, method = "silhouette")
  km_Heart_wss <- kmeans(CIBERSOTx_DiscoverySet, 5, nstart = 50)
  km_Heart_silhouette <- kmeans(CIBERSOTx_DiscoverySet, 2, nstart = 50)
  
  save(km_Heart_silhouette, file = "DiscoverySet_KMeans.Rdata")
  CIBERSOTx_DiscoverySet$kmeans_Heart_wss = km_Heart_wss$cluster
  CIBERSOTx_DiscoverySet$kmeans_Heart_silhouette = km_Heart_silhouette$cluster
  Cluster = CIBERSOTx_DiscoverySet$kmeans_Heart_silhouette
  Cluster1 = Cluster
  Cluster1[which(Cluster == 1)] = 2
  Cluster1[which(Cluster == 2)] = 1
  CIBERSOTx_DiscoverySet$kmeans_Heart_silhouette = Cluster1
  set.seed(2)
  Heart_tsne_out <- Rtsne(CIBERSOTx_DiscoverySet[,1:(ncol(CIBERSOTx_DiscoverySet)-2)],pca=FALSE,perplexity=30,theta=0.0)
  str(Heart_tsne_out)
  Heart_tsnes=Heart_tsne_out$Y
  colnames(Heart_tsnes) <- c("tSNE1", "tSNE2")
  Heart_tsnes=as.data.frame(Heart_tsnes)
  Heart_tsnes$cluster_silu = as.character(CIBERSOTx_DiscoverySet$kmeans_Heart_silhouette)
  Heart_tsnes$cluster_wss = as.character(CIBERSOTx_DiscoverySet$kmeans_Heart_wss)
  ggplot(Heart_tsnes, aes(x = tSNE1, y = tSNE2, colour = cluster_silu))+ geom_point() + theme_bw()
  ggplot(Heart_tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col = cluster_wss))
  
  Heart_clustered = CIBERSOTx_DiscoverySet
  Heart_clustered$kmeans_Heart_wss = as.character(Heart_clustered$kmeans_Heart_wss)
  Heart_clustered$kmeans_Heart_silhouette = as.character(Heart_clustered$kmeans_Heart_silhouette)
  
  Heart_clustered$'Fibrosis_score' = Heart_UNorm_combine$`Fibrosis score`
  save(Heart_clustered,file = "Heart_clustered.Rdata")
}

