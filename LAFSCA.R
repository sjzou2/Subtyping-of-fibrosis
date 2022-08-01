rm(list = ls())

library(singscore)
library(factoextra)
library("caret")
library("randomForest")
library("FactoMineR")
library(plyr)
library(pheatmap)
library(pROC)
library(lsa)
library(dplyr)
library(reshape2)
library(WGCNA)
library(glmnet)
library(ggplot2)
library(ggpubr)


Singscore_mat = function(Expr, GO_geneset, cluster_frame){
  Signature_score = cluster_frame
  Signature_score$Sam_ID = rownames(cluster_frame)
  Signature_score = as.data.frame(Signature_score$Sam_ID)
  rownames(Signature_score) = Signature_score$`Signature_score$Sam_ID`
  colnames(Signature_score) = "Sam_ID"
  Geneset_name = rownames(GO_geneset)
  
  GO_geneset = t(GO_geneset)
  Expr = Expr[,rownames(cluster_frame)]
  rank_Heart_UNorm = rankGenes(Expr)
  pb <- txtProgressBar(style=3)
  for (i in 1:length(Geneset_name)) {
    Geneset = GO_geneset[,Geneset_name[i]]
    Geneset = Geneset[!duplicated(Geneset)]
    singscore = simpleScore(rank_Heart_UNorm, upSet = Geneset)
    k = as.data.frame(singscore$TotalScore)
    rownames(k) = rownames(singscore)
    colnames(k) = Geneset_name[i]
    k$Sam_ID = rownames(k)
    Signature_score = merge(Signature_score, k, by = "Sam_ID")
    setTxtProgressBar(pb, i/length(Geneset_name))
  }
  rownames(Signature_score) = Signature_score$Sam_ID
  Signature_score = Signature_score[,-which(colnames(Signature_score) == c("Sam_ID"))]
  return(Signature_score)
}


# Calculate correlation of the 10452 gene sets with EMT 
if(F){
  load("DiscoverySet_RPKM.Rdata")
  load("DiscoverySet_clustered.Rdata")
  Expr = DiscoverySet_RPKM[,rownames(DiscoverySet_clustered)]
  rank_Heart_UNorm = rankGenes(Expr)
  
  colnames(Expr) == rownames(DiscoverySet_clustered)
  
  
  Hall_mark_EMT = read.table("HALLMARK geneset EMT.txt", sep = "\t")
  Hall_mark_EMT = Hall_mark_EMT[-c(1:2),]
  
  singscore_HeartUNorm_EMT <- simpleScore(rank_Heart_UNorm, upSet = Hall_mark_EMT)
  
  
  load("HALLMARK_GO_genesets.Rdata")
  Functional_genesets = rbind.fill(GO_BP_genesets, GO_CC_genesets, GO_MF_genesets, HALLMARK_genesets)
  GN = c(rownames(GO_BP_genesets), c(rownames(GO_CC_genesets)), c(rownames(GO_MF_genesets)), c(rownames(HALLMARK_genesets)))
  rownames(Functional_genesets) = GN

  Signature_scores = Singscore_mat(DiscoverySet_RPKM, Functional_genesets, DiscoverySet_clustered)
  save(Signature_scores, file = "Discovery_Signature_scores.Rdata")
  
  
  rownames(singscore_HeartUNorm_EMT) == rownames(DiscoverySet_clustered)
  singscore_HeartUNorm_EMT$Cluster = DiscoverySet_clustered$kmeans_Heart_silhouette
  
  Signature_scores = Signature_scores[rownames(singscore_HeartUNorm_EMT),]
  rownames(Signature_scores) == rownames(singscore_HeartUNorm_EMT)
  Signature_scores$EMT_score = singscore_HeartUNorm_EMT$TotalScore
  
  Signature_scores$Cluster = singscore_HeartUNorm_EMT$Cluster
  Signature_scores_cluster1 = Signature_scores[which(Signature_scores$Cluster == "1"),]
  Signature_scores_cluster2 = Signature_scores[which(Signature_scores$Cluster == "2"),]
  
  Genesets = colnames(Signature_scores)
  Genesets = Genesets[-which(Genesets %in% c("EMT_score", "Cluster"))]
  Pvalue = c()
  Coef = c()
  
  Signature_scores_cluster1 = as.data.frame(t(na.omit(t(Signature_scores_cluster1))))
  tmp = Signature_scores_cluster1
  tmp = apply(tmp, 2, as.numeric)
  rownames(tmp) = rownames(Signature_scores_cluster1)
  Signature_scores_cluster1 = as.data.frame(tmp)
  
  Genesets = Genesets[(Genesets %in% colnames(Signature_scores_cluster1))]
  pb <- txtProgressBar(style=3)
  
  for (i in 1:length(Genesets)) {
    res = cor.test(Signature_scores_cluster1[,Genesets[i]], Signature_scores_cluster1$EMT_score, method = "spearman")
    Pvalue[i] = res$p.value
    Coef[i] = res$estimate
    setTxtProgressBar(pb, i/length(Genesets))
  }
  close(pb)
  names(Pvalue) = Genesets
  names(Coef) = Genesets
  
  All_Cluster1_cor = as.data.frame(Genesets)
  All_Cluster1_cor$Pvalue = Pvalue
  All_Cluster1_cor$Coef = Coef
  rownames(All_Cluster1_cor) = All_Cluster1_cor$Genesets
  
  Significant_gene_cluster1 = Genesets[which(Pvalue < 0.05)]
  Significant_gene_cluster1 = as.data.frame(Significant_gene_cluster1)
  Significant_gene_cluster1$Pvalue = Pvalue[Significant_gene_cluster1[,1]]
  Significant_gene_cluster1$Coeff = Coef[Significant_gene_cluster1[,1]]
  rownames(Significant_gene_cluster1) = Significant_gene_cluster1$Significant_gene_cluster1
  Significant_gene_cluster1 = Significant_gene_cluster1[order(Significant_gene_cluster1$Coeff, decreasing = T),]
  
  Positive_gene_Cluster1 = rownames(Significant_gene_cluster1)[which(Significant_gene_cluster1$Coeff >0)]
  Negative_gene_Cluster1 = rownames(Significant_gene_cluster1)[which(Significant_gene_cluster1$Coeff <0)]
  
  
  Genesets = colnames(Signature_scores)
  Genesets = Genesets[-which(Genesets %in% c("EMT_score", "Cluster"))]
  Pvalue = c()
  Coef = c()
  
  Signature_scores_cluster2 = as.data.frame(t(na.omit(t(Signature_scores_cluster2))))
  tmp = Signature_scores_cluster2
  tmp = apply(tmp, 2, as.numeric)
  rownames(tmp) = rownames(Signature_scores_cluster2)
  Signature_scores_cluster2 = as.data.frame(tmp)
  
  Genesets = Genesets[(Genesets %in% colnames(Signature_scores_cluster2))]
  pb <- txtProgressBar(style=3)
  
  for (i in 1:length(Genesets)) {
    res = cor.test(Signature_scores_cluster2[,Genesets[i]], Signature_scores_cluster2$EMT_score)
    Pvalue[i] = res$p.value
    Coef[i] = res$estimate
    setTxtProgressBar(pb, i/length(Genesets))
  }
  close(pb)
  names(Pvalue) = Genesets
  names(Coef) = Genesets
  
  All_cluster2_cor = as.data.frame(Genesets)
  All_cluster2_cor$Pvalue = Pvalue
  All_cluster2_cor$Coef = Coef
  rownames(All_cluster2_cor) = All_cluster2_cor$Genesets
  
  Significant_gene_cluster2 = Genesets[which(Pvalue < 0.05)]
  Significant_gene_cluster2 = as.data.frame(Significant_gene_cluster2)
  Significant_gene_cluster2$Pvalue = Pvalue[Significant_gene_cluster2[,1]]
  Significant_gene_cluster2$Coeff = Coef[Significant_gene_cluster2[,1]]
  rownames(Significant_gene_cluster2) = Significant_gene_cluster2$Significant_gene_cluster2
  Significant_gene_cluster2 = Significant_gene_cluster2[order(Significant_gene_cluster2$Coeff, decreasing = T),]
  
  Positive_gene_cluster2 = rownames(Significant_gene_cluster2)[which(Significant_gene_cluster2$Coeff >0)]
  Negative_gene_cluster2 = rownames(Significant_gene_cluster2)[which(Significant_gene_cluster2$Coeff <0)]
  
  
  Cluster1_specific = Positive_gene_Cluster1[(which(Positive_gene_Cluster1 %in% Positive_gene_cluster2 == F))]
  Cluster2_specific = Positive_gene_cluster2[(which(Positive_gene_cluster2 %in% Positive_gene_Cluster1 == F))]
  
  High_coeff_1 = rownames(Significant_gene_cluster1)[which(Significant_gene_cluster1$Coeff >0.5)]
  High_coeff_2 = rownames(Significant_gene_cluster2)[which(Significant_gene_cluster2$Coeff >0.5)]
  
  Cluster1_specific = Cluster1_specific[which(Cluster1_specific %in% High_coeff_1)]
  Cluster2_specific = Cluster2_specific[which(Cluster2_specific %in% High_coeff_2)]
  
  save(All_Cluster1_cor, All_cluster2_cor, Cluster1_specific, Cluster2_specific, file = "DiscoverySet_Subtype_specific_signatures.Rdata")
  
}
# Common significant genesets GSE57338 and Discovery set
if(F){
  load("DiscoverySet_Subtype_specific_signatures.Rdata")
  # Correlation of the functional gene sets with EMT scores and subtype-specific gene sets
  
  Discovery = cbind(All_Subtype1_cor, All_Subtype2_cor)
  Discovery = Discovery[union(Subtype1_specific, Subtype2_specific),]
  Discovery = Discovery[,-which(colnames(Discovery) %in% c("Genesets"))]
  colnames(Discovery) = c("Pvalue_1","Coef_1","Pvalue_2","Coef_2")
  Subtype_enrich = rownames(Discovery)
  Subtype_enrich[(Subtype_enrich %in% Subtype1_specific)] = "Subtype1"
  Subtype_enrich[(Subtype_enrich %in% Subtype2_specific)] = "Subtype2"
  Discovery$Subtype_enrich = Subtype_enrich
  
  
  load("GSE57338_geneset_Cor_val.Rdata")
  rownames(GSE57338_Correlation_Subtype1) == rownames(GSE57338_Correlation_Subtype2)
  GSE57338 = cbind(GSE57338_Correlation_Subtype1, GSE57338_Correlation_Subtype2)
  GSE57338 = GSE57338[,-3]
  colnames(GSE57338) = c("Coef_1", "Pvalue_1","Coef_2", "Pvalue_2","Subtype_enrich")
  
  
  Common_Subtype1_P = rownames(Discovery)[which(Discovery$Pvalue_1<0.05 & GSE57338$Pvalue_1<0.05)]
  Common_Subtype2_P = rownames(Discovery)[which(Discovery$Pvalue_2<0.05 & GSE57338$Pvalue_2<0.05)]

  save(Common_Subtype1_P, Common_Subtype2_P, file = "Common_Discovery_GSE57338.Rdata")
}


# CSEA of key genesets from EMT-cor-------------
if(F){
  load("DEG_Subtype1_highlow_CSEA.Rdata")
  load("DEG_Subtype2_highlow_CSEA.Rdata")
  
  load("Common_Discovery_GSE57338.Rdata")
  
  rownames(Discovery.CSEA.Subtype_1_high) = Discovery.CSEA.Subtype_1_high$names
  rownames(Discovery.CSEA.Subtype_2_high) = Discovery.CSEA.Subtype_2_high$names
  
  Subtype_1_NES1 = Discovery.CSEA.Subtype_1_high[Subtype1_specific,]
  Subtype_1_NES2 = Discovery.CSEA.Subtype_2_high[Subtype1_specific,]
  Subtype_2_NES1 = Discovery.CSEA.Subtype_1_high[Subtype2_specific,]
  Subtype_2_NES2 = Discovery.CSEA.Subtype_2_high[Subtype2_specific,]
  
  Subtype_1_NES1 = na.omit(Subtype_1_NES1)
  Subtype_2_NES2 = na.omit(Subtype_2_NES2)
  
  Subtype_1_NESScreened = rownames(Subtype_1_NES1)[which(Subtype_1_NES1$pValue<0.05)]
  Subtype_2_NESScreened = rownames(Subtype_2_NES2)[which(Subtype_2_NES2$pValue<0.05)] 
}

FurtherScreened_1 = intersect(Subtype_1_NESScreened, Common_Subtype1_P)
FurtherScreened_2 = intersect(Subtype_2_NESScreened, Common_Subtype2_P)
save(FurtherScreened_1, FurtherScreened_2, file = "Further_screened.Rdata")
