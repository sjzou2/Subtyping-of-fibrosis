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
