rm(list = ls())

library(pheatmap)
library(ggplot2)
library(org.Hs.eg.db)
library(edgeR)
library(singscore)
library(factoextra)
library(reshape2)
library(ggpubr)
library(pROC)
library(caret)
library(e1071)
library(glmnet)
library(uniConSig)
library(clusterProfiler)
library(plyr)


load("DiscoverySet_RPKM.Rdata")
load("DiscoverySet_clustered.Rdata")
load("TTEST_Subtype_different_secreted_proteins.Rdata")


Expr = DiscoverySet_RPKM[,rownames(DiscoverySet_clustered)]
rank_Heart_UNorm = rankGenes(Expr)


Combine = union(Subtype_1_high, Subtype_2_high)
# TTEST Combine high LASSO
if(F){
  Combine = Combine[which(Combine %in% rownames(Expr))]
  Combine_expr = t(Expr[Combine,])
  Cluster = DiscoverySet_clustered$kmeans_Heart_silhouette
  Cluster[which(Cluster == "1")] = 0
  Cluster[which(Cluster == "2")] = 1
  
  fit_cyt = glmnet(Combine_expr, Cluster, family = "binomial")
  plot(fit_cyt, xvar = "lambda", label = TRUE)
  plot(fit_cyt, xvar = "dev", label = TRUE)
  
  cv.lasso <- cv.glmnet(Combine_expr, Cluster,type.measure = "class", nfolds = nrow(Combine_expr), alpha = 1, family = "binomial")
  plot(cv.lasso)
  
  
  print(cv.lasso)
  coef(cv.lasso,s=cv.lasso$lambda.min)
  coef(cv.lasso,s=cv.lasso$lambda.1se)
  Selected_Combinehigh_min = rownames(coef(cv.lasso,s=cv.lasso$lambda.min))[which(coef(cv.lasso,s=cv.lasso$lambda.min) != 0)]
  Selected_Combinehigh_min = Selected_Combinehigh_min[-which(Selected_Combinehigh_min == "(Intercept)")]
  
  Selected_Combinehigh_1se = rownames(coef(cv.lasso,s=cv.lasso$lambda.1se))[which(coef(cv.lasso,s=cv.lasso$lambda.1se) != 0)]
  Selected_Combinehigh_1se = Selected_Combinehigh_1se[-which(Selected_Combinehigh_1se == "(Intercept)")]
}
Up_combine = Selected_Combinehigh_1se[which(Selected_Combinehigh_1se %in% Subtype_1_high)]
Down_combine = Selected_Combinehigh_1se[which(Selected_Combinehigh_1se %in% Subtype_2_high)]

save(Up_combine, Down_combine, file = "Subtype_classification_Biomarkers.Rdata")

# AUC Discovery
if(F){
  rank_Heart_UNorm = rankGenes(Expr)
  singscore_HeartUNorm_cytokine_cluster <- simpleScore(rank_Heart_UNorm, upSet = Up_combine, downSet = Down_combine)
  
  rownames(singscore_HeartUNorm_cytokine_cluster) == rownames(DiscoverySet_clustered)
  Cluster = DiscoverySet_clustered$kmeans_Heart_silhouette
  Cluster[which(Cluster == "1")] = 0
  Cluster[which(Cluster == "2")] = 1
  Cluster = as.numeric(Cluster)
  singscore_HeartUNorm_cytokine_cluster$Cluster = Cluster
  
  model1 <- glm(singscore_HeartUNorm_cytokine_cluster$Cluster~singscore_HeartUNorm_cytokine_cluster$TotalScore, data=singscore_HeartUNorm_cytokine_cluster, family='binomial')
  pre <- predict(model1,type='response')
  modelroc <- roc(singscore_HeartUNorm_cytokine_cluster$Cluster,pre)
  plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
       grid.col=c("green", "red"), max.auc.polygon=TRUE,
       auc.polygon.col="springgreen", print.thres=TRUE, main="TTEST secret Discovery")
  # c("CLEC11A","MADCAM1","GGACT") = 0.748
  
  
  p <- ggboxplot(singscore_HeartUNorm_cytokine_cluster, x = "Cluster", y = "TotalScore",
                 color = "Cluster", palette = "jama",
                 add = "jitter")
  p+stat_compare_means(method = "t.test")
}
# AUC GSE57338
if(F){
  load("Subtype_GSE57338.Rdata")
  
  rankGSE57338 = rankGenes(GSE57338_expr)
  Score_28_GSE57338 = simpleScore(rankGSE57338, upSet = Up_combine, downSet = Down_combine)
  rownames(Score_28_GSE57338) == rownames(Subtype_GSE57338)
  Subtype = Subtype_GSE57338$Subtype
  Subtype[(Subtype == "1")] = 0
  Subtype[(Subtype == "2")] = 1
  
  Score_28_GSE57338$Subtype = as.numeric(Subtype)
  
  model1 <- glm(Score_28_GSE57338$Subtype~Score_28_GSE57338$TotalScore, data=Score_28_GSE57338, family='binomial')
  pre <- predict(model1,type='response')
  modelroc <- roc(Score_28_GSE57338$Subtype,pre)
  plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
       grid.col=c("green", "red"), max.auc.polygon=TRUE,
       auc.polygon.col="springgreen", print.thres=TRUE, main="TTEST secret Discovery")
  
  p <- ggboxplot(Score_28_GSE57338, x = "Subtype", y = "TotalScore",
                 color = "Subtype", palette = "jama",
                 add = "jitter")
  p+stat_compare_means(method = "t.test")
}


