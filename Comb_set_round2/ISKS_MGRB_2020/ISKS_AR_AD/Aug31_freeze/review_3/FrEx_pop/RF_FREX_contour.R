---
title: "rev3_FrEx"
author: "Swetansu Pattnaik"
date: "27/03/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

##libraries
```{r echo=FALSE, message=FALSE}
library(ggplot2)
library(readxl)
library(mclust)
library(gg3D)
library(dplyr)
library(kableExtra)
library(FactoMineR)

`%nin%` = Negate(`%in%`)

source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")

```

##2D plots

```{r echo=FALSE, message=FALSE}

#for case control
plot_case_cont_2D = function(df_inp){
plot12_rect_cc <- ggplot(df_inp, aes(x = PC1, y = PC2, colour = superPopulation)) + geom_point(alpha = 0.35) + theme_bw() + theme(legend.position="bottom", legend.title = element_blank()) 
plot12_rect_pop <- ggplot(df_inp, aes(x = PC1, y = PC2, colour = pred.superPop)) + geom_point() + theme_bw() + theme(legend.position="bottom", legend.title = element_blank()) 
plot12_rect_country <- ggplot(df_inp[df_inp$Country %nin% "Unknown",], aes(x = PC1, y = PC2, colour = Country)) + geom_point() + theme_bw() + theme(legend.position="bottom", legend.title = element_blank()) 

multiplot(plot12_rect_cc, plot12_rect_pop, plot12_rect_country, cols = 3)

}

#for case control
plot_case_cont_2D_PC123 = function(df_inp, alpha = NULL){
plot12_rect_cc12 <- ggplot(df_inp, aes(x = PC1, y = PC2, colour = superPopulation)) + geom_point(alpha = alpha) +
  xlim(-0.6, 0.15) + theme_bw() + theme(legend.position="bottom", legend.title = element_blank()) 
plot12_rect_cc13 <- ggplot(df_inp, aes(x = PC1, y = PC3, colour = superPopulation)) + geom_point(alpha = alpha) + xlim(-0.6, 0.15) + theme_bw() + theme(legend.position="bottom", legend.title = element_blank()) 
plot12_rect_cc23 <- ggplot(df_inp, aes(x = PC2, y = PC3, colour = superPopulation)) + geom_point(alpha = alpha) + xlim(-0.225, 0.52) + theme_bw() + theme(legend.position="bottom", legend.title = element_blank()) 

multiplot(plot12_rect_cc12, plot12_rect_cc13, plot12_rect_cc23, cols = 3)

}

#for case control QC
plot_case_cont_2D_PC123_QC = function(df_inp, alpha = NULL){
plot12_rect_cc12 <- ggplot(df_inp, aes(x = PC1, y = PC2, colour = vis_ret)) + geom_point(alpha = alpha) +
  xlim(-0.6, 0.15) + theme_bw() + theme(legend.position="bottom", legend.title = element_blank()) 
plot12_rect_cc13 <- ggplot(df_inp, aes(x = PC1, y = PC3, colour = vis_ret)) + geom_point(alpha = alpha) + xlim(-0.6, 0.15) + theme_bw() + theme(legend.position="bottom", legend.title = element_blank()) 
plot12_rect_cc23 <- ggplot(df_inp, aes(x = PC2, y = PC3, colour = vis_ret)) + geom_point(alpha = alpha) + xlim(-0.225, 0.52) + theme_bw() + theme(legend.position="bottom", legend.title = element_blank()) 

multiplot(plot12_rect_cc12, plot12_rect_cc13, plot12_rect_cc23, cols = 3)

}

#for case control (Eur population)
plot_eur_case_cont_2D = function(df_inp,alpha = NULL){
plot12_rect_cc <- ggplot(df_inp, aes(x = PC1, y = PC2, colour = superPopulation)) + geom_point(alpha = alpha) + theme_bw() + theme(legend.position="bottom", legend.title = element_blank()) 
plot12_rect_eurpop <- ggplot(df_inp, aes(x = PC1, y = PC2, colour = pred.eurPop)) + geom_point(alpha = alpha) + theme_bw() + theme(legend.position="bottom", legend.title = element_blank()) 
plot12_rect_country <- ggplot(df_inp[df_inp$Country %nin% "Unknown",], aes(x = PC1, y = PC2, colour = Country)) + geom_point() + theme_bw() + theme(legend.position="bottom", legend.title = element_blank()) 

multiplot(plot12_rect_cc, plot12_rect_eurpop, plot12_rect_country, cols = 3)

}

#for clustered output
plot_clust_2D = function(df_inp){
clust_plot12 <- ggplot(df_inp, aes(x = PC1, y = PC2, colour = clust_code)) + geom_point() +
  xlim(-0.6, 0.15) + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) 
#+ geom_vline(xintercept = -0.2, lty = 2)
clust_plot13 <- ggplot(df_inp, aes(x = PC1, y = PC3, colour = clust_code)) + geom_point() +
  xlim(-0.6, 0.15) + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) 
#+ geom_vline(xintercept = -0.2, lty = 2)
clust_plot23 <- ggplot(df_inp, aes(x = PC2, y = PC3, colour = clust_code)) + geom_point() +
  xlim(-0.225, 0.52) + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) 
#+  geom_vline(xintercept = -0.2, lty = 2)

multiplot(clust_plot12, clust_plot13, clust_plot23, cols = 3)

}

```


##Prepare FrEx PCA output

```{r}
frex_score = read.delim("~/RVAS/comb_set_2020/FrEx_pop_PCA/merge.frex.isks.pca.vs.1000g.eigenvec", 
                        sep = " ", header = F, stringsAsFactors = F)
frex_score = frex_score[,-1]
cnames = c("sample", paste("PC",1:20, sep = ""), "population", "superPopulation")
colnames(frex_score) = cnames

##sample map
map_sam = read.table("~/RVAS/comb_set_2020/FrEx_pop_PCA/samp_old_new.txt", header = F, stringsAsFactors = F, sep = "")
rect_sam_fr = frex_score$sample[1:560]
ind_1K = which(map_sam$V1 %in% frex_score$sample)
ind_coh = which(map_sam$V2 %in% frex_score$sample)
ind_coh = ind_coh[1:4665]

rect_sam = c(rect_sam_fr, map_sam[ind_1K,1], map_sam[ind_coh,1])
frex_score = frex_score[frex_score$sample %nin% "SAMP4850",]
frex_score$rect_sam = rect_sam

frex_score$population = ifelse(grepl("SAMP", frex_score$sample) & grepl("^[ABZ]", frex_score$rect_sam), "MGRB", frex_score$population)
frex_score$superPopulation = ifelse(grepl("SAMP", frex_score$sample) & grepl("^[ABZ]", frex_score$rect_sam), "MGRB", frex_score$superPopulation)

##add country
comb_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/review_add/R2Q6_popstruc/Ethnicity.xlsx",
                         sheet = 1, col_types = c("list"))
comb_pheno <- as.data.frame(comb_pheno)
comb_pheno1 <- sapply(comb_pheno, unlist)
colnames(comb_pheno1) <- colnames(comb_pheno)
comb_pheno <- comb_pheno1
comb_pheno <- as.data.frame(comb_pheno, stringsAsFactors = F)
comb_pheno <- unique(comb_pheno)
comb_pheno <- comb_pheno[!is.na(comb_pheno$pid),]
comb_pheno$`age at dateExtracted` <- as.numeric(comb_pheno$`age at dateExtracted`)
comb_pheno$AgeatSarcoma <- as.numeric(comb_pheno$AgeatSarcoma)
comb_pheno$SubjectAgeCancer <- as.numeric(comb_pheno$SubjectAgeCancer)

frex_score_fin = frex_score[frex_score$rect_sam %nin% map_sam[ind_1K,1],] ##remove 1000 genome
frex_score_fin$Country <- comb_pheno[match(frex_score_fin$rect_sam, comb_pheno$pmn),64]
frex_score_fin$Country <- ifelse(frex_score_fin$superPopulation %nin% c("MGRB", "ISKS"), "France",
                                 ifelse(frex_score_fin$superPopulation %in% c("MGRB"), "Aus", frex_score_fin$Country))

ggplot(frex_score_fin, aes(x = PC1, y = PC2, colour = superPopulation)) + geom_point(alpha = 1) + theme_bw() + theme(legend.position="bottom", legend.title = element_blank())

ggplot(frex_score_fin, aes(x = PC1, y = PC2, colour = Country)) + geom_point(alpha = 1) + theme_bw() + theme(legend.position="bottom", legend.title = element_blank())

frex_score_fin_chk = frex_score_fin[frex_score_fin$Country %in% "France",]

ggplot(frex_score_fin_chk, aes(x = PC1, y = PC2, colour = superPopulation)) + geom_point(alpha = 1) + theme_bw() + theme(legend.position="bottom", legend.title = element_blank())


```

##check if the samples have hits in discovered complexes
```{r}
##Variant file used for SKAT
scores = read.table("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_combset_pca.scores_clustered_rect_sam.tsv", header = T, sep = "\t", stringsAsFactors = F)

fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv", 
                      sep = "\t", header = T, stringsAsFactors = F)

fil_tab <- fil_tab[fil_tab$SAMPLE %in% scores$rect_sam,]
#fil_tab <- fil_tab[fil_tab$SAMPLE %in% scores_PCnew$sample,]


##chk_gene function

chk_gene = function(samp_inp){
  genes_outlier <- unique(fil_tab[fil_tab$SAMPLE %in% samp_inp,]$gene_symbol)
  Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")
#  Shelterin <- c("POT1", "TINF2", "TERF1", "SMARCAL1", "STAG3")
CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
MPNST_pos <- c("NF1", "LZTR1", "SDHA", "SDHB", "SDHD")
sarcoma_genes <- c("TP53", "SDHA", "SDHB", "SDHD")

genes_outlier_df <- unique(fil_tab[fil_tab$SAMPLE %in% samp_inp & 
                                     fil_tab$gene_symbol %in% c(Shelterin, CEP_HAUS_core, MPNST_pos, "TP53"),c(1,9,127)])

return(genes_outlier_df)
}

```





##With RF classifier output

```{r}

##RF classifier with 6 PCs
##Random Forest classification (multiclass)
#setwd("~/RVAS/comb_set_2020/pop_PCA/review_2/") ##copied code from ~/RVAS/comb_set_2020/pop_PCA/review_2/pca_pop_cluster.R

library(randomForest)
library(caret)
library(e1071)

scores = frex_score
scores$superPopulation = ifelse(scores$superPopulation %in% c("BORDEAUX", "BREST", "DIJON","LILLE", "NANTES", "ROUEN"), "FrEx", scores$superPopulation)
#scores_subPC1 = scores[scores$PC1 >= -0.2,]
#scores_ISKS_MGRB = scores_subPC1[scores_subPC1$superPopulation %in% c("ISKS","MGRB", "RISC", "LIONS"),]

##global population: Use first 4 PCs
train_pop <- scores[scores$superPopulation %nin% c("ISKS","MGRB", "FrEx"),c(23,2:7)]
train_pop$superPopulation = as.factor(train_pop$superPopulation)

test_pop <- scores[scores$superPopulation %in% c("ISKS","MGRB", "FrEx"),c(23,2:7)]
test_pop$superPopulation = as.factor(test_pop$superPopulation)
set.seed(94162240)
rf <- randomForest(superPopulation~., data=train_pop, importance = TRUE, proximity = TRUE, ntree = 10000) 
print(rf)
varImp(rf)
##prediction
##evaluate model

##superPopulation
prediction <-predict(rf, test_pop[,-1])
prob_pop = predict(rf, test_pop[,2:7], type="prob")
colnames(prob_pop) = paste("pred.superPop.", colnames(prob_pop), sep = "")
test_pop$pred.superPop = prediction
scores_ISKS_MGRB_FrEx = scores[scores$superPopulation %in% c("ISKS","MGRB", "FrEx"), 1:24]
scores_ISKS_MGRB_FrEx_pop = cbind.data.frame(scores_ISKS_MGRB_FrEx, "pred.superPop" = prediction, prob_pop)

##EUR population
##European population

train_pop_eur <- scores[scores$superPopulation %in% "EUR",c(22,2:7)]
train_pop_eur$population = as.factor(train_pop_eur$population)

test_pop_eur <- scores[scores$superPopulation %in% c("ISKS","MGRB", "FrEx"),c(22,2:7)]
test_pop_eur$population = as.factor(test_pop_eur$population)

set.seed(94162240)
rf_eur <- randomForest(population~., data=train_pop_eur, importance = TRUE, proximity = TRUE, ntree = 10000) 
print(rf_eur)
varImp(rf_eur)

##prediction
##evaluate model

prediction_eur <- predict(rf_eur, test_pop_eur[,-1])
prob_pop_eur = predict(rf_eur, test_pop_eur[,2:7], type="prob")
colnames(prob_pop_eur) = paste("pred.eurPop.", colnames(prob_pop_eur), sep = "")
test_pop_eur$pred.eurPop = prediction_eur
scores_ISKS_MGRB_pop_eur = cbind.data.frame(scores_ISKS_MGRB_pop, "pred.eurPop" = prediction_eur, prob_pop_eur)

scores_ISKS_MGRB_pop_eur$pred.NFE = FALSE
scores_ISKS_MGRB_pop_eur$pred.NFE[scores_ISKS_MGRB_pop_eur$pred.superPop == "EUR" & scores_ISKS_MGRB_pop_eur$pred.eurPop != "FIN"] = TRUE
scores_ISKS_MGRB_pop_eur$prob.NFE = 1 - scores_ISKS_MGRB_pop_eur$pred.eurPop.FIN


#ggplot(scores_ISKS_MGRB_pop, aes(x = PC1, y = PC2, colour = pred.superPop)) + geom_point() + theme_bw() + theme(legend.position="bottom", legend.title = element_blank())
#make_3Dplot_pop(scores_ISKS_MGRB_pop, -120, 0)

scores_ISKS_MGRB_pop_eur$superPopulation = ifelse(scores_ISKS_MGRB_pop_eur$superPopulation %in% "MGRB", "controls", "sarcoma")
##Add country
scores_ISKS_MGRB_pop_eur$Country <- comb_pheno[match(scores_ISKS_MGRB_pop_eur$rect_sam, comb_pheno$pmn),64]
scores_ISKS_MGRB_pop_eur$Country[is.na(scores_ISKS_MGRB_pop_eur$Country)] = "Aus_MGRB"

plot_case_cont_2D(scores_ISKS_MGRB_pop_eur)
case_cont_ratio_pop(scores_ISKS_MGRB_pop_eur)
#make_3Dplot_pop(scores_ISKS_MGRB_pop_eur, -120, 0)
#make_3Dplot_eurpop(scores_ISKS_MGRB_pop_eur, -120, 0)

write.table(scores_ISKS_MGRB_pop_eur, "~/RVAS/comb_set_2020/pop_PCA/review_2/MGRB_ISKS_combset_pca.scores_RF_clustered.tsv", col.names = TRUE, row.names = FALSE, quote = F, sep = "\t")


##EUR only based on Epi study
scores_ISKS_MGRB_Eur = scores_ISKS_MGRB_pop_eur[scores_ISKS_MGRB_pop_eur$pred.superPop %in% "EUR",]
plot_eur_case_cont_2D(scores_ISKS_MGRB_Eur,alpha = 1)

```

