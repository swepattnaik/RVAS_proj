#!/usr/bin/env Rscript
#.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

library(ggplot2)
library(readxl)
`%nin%` = Negate(`%in%`)
#setwd("~/RVAS/comb_set_2020/pop_PCA/")
#comb_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/review_add/SupplementaryPIDfile129Nov2020.xlsx",
#                      sheet = 1, col_types = c("list"))
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

##QC pass filter
QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv", 
                      header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %in% "TRUE",]
QC2_dat_pass$isFemale <- ifelse(QC2_dat_pass$f_stat < 0.2, 1, 
                                ifelse(QC2_dat_pass$f_stat > 0.8, 0, 2))


p_Data <- read.table("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_combset_pca.scores_clustered.tsv", 
                     header = T, sep = "\t", stringsAsFactors = F)
#p_Data_1000G <- p_Data[p_Data$superPopulation %nin% c("ISKS", "RISC", "LIONS", "MGRB"),]
#write.table(p_Data_1000G$sample, "~/RVAS/comb_set_2020/pop_PCA/samp_1000G.txt",
#                                 row.names = F, quote = F, sep = "")

p_Data <- p_Data[p_Data$superPopulation %in% c("ISKS", "RISC", "LIONS", "MGRB"),]
p_Data_MGRB <-  p_Data[p_Data$superPopulation %in% c("MGRB"),]

#remove duplicates
dup_samp <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/fin_samp_dup_drop.txt", sep = "", header = T,
                       stringsAsFactors = F)
dup_samp <- dup_samp$x[grepl("^[ABZ]", dup_samp$x)]
p_Data_MGRB <- p_Data_MGRB[p_Data_MGRB$sample %nin% dup_samp,]
p_Data_MGRB$rect_sam <- p_Data_MGRB$sample
p_Data_MGRB <- p_Data_MGRB[p_Data_MGRB$rect_sam %in% QC2_dat_pass$new_sampleid,]
##rename p_Data
p_Data_ISKS <-  p_Data[p_Data$superPopulation %in% c("ISKS", "RISC", "LIONS"),]
rect_sam_dat <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/ISKS_RISC_LIONS_final_freeze.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
p_Data_ISKS$rect_sam <- rect_sam_dat[match(p_Data_ISKS$sample,rect_sam_dat$JCInputID),9]
p_Data_ISKS <- p_Data_ISKS[!is.na(p_Data_ISKS$rect_sam),]
# write.table(p_Data_ISKS, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/p_Data_ISKS.tsv",
#             sep = "\t", quote = F, row.names = F)
##two samples missed due to mislabelling and QC2
#rect_sam_dat$JCInputRecID[rect_sam_dat$JCInputRecID %nin% p_Data_ISKS$sample]
p_Data_ISKS <- p_Data_ISKS[p_Data_ISKS$rect_sam %in% QC2_dat_pass$new_sampleid,] ##3105 is lost(QC2 fail)
p_Data_noCH <- rbind.data.frame(p_Data_ISKS, p_Data_MGRB)
# write.table(p_Data_noCH, "~/RVAS/comb_set_2020/pop_PCA/rect_sam_MGRB_ISKS_no1000G_combset_pca.scores_clustered.tsv", sep = "\t", row.names = F,
#             quote = F)
scores <- p_Data_noCH

#scores$Country <- comb_pheno[match(scores$sample, comb_pheno$pmn),120]
scores$Country <- comb_pheno[match(scores$sample, comb_pheno$pmn),64]
scores$Country <- ifelse(is.na(scores$Country), "Aus_MGRB", scores$Country)
scores$Country <- ifelse((!grepl("^[ABZ]",scores$sample) & (scores$Country %in% "Aus_MGRB")), "Unknown", scores$Country)

ggplot(scores, aes(x = PC1, y = PC2, colour = Country)) + geom_point() + theme_bw()

scores = scores[scores$Country %nin% "Unknown",]

scores$superPopulation <- ifelse(scores$superPopulation %in% "MGRB", "MGRB", "ISKS")
scores$superPopulation <- ifelse(scores$superPopulation %in% "MGRB", "controls", "sarcoma")

# g12 <- ggplot(scores, aes(x = PC1, y = PC2, colour = Country)) + geom_point() + theme_bw() + 
#   theme(legend.position="bottom", legend.title = element_blank()) + 
#   geom_vline(xintercept = -0.2, lty = 2) + geom_hline(yintercept = -0.1, lty = 2)
g12 <- ggplot(scores, aes(x = PC1, y = PC2, colour = Country)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = -0.2, lty = 2) 
g13 <- ggplot(scores, aes(x = PC1, y = PC3, colour = Country)) + geom_point() + theme_bw() +
  geom_vline(xintercept = -0.2, lty = 2) + 
  theme(legend.position="bottom", legend.title = element_blank()) 
g23 <- ggplot(scores, aes(x = PC2, y = PC3, colour = Country)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) 
# g23 <- ggplot(scores, aes(x = PC2, y = PC3, colour = Country)) + geom_point() + theme_bw() + 
#   theme(legend.position="bottom", legend.title = element_blank()) + 
#   geom_vline(xintercept = -0.1, lty = 2) + geom_hline(yintercept = -0.05, lty = 2)

source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")
#multiplot(g12, g13, g23, cols = 1)
##stat ellipse
##factor: predicted superPopulation

ggplot(scores, aes(x = PC1, y = PC2, colour = pred.superPop)) + 
  geom_point() + theme_bw() + stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(level = 0.95) +
  coord_fixed()

p12 <- ggplot(scores, aes(x = PC1, y = PC2, colour = pred.superPop)) + geom_point() + theme_bw() + 
  stat_ellipse(type = "norm", linetype = 2) + theme(legend.position="bottom") +
  theme(legend.title = element_blank()) + 
  stat_ellipse(level = 0.95) + geom_vline(xintercept = -0.2, lty = 2) 
p13 <- ggplot(scores, aes(x = PC1, y = PC3, colour = pred.superPop)) + geom_point() + theme_bw() + 
  stat_ellipse(type = "norm", linetype = 2) + theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  stat_ellipse(level = 0.95) + geom_vline(xintercept = -0.2, lty = 2) 
p23 <- ggplot(scores, aes(x = PC2, y = PC3, colour = pred.superPop)) + geom_point() + theme_bw() + 
  stat_ellipse(type = "norm", linetype = 2) + theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  stat_ellipse(level = 0.95) 

#multiplot(p12, p13, p23, cols = 3)

##factor: MGRB or ISKS
#scores$superPopulation <- ifelse(scores$superPopulation %in% "MGRB", "MGRB", "ISKS")
#scores$superPopulation <- ifelse(scores$superPopulation %in% "MGRB", "controls", "sarcoma")
# k12 <- ggplot(scores, aes(x = PC1, y = PC2, colour = superPopulation)) + geom_point() + theme_bw() + 
#   theme(legend.position="bottom", legend.title = element_blank()) + 
#   geom_vline(xintercept = -0.2, lty = 2) + geom_hline(yintercept = -0.1, lty = 2)
k12 <- ggplot(scores, aes(x = PC1, y = PC2, colour = superPopulation)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = -0.2, lty = 2)
k13 <- ggplot(scores, aes(x = PC1, y = PC3, colour = superPopulation)) + geom_point() + theme_bw() +
  geom_vline(xintercept = -0.2, lty = 2) + 
  theme(legend.position="bottom", legend.title = element_blank()) 
k23 <- ggplot(scores, aes(x = PC2, y = PC3, colour = superPopulation)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) 
# k23 <- ggplot(scores, aes(x = PC2, y = PC3, colour = superPopulation)) + geom_point() + theme_bw() + 
#   theme(legend.position="bottom", legend.title = element_blank()) + 
#   geom_vline(xintercept = -0.1, lty = 2) + geom_hline(yintercept = -0.05, lty = 2)

#multiplot(k12, k13, k23, cols = 3)

##combine plots
#multiplot(p12, k12, p13,  k13, p23, k23, cols = 3)

#multiplot(p12, k12, g12, p13, k13, g13, p23, k23, g23, cols = 3) 
multiplot(g12, k12, p12, g13, k13, p13, g23, k23, p23, cols = 3)


#source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")
#multiplot(g12, g13, g23, cols = 1)

##Subset based on PC1

scores_1000sub_PC1 <- scores[scores$PC1 < -0.2,]
scores_1000sub_PC23 <- scores[scores$PC2 < -0.07 & scores$PC3 < 0.05,] #use this
table(scores_sub_PC2$rect_sam %in% scores_sub_PC1$rect_sam)
amr_samp = scores[scores$pred.superPop %in% "AMR",]$rect_sam
scores_sub_PC1_AMR = scores[scores$rect_sam %in% amr_samp & scores$rect_sam %nin% scores_1000sub_PC1$rect_sam,]

write.table(scores[scores$PC1 < -0.2,]$rect_sam, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/review_add/R2Q6_popstruc/Results/rm_13_outlier_list.txt", sep = "", row.names = F, quote = F)

##Variant file used for SKAT
fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv", 
                      sep = "\t", header = T, stringsAsFactors = F)

fil_tab_PC1_sub <- fil_tab[fil_tab$SAMPLE %nin% scores_1000sub_PC1$rect_sam,]

write.table(fil_tab_PC1_sub, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze_rev_PC1_sub.tsv",
            sep = "\t", quote = F, row.names = F)


genes_outlier <- unique(fil_tab[fil_tab$SAMPLE %in% scores_1000sub_PC1$rect_sam,]$gene_symbol)

Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")
CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
MPNST_pos <- c("NF1", "LZTR1", "SDHA", "SDHB", "SDHD")
table(genes_outlier %in% Shelterin)
genes_outlier[genes_outlier %in% Shelterin]  #"TIMELESS", "TERF2"; both C3
genes_outlier[genes_outlier %in% CEP_HAUS_core] ##None
genes_outlier[genes_outlier %in% MPNST_pos] #"NF1"; C4

fil_tab[fil_tab$SAMPLE %in% scores_sub_PC1$rect_sam & fil_tab$gene_symbol %in% "NF1",c(1:9,11,127:130)] #1833
fil_tab[fil_tab$SAMPLE %in% scores_sub_PC1$rect_sam & fil_tab$gene_symbol %in% c("TIMELESS", "TERF2"),c(1:9,11,127:130)]
#1836, 2091

##For AMR sample
table(fil_tab[fil_tab$SAMPLE %in% amr_samp & fil_tab$gene_symbol %in% Shelterin,]$auto_call)
table(fil_tab[fil_tab$SAMPLE %in% amr_samp & fil_tab$gene_symbol %in% CEP_HAUS_core,]$auto_call)

##All SKAT genes
ppi_res_fil_final <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/ppi_res_fil_final_SKATbin_Aug31.tsv", sep = "\t", header = T)
table(ppi_res_fil_final$gene %in% genes_outlier)
skat_outlier <- as.character(ppi_res_fil_final$gene[ppi_res_fil_final$gene %in% genes_outlier])
fil_tab_skat_outlier <- fil_tab[fil_tab$SAMPLE %in% scores_sub_PC1$rect_sam & fil_tab$gene_symbol %in% skat_outlier, ]


##cases with Shelterin

fil_tab_Shel <- fil_tab[fil_tab$gene_symbol %in% Shelterin & fil_tab$set %in% "ISKS_AR_AD",]$SAMPLE

table(scores[scores$rect_sam %in% fil_tab_Shel, ]$Country)
nonEUR_cases <- scores[scores$rect_sam %in% fil_tab_Shel & scores$pred.NFE %in% "FALSE", ]$rect_sam
fil_tab[fil_tab$SAMPLE %in% nonEUR_cases & fil_tab$gene_symbol %in% Shelterin,]$gene_symbol

##cases with CEP HAUS
fil_tab_CEP <- fil_tab[fil_tab$gene_symbol %in% CEP_HAUS_core & fil_tab$set %in% "ISKS_AR_AD",]$SAMPLE
table(scores[scores$rect_sam %in% fil_tab_CEP, ]$Country)
table(scores[scores$rect_sam %in% fil_tab_CEP, ]$pred.superPop)
nonEUR_cases <- scores[scores$rect_sam %in% fil_tab_CEP & scores$pred.NFE %in% "FALSE", ]$rect_sam
fil_tab[fil_tab$SAMPLE %in% nonEUR_cases & fil_tab$gene_symbol %in% CEP_HAUS_core,]$gene_symbol


##Predicted superfamilies in ISKS versus MGRB
table(scores[scores$Country %in% "Aus_MGRB",]$pred.superPop) #MGRB
table(scores[scores$Country %nin% "Aus_MGRB",]$pred.superPop) #ISKS


###Proportion of variance explained
##With 1000G
eigenval_1000G <- read.delim("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_pca_combset.evals.txt",
                             sep = ":", header = F, stringsAsFactors = F)
eigenval_1000G$V1 <- gsub("^.*PC", "", eigenval_1000G$V1)
eigenval_1000G$V1 <- as.numeric(gsub("'", "", eigenval_1000G$V1))
eigenval_1000G$V2 <- gsub("}", "", eigenval_1000G$V2)
eigenval_1000G$V2 <- as.numeric(gsub(",", "", eigenval_1000G$V2))
eigenval_1000G <- eigenval_1000G[order(eigenval_1000G$V1, decreasing = F),]
eigenval_1000G$prop_var <- eigenval_1000G$V2/sum(eigenval_1000G$V2)

colnames(eigenval_1000G)[1:2] <- c("PC", "Eigen_val") 

ggplot(eigenval_1000G, aes(PC, prop_var)) +
  geom_bar(stat = "Identity") +
  xlab("Principal Components") +
  ylab("Variance Explained")






for (dest in c("svg", "png"))
{
    message(dest)
    if (dest == "svg")
        svg("combset_isks_mgrb_1000G_scores_%02d.svg", height = 8, width = 8)
    else if (dest == "png")
        png("combset_isks_mgrb_1000G_scores_%02d.png", height = 800, width = 800)
    print(ggplot(scores, aes(x = PC1, y = PC2, colour = superPopulation)) + geom_point() + theme_bw())
    print(ggplot(scores, aes(x = PC1, y = PC3, colour = superPopulation)) + geom_point() + theme_bw())
    print(ggplot(scores, aes(x = PC1, y = PC4, colour = superPopulation)) + geom_point() + theme_bw())
    print(ggplot(scores, aes(x = PC2, y = PC3, colour = superPopulation)) + geom_point() + theme_bw())
    print(ggplot(scores, aes(x = PC2, y = PC4, colour = superPopulation)) + geom_point() + theme_bw())
    print(ggplot(scores, aes(x = PC3, y = PC4, colour = superPopulation)) + geom_point() + theme_bw())
    dev.off()
}




################################
##Hierarchical clustering with bootstrapping to check the health of each cluster

scores_sub_PC1 = scores[scores$PC1 >= -0.2,]

scores_sub_dist = scores_sub_PC1[,19:21]
rownames(scores_sub_dist) = scores_sub_PC1$rect_sam
scores_sub_PC1_hc = hclust(dist(scores_sub_PC1[,19:21],method = "euclidean", diag = T), method = "centroid")
#plot(scores_sub_PC1_hc)

memb <- cutree(scores_sub_PC1_hc, k = 5) ## as there are 5 superpopulations

##bootstrap check
library(fpc)   
# set the desired number of clusters                               
kbest.p<-5       

#   Run clusterboot() with hclust 
#   ('clustermethod=hclustCBI') using Ward's method 
#   ('method="ward"') and kbest.p clusters 
#   ('k=kbest.p'). Return the results in an object 
#   called cboot.hclust.
cboot.hclust <- clusterboot(dist(scores_sub_PC1[,19:21], method = "euclidean", diag = T),clustermethod=disthclustCBI,
                            method="centroid", k=kbest.p, B = 1000)
saveRDS(cboot.hclust, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_2/cboot.hclust.rds",
        compress = T)
#   The results of the clustering are in 
#   cboot.hclust$result. The output of the hclust() 
#   function is in cboot.hclust$result$result. 
#
#   cboot.hclust$result$partition returns a 
#   vector of clusterlabels. 
#groups <- cboot.hclust$result$partition  

# Values close to 1 indicate stable clusters
cboot.hclust$bootmean 

##p value for clusters
# library(pvclust)
# set.seed(455123)
# 
# res.pv <- pvclust(scores_sub_PC1[,19:21], method.dist="euclidean", 
#                   method.hclust="centroid", nboot = 100)
# pvpick(res.pv)
table(memb)
library(tidyverse)
scores_sub_PC1_mem = scores_sub_PC1 %>%
  mutate(cluster = memb)

case_cont_df = list()
for(i in 1:length(unique(scores_sub_PC1_mem$cluster))){
  #print(i)
  clust_prop = as.data.frame(table(scores_sub_PC1_mem[scores_sub_PC1_mem$cluster == i,]$pred.superPop))
  clust_coh_prop = as.data.frame(table(scores_sub_PC1_mem[scores_sub_PC1_mem$cluster == i,]$superPopulation))
  cont = sum(clust_coh_prop[clust_coh_prop$Var1 %in% "MGRB", ]$Freq)
  case = sum(clust_coh_prop[clust_coh_prop$Var1 %nin% "MGRB", ]$Freq)
  
  #print(case)
  #print(cont)
  ratio = (case/1631)/(cont/3205)
  #print(ratio)
  
  case_cont_df[[i]] = cbind.data.frame("case" = case, "control" =  cont, "ratio" = ratio, "cluster" = i)
}

case_cont_df_fin = do.call("rbind.data.frame", case_cont_df)

case_cont_df_fin$bootmean = cboot.hclust$bootmean
View(case_cont_df_fin)
write.table(case_cont_df_fin, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/review_2/Results/case_cont_df_fin_clust.tsv",
            sep = "\t", row.names = F, quote = F)
scores_sub_PC1_mem$clust_code = paste("C", scores_sub_PC1_mem$cluster, sep = "_")
clust_plot12 <- ggplot(scores_sub_PC1_mem, aes(x = PC1, y = PC2, colour = clust_code)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = -0.2, lty = 2)
clust_plot13 <- ggplot(scores_sub_PC1_mem, aes(x = PC1, y = PC3, colour = clust_code)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = -0.2, lty = 2)
clust_plot23 <- ggplot(scores_sub_PC1_mem, aes(x = PC2, y = PC3, colour = clust_code)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = -0.2, lty = 2)

multiplot(clust_plot12, clust_plot13, clust_plot23, cols = 3)

scores_sub_PC1_mem$superPopulation <- ifelse(scores_sub_PC1_mem$superPopulation %in% "MGRB", "MGRB", "ISKS")
scores_sub_PC1_mem$superPopulation <- ifelse(scores_sub_PC1_mem$superPopulation %in% "MGRB", "controls", "sarcoma")
write.table(scores_sub_PC1_mem, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/review_2/Results/scores_sub_PC1_mem_clust.tsv",
            sep = "\t", row.names = F, quote = F)
scores_sub_PC1_mem = read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/review_2/Results/scores_sub_PC1_mem_clust.tsv",
           sep = "\t", header = T, stringsAsFactors = F)
##For AMR sample from cluster2 (K= 5)
amr_samp_c2 = scores_sub_PC1_mem[scores_sub_PC1_mem$cluster == 2,]$rect_sam
amr_samp_c5 = scores_sub_PC1_mem[scores_sub_PC1_mem$cluster == 5,]$rect_sam
table(fil_tab[fil_tab$SAMPLE %in% amr_samp_c5 & fil_tab$gene_symbol %in% Shelterin,]$auto_call)
table(fil_tab[fil_tab$SAMPLE %in% amr_samp_c5 & fil_tab$gene_symbol %in% CEP_HAUS_core,]$auto_call)
table(fil_tab[fil_tab$SAMPLE %in% amr_samp_c5 & fil_tab$gene_symbol %in% MPNST_pos,]$auto_call)

clust_plot12_rect <- ggplot(scores_sub_PC1_mem[scores_sub_PC1_mem$sample %nin% amr_samp_c2,], aes(x = PC1, y = PC2, colour = clust_code)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = -0.2, lty = 2)

clust_plot12_rect_pop <- ggplot(scores_sub_PC1_mem[scores_sub_PC1_mem$sample %nin% amr_samp_c2,], aes(x = PC1, y = PC2, colour = superPopulation)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = -0.2, lty = 2)

clust_plot12_pop <- ggplot(scores_sub_PC1_mem, aes(x = PC1, y = PC2, colour = superPopulation)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = -0.2, lty = 2)
multiplot(clust_plot12, clust_plot12_pop, clust_plot12_rect, clust_plot12_rect_pop, cols = 2)

##3D plots
source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")
library("gg3D")

make_3Dplot_clust <- function(df_inp, theta=0, phi=0){
  ggplot(df_inp, aes(x=PC1, y=PC2, z=PC3, colour=clust_code)) +
    axes_3D(theta=theta, phi=phi) +
    stat_3D(theta=theta, phi=phi, geom="point") +
    labs_3D(theta=theta, phi=phi, 
            labs=c("PC1", "PC2", "PC3"), 
            angle=c(0,0,0),
            hjust=c(1,1,1), 
            vjust=c(2,2,-1)) +
    theme_void() 
}

clust3D_12 = make_3Dplot_clust(scores_sub_PC1_mem, theta=-90, phi=60)

make_3Dplot_coh <- function(df_inp, theta=0, phi=0){
  ggplot(df_inp, aes(x=PC1, y=PC2, z=PC3, colour=superPopulation)) +
    axes_3D(theta=theta, phi=phi) +
    stat_3D(theta=theta, phi=phi, geom="point") +
    labs_3D(theta=theta, phi=phi, 
            labs=c("PC1", "PC2", "PC3"), 
            angle=c(0,0,0),
            hjust=c(1,1,1), 
            vjust=c(2,2,-1)) +
    theme_void() 
}

coh3D_12 = make_3Dplot_coh(scores_sub_PC1_mem, theta=-90, phi=60)

multiplot(clust3D_12, coh3D_12, cols = 2)



#######
##With EUR stratification

scores_EUR = scores[!is.na(scores$pred.eurPop),]
eur_12 = ggplot(scores_EUR, aes(x = PC1, y = PC2, colour = pred.eurPop)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = -0.2, lty = 2)

all_12 = ggplot(scores_EUR, aes(x = PC1, y = PC2, colour = superPopulation)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = -0.2, lty = 2)

multiplot(eur_12, all_12, cols = 2)

###3D plots

make_3Dplot_eurpop <- function(df_inp, theta=0, phi=0){
  ggplot(df_inp, aes(x=PC1, y=PC2, z=PC3, colour=pred.eurPop)) +
    axes_3D(theta=theta, phi=phi) +
    stat_3D(theta=theta, phi=phi, geom="point") +
    labs_3D(theta=theta, phi=phi, 
            labs=c("PC1", "PC2", "PC3"), 
            angle=c(0,0,0),
            hjust=c(1,1,1), 
            vjust=c(2,2,-1)) +
    theme_void() 
}

eurpop3D_12 = make_3Dplot_eurpop(scores_EUR, theta=180, phi=0)

eurcoh3D_12 = make_3Dplot_coh(scores_EUR, theta=180, phi=0)

multiplot(eurpop3D_12, eurcoh3D_12, cols = 2)


eurpop = unique(scores_EUR$pred.eurPop) ##based on 4 PCs
case_cont_eur_df = list()
for(i in 1:length(unique(eurpop))){
  #print(i)
  #clust_prop = as.data.frame(table(scores_EUR[scores_EUR$pred.eurPop %in% eurpop[i],]$pred.superPop))
  clust_coh_prop = as.data.frame(table(scores_EUR[scores_EUR$pred.eurPop %in% eurpop[i],]$superPopulation))
  cont = sum(clust_coh_prop[clust_coh_prop$Var1 %in% "controls", ]$Freq)
  case = sum(clust_coh_prop[clust_coh_prop$Var1 %in% "sarcoma", ]$Freq)
  
  #print(case)
  #print(cont)
  ratio = (case/1514)/(cont/3116)
  #print(ratio)
  
  case_cont_eur_df[[i]] = cbind.data.frame("case" = case, "control" =  cont, "ratio" = ratio, "eurpop" = eurpop[i])
}

eur_case_cont_df_fin = do.call("rbind.data.frame", case_cont_eur_df)

##Run EDDA with 3-PCs
library(mclust)

score_rect = read.table("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_combset_pca.scores_clustered_rect_sam.tsv", 
           header = T, sep = "\t", stringsAsFactors = F)
#score_rect = score_rect[score_rect$superPopulation %in% c("ISKS", "RISC", "LIONS", "MGRB"),]
##model building
pred_vars_3 = c("PC1", "PC2", "PC3")
pred_vars_4 = c("PC1", "PC2", "PC3", "PC4")
`%nin%` = Negate(`%in%`)

model_superpop_3 = MclustDA(score_rect[score_rect$superPopulation %nin% c("ISKS","MGRB", "RISC", "LIONS"),pred_vars_3], 
                            score_rect$superPopulation[score_rect$superPopulation %nin% c("ISKS","MGRB", "RISC", "LIONS")], modelType = "EDDA")
summary(model_superpop_3)
model_superpop_4 = MclustDA(score_rect[score_rect$superPopulation %nin% c("ISKS","MGRB", "RISC", "LIONS"),pred_vars_4], 
                            score_rect$superPopulation[score_rect$superPopulation %nin% c("ISKS","MGRB", "RISC", "LIONS")], modelType = "EDDA")
summary(model_superpop_4)
plot(model_superpop, "scatterplot")
model_eur_3 = MclustDA(score_rect[score_rect$superPopulation == "EUR",pred_vars_3], score_rect$population[score_rect$superPopulation == "EUR"], modelType = "EDDA")
model_eur_4 = MclustDA(score_rect[score_rect$superPopulation == "EUR",pred_vars_4], score_rect$population[score_rect$superPopulation == "EUR"], modelType = "EDDA")

summary(model_eur_3)
summary(model_eur_4)
#png("Eur_model_EDDA.png", height = 800, width = 800)
#plot(model_eur, "scatterplot")
#dev.off()
pred_mod = function(mod_inp, pred_var, nvar, population){
  temp = predict(mod_inp, score_rect[,pred_var])
#  score_rect$pred.superPop = temp$classification
  te1 = cbind.data.frame(temp$z, "class" = temp$classification)
  colnames(te1) = paste(population, colnames(te1), ".", nvar, sep = "")
  score_rect_mod = cbind(score_rect, te1)
  return(score_rect_mod)
}
score_rect_3 = pred_mod(model_superpop_3, pred_vars_3, 3, "pred.superPop.")
score_rect_4 = pred_mod(model_superpop_4, pred_vars_4, 4) ##same as existing pred.* columns ##not needed

##Eur model with score_rect_3

pred_mod_eur = function(mod_inp, pred_var, nvar, population){
  temp = predict(mod_inp, score_rect_3[,pred_var])
  predclass = as.character(temp$classification)
  predclass[score_rect_3$pred.superPop.class.3 != "EUR"] = NA
  te = temp$z
  te[score_rect_3$pred.superPop.class.3 != "EUR",] = NA
  te1 = cbind.data.frame(te, "class" = predclass)
  colnames(te1) = paste(population, colnames(te1), ".", nvar, sep = "")
  score_rect_mod = cbind(score_rect_3, te1)
  return(score_rect_mod)
}

score_rect_eur_3 = pred_mod_eur(model_eur_3, pred_vars_3, 3, "pred.eurPop.")

##case control ratio

score_rect_eur_3_no1000G = score_rect_eur_3[score_rect_eur_3$superPopulation %in% c("ISKS","MGRB", "RISC", "LIONS"),]
score_rect_eur_3_no1000G = score_rect_eur_3_no1000G[!is.na(score_rect_eur_3_no1000G$pred.eurPop.class.3),]
score_rect_eur_3_no1000G_sub_PC1 = score_rect_eur_3_no1000G[score_rect_eur_3_no1000G$PC1 >= -0.2, ]

dim(score_rect_eur_3_no1000G_sub_PC1)
score_rect_eur_3_no1000G_sub_PC1$superPopulation = ifelse(score_rect_eur_3_no1000G_sub_PC1$superPopulation %in% "MGRB", "controls", "sarcoma")

eurpop3 = as.character(unique(score_rect_eur_3_no1000G_sub_PC1$pred.eurPop.class.3)) ##based on 3 PCs
case_cont_eur3_df = list()
for(i in 1:length(unique(eurpop3))){
  #print(i)
  #clust_prop = as.data.frame(table(scores_EUR[scores_EUR$pred.eurPop %in% eurpop[i],]$pred.superPop))
  clust_coh_prop = as.data.frame(table(score_rect_eur_3_no1000G_sub_PC1[score_rect_eur_3_no1000G_sub_PC1$pred.eurPop.class.3 %in% eurpop3[i],]$superPopulation))
  cont = sum(clust_coh_prop[clust_coh_prop$Var1 %in% "controls", ]$Freq)
  case = sum(clust_coh_prop[clust_coh_prop$Var1 %in% "sarcoma", ]$Freq)
  
  #print(case)
  #print(cont)
  ratio = (case/1523)/(cont/3121)
  #print(ratio)
  
  case_cont_eur3_df[[i]] = cbind.data.frame("case" = case, "control" =  cont, "ratio" = ratio, "eurpop3" = eurpop3[i])
}

eur3_case_cont_df_fin = do.call("rbind.data.frame", case_cont_eur3_df)

##plots

eurpop3D_12 = make_3Dplot_eurpop(score_rect_eur_3_no1000G_sub_PC1[score_rect_eur_3_no1000G_sub_PC1$pred.eurPop.class.3], theta=180, phi=0)

eurcoh3D_12 = make_3Dplot_coh(scores_EUR, theta=180, phi=0)

multiplot(eurpop3D_12, eurcoh3D_12, cols = 2)

##modelling with tsne

library(Rtsne)
