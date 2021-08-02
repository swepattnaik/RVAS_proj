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


p_Data <- read.table("~/RVAS/comb_set_2020/pop_PCA/no1000G/MGRB_ISKS_no1000G_pca_rnd2.scores.tsv", 
                     header = T, sep = "\t", stringsAsFactors = F)
# p_Data_1000G <- p_Data[p_Data$superPopulation %nin% c("ISKS", "RISC", "LIONS", "MGRB"),]
# write.table(p_Data_1000G$sample, "~/RVAS/comb_set_2020/pop_PCA/samp_1000G.txt",
#             row.names = F, quote = F, sep = "")
p_Data$superPopulation <- ifelse(grepl("^[ABZ]", p_Data$sample), "MGRB", "ISKS")
#p_Data <- p_Data[p_Data$superPopulation %in% c("ISKS", "RISC", "LIONS", "MGRB"),]
#p_Data <- p_Data[p_Data$superPopulation %in% c("ISKS"),]
p_Data_MGRB <-  p_Data[p_Data$superPopulation %in% c("MGRB"),]

#remove duplicates
dup_samp <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/fin_samp_dup_drop.txt", sep = "", header = T,
                       stringsAsFactors = F)
dup_samp <- dup_samp$x[grepl("^[ABZ]", dup_samp$x)]
p_Data_MGRB <- p_Data_MGRB[p_Data_MGRB$sample %nin% dup_samp,]
p_Data_MGRB$rect_sam <- p_Data_MGRB$sample
p_Data_MGRB <- p_Data_MGRB[p_Data_MGRB$rect_sam %in% QC2_dat_pass$new_sampleid,]
##rename p_Data
#p_Data_ISKS <-  p_Data[p_Data$superPopulation %in% c("ISKS", "RISC", "LIONS"),]
p_Data_ISKS <-  p_Data[p_Data$superPopulation %in% c("ISKS"),]
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
# write.table(p_Data_noCH, "~/RVAS/comb_set_2020/pop_PCA/rect_sam_MGRB_ISKS_no1000G_combset_pca.scores_no1000G_clustered.tsv", sep = "\t", row.names = F,
#             quote = F)
scores_no1000G <- p_Data_noCH

##Add 1000G predicted superPOP cluster
superpop_pred <- read.delim("~/RVAS/comb_set_2020/pop_PCA/rect_sam_MGRB_ISKS_no1000G_combset_pca.scores_clustered.tsv",
                            sep = "\t", header = T, stringsAsFactors = F)
scores_no1000G$pred.superPop <- superpop_pred[match(scores_no1000G$rect_sam, superpop_pred$rect_sam), 39]

#scores_no1000G$Country <- comb_pheno[match(scores_no1000G$sample, comb_pheno$pmn),120]
scores_no1000G$Country <- comb_pheno[match(scores_no1000G$rect_sam, comb_pheno$pmn),64]
scores_no1000G$Country <- ifelse(is.na(scores_no1000G$Country), "Aus_MGRB", scores_no1000G$Country)
scores_no1000G$Country <- ifelse((!grepl("^[ABZ]",scores_no1000G$sample) & (scores_no1000G$Country %in% "Aus_MGRB")), "Unknown", scores_no1000G$Country)

ggplot(scores_no1000G, aes(x = PC1, y = PC2, colour = Country)) + geom_point() + theme_bw()
g12 <- ggplot(scores_no1000G, aes(x = PC1, y = PC2, colour = Country)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = 0.1, lty = 2) + geom_hline(yintercept = -0.16, lty = 2)
g13 <- ggplot(scores_no1000G, aes(x = PC1, y = PC3, colour = Country)) + geom_point() + theme_bw() +
  theme(legend.position="bottom", legend.title = element_blank()) 
g23 <- ggplot(scores_no1000G, aes(x = PC2, y = PC3, colour = Country)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = -0.16, lty = 2) + geom_hline(yintercept = -0.05, lty = 2)


source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")
multiplot(g12, g13, g23, cols = 1)

##factor: predicted superPopulation
ggplot(scores_no1000G, aes(x = PC1, y = PC2, colour = pred.superPop)) + 
  geom_point() + theme_bw() + stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(level = 0.95) +
  coord_fixed()

p12 <- ggplot(scores_no1000G, aes(x = PC1, y = PC2, colour = pred.superPop)) + geom_point() + theme_bw() + 
  stat_ellipse(type = "norm", linetype = 2) + theme(legend.position="bottom") +
  theme(legend.title = element_blank()) + 
  stat_ellipse(level = 0.95)
p13 <- ggplot(scores_no1000G, aes(x = PC1, y = PC3, colour = pred.superPop)) + geom_point() + theme_bw() + 
  stat_ellipse(type = "norm", linetype = 2) + theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  stat_ellipse(level = 0.95)
p23 <- ggplot(scores_no1000G, aes(x = PC2, y = PC3, colour = pred.superPop)) + geom_point() + theme_bw() + 
  stat_ellipse(type = "norm", linetype = 2) + theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  stat_ellipse(level = 0.95) 

multiplot(p12, p13, p23, cols = 3)

##factor: MGRB or ISKS
k12 <- ggplot(scores_no1000G, aes(x = PC1, y = PC2, colour = superPopulation)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = 0.1, lty = 2) + geom_hline(yintercept = -0.16, lty = 2)
k13 <- ggplot(scores_no1000G, aes(x = PC1, y = PC3, colour = superPopulation)) + geom_point() + theme_bw() +
  theme(legend.position="bottom", legend.title = element_blank()) 
k23 <- ggplot(scores_no1000G, aes(x = PC2, y = PC3, colour = superPopulation)) + geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = -0.16, lty = 2) + geom_hline(yintercept = -0.05, lty = 2)

#multiplot(k12, k13, k23, cols = 3)

##combine plots
#multiplot(p12, k12, p13,  k13, p23, k23, cols = 3)

#multiplot(p12, k12, g12, p13, k13, g13, p23, k23, g23, cols = 3) 
multiplot(g12, k12, p12, g13, k13, p13, g23, k23, p23, cols = 3) 

###Proportion of variance explained
##Without 1000G
eigenval_no1000G <- read.delim("~/RVAS/comb_set_2020/pop_PCA/no1000G/MGRB_ISKS_no1000G_pca_rnd2.evals.txt",
                               sep = ":", header = F, stringsAsFactors = F)
eigenval_no1000G$V1 <- gsub("^.*PC", "", eigenval_no1000G$V1)
eigenval_no1000G$V1 <- as.numeric(gsub("'", "", eigenval_no1000G$V1))
eigenval_no1000G$V2 <- gsub("}", "", eigenval_no1000G$V2)
eigenval_no1000G$V2 <- as.numeric(gsub(",", "", eigenval_no1000G$V2))
eigenval_no1000G <- eigenval_no1000G[order(eigenval_no1000G$V1, decreasing = F),]
eigenval_no1000G$prop_var <- eigenval_no1000G$V2/sum(eigenval_no1000G$V2)
colnames(eigenval_no1000G)[1:2] <- c("PC", "Eigen_val") 

ggplot(eigenval_no1000G, aes(PC, prop_var)) +
  geom_bar(stat = "Identity") +
  xlab("Principal Components") +
  ylab("Variance Explained")

##samples that are not covered by MGRB
scores_sub_PC12 = scores_no1000G[scores_no1000G$PC2 < -0.16 & scores_no1000G$PC1 > 0.1,]

scores_sub_PC23 = scores_no1000G[scores_no1000G$PC2 < -0.16 & scores_no1000G$PC3 > -0.05,]
##Variant file used for SKAT
fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv", 
                      sep = "\t", header = T, stringsAsFactors = F)

genes_outlier <- unique(fil_tab[fil_tab$SAMPLE %in% scores_sub_PC23$rect_sam,]$gene_symbol)

Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")
CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
MPNST_pos <- c("NF1", "LZTR1", "SDHA", "SDHB", "SDHC", "SDHD")
Breast_cancer_genes <- c("BRCA1", "BRCA2", "PALB2", "RAD51C", "CDH1")
table(genes_outlier %in% Shelterin)
genes_outlier[genes_outlier %in% Shelterin]  #"TIMELESS" (3; C3), "TERF2" (2; C3); All C3
genes_outlier[genes_outlier %in% CEP_HAUS_core] ##CEP72 (1; C4)
genes_outlier[genes_outlier %in% MPNST_pos] #"NF1"(1; C4), LZTR1 (2; C3)
genes_outlier[genes_outlier %in% "TP53"] #0
genes_outlier[genes_outlier %in% Breast_cancer_genes] #"RAD51C"(1; C3)

fil_tab[fil_tab$SAMPLE %in% scores_sub_PC23$rect_sam & fil_tab$gene_symbol %in% c("NF1", "LZTR1"),c(1:9,11,127:130)] #1833
fil_tab[fil_tab$SAMPLE %in% scores_sub_PC23$rect_sam & fil_tab$gene_symbol %in% c("TIMELESS", "TERF2"),c(1:9,11,127:130)]
fil_tab[fil_tab$SAMPLE %in% scores_sub_PC23$rect_sam & fil_tab$gene_symbol %in% "CEP72", c(1:9,11,127:130)]
fil_tab[fil_tab$SAMPLE %in% scores_sub_PC23$rect_sam & fil_tab$gene_symbol %in% "RAD51C", c(1:9,11,127:130)]

##complex diagnostic

source("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/review_add/R2Q6_popstruc/complex_diagnostic.R")
sens_list <- list()
sens_list[[1]] <- cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 0, 0, 0, 0, "ISKSvsMGRB")
sens_list[[2]] <- cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 5, 0, 43, 0, "ISKSvsMGRB")
#cpx_OR_fisher_one(isks_mgrb_genes_noC3, 1644, 3205, Shelterin, 0, 0, 0, 0, "ISKSvsMGRBnoC3") ##no change as no C4/5 variant were lost
#cpx_OR_fisher_one(isks_mgrb_genes_noC3, 1644, 3205, Shelterin, 0, 0, 0, 0, "ISKSvsMGRBnoC3")
sens_list[[3]] <- cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 0, 0, 0, 0, "ISKSvsMGRB")
sens_list[[4]] <- cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 1, 0, 43, 0, "ISKSvsMGRB")
sens_list[[5]] <- cpx_OR_fisher_one(isks_mgrb_genes_noC3, 1644, 3205, CEP_HAUS_core, 0, 0, 0, 0, "ISKSvsMGRBnoC3")
sens_list[[6]] <- cpx_OR_fisher_one(isks_mgrb_genes_noC3, 1644, 3205, CEP_HAUS_core, 1, 0, 43, 0, "ISKSvsMGRBnoC3")

##overlap with MPNST group (scores_sub_PC23$rect_sam[scores_sub_PC23$rect_sam %in% samp_MPNST])
##from MPNST folder
##3 MPNST samples have been excluded: "1833" "2664" "2792"
sens_list[[7]] <- cpx_OR_fisher_one(isks_mpnst_mgrb_genes, 157, 3205, MPNST_pos, 0, 0, 0, 0, "MPNSTvsMGRB")
sens_list[[8]] <- cpx_OR_fisher_one(isks_mpnst_mgrb_genes, 157, 3205, MPNST_pos, 3, 0, 3, 0, "MPNSTvsMGRB")
sens_list[[9]] <- cpx_OR_fisher_one(isks_mpnst_mgrb_genes_noC3, 157, 3205, MPNST_pos, 0, 0, 0, 0, "MPNSTvsMGRB_noC3")
sens_list[[10]] <- cpx_OR_fisher_one(isks_mpnst_mgrb_genes_noC3, 157, 3205, MPNST_pos, 1, 0, 3, 0, "MPNSTvsMGRB_noC3")

##
sens_list[[11]] <- cpx_OR_fisher_one(isks_mpnst_mgrb_genes, 157, 3205, CEP_HAUS_core, 0, 0, 0, 0, "MPNSTvsMGRB")
sens_list[[12]] <- cpx_OR_fisher_one(isks_mpnst_mgrb_genes, 157, 3205, CEP_HAUS_core, 1, 0, 3, 0, "MPNSTvsMGRB")
sens_list[[13]] <- cpx_OR_fisher_one(isks_mpnst_mgrb_genes_noC3, 157, 3205, CEP_HAUS_core, 0, 0, 0, 0, "MPNSTvsMGRB_noC3")
sens_list[[14]] <- cpx_OR_fisher_one(isks_mpnst_mgrb_genes_noC3, 157, 3205, CEP_HAUS_core, 1, 0, 3, 0, "MPNSTvsMGRB_noC3")

sens_list_df <- do.call("rbind.data.frame", sens_list)
write.table(sens_list_df, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/review_add/R2Q6_popstruc/Results/complex_sensitivity.tsv",
            row.names = F, sep = "\t", quote = F)
#View(sens_list_df)
