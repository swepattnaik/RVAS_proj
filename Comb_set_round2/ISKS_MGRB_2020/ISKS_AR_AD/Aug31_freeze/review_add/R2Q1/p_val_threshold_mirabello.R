

ISKS_ASRB_sub_old = read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/isks_cyto_ppifisher_ASRB_comb_genes_top_new_rand_comp_score_Aug31_sep19.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)

ISKS_ASRB_sub = read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/isks_cyto_ppifisher_ASRB_comb_genes_totex_2021_subPC1.tsv",
           sep = "\t", header = T, stringsAsFactors = F)
ppi_res_fil_final_comb <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/ppi_res_fil_final_SKATbin_comb_str_biog_totex_subPC1.tsv", sep = "\t", header = T,
                                     stringsAsFactors = F)
ppi_res_fil_final_comb_asrb_sub = ppi_res_fil_final_comb[ppi_res_fil_final_comb$gene %in% ISKS_ASRB_sub$gene,]
ppi_res_fil_final_comb_asrb_sub = ppi_res_fil_final_comb_asrb_sub[ppi_res_fil_final_comb_asrb_sub$MGRB <= 14,]

Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")
Shelterin_extn = c("ACD", "POT1", "TERF1", "TERF2", "TERF2IP", "TINF2", "ATM", 
                   "BAG3", "BLM", "BRCA1", "CALD1", "CLK3", "DCLRE1B", "FANCD2", 
                   "FBXO4", "HSPA4", "KIAA1191", "MRE11A", "NBN", "PINX1", "PRKDC", 
                   "RAD50", "SLX4", "STUB1", "TNKS", "TNKS2", "U2AF2", "UCHL1", 
                   "WRN", "XRCC5", "XRCC6", "SMARCAL1", "STAG3")
CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
mirabello_genes <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/mirabello_genes.txt",
                              sep = "", header = F, stringsAsFactors = F)
##remove trailing and leading white spaces
library(dplyr)
library(stringr)
mira_genes <- mirabello_genes %>% 
  mutate(V1 = str_trim(mirabello_genes$V1, side = "both"))


ISKS_ASRB_sub_old[ISKS_ASRB_sub_old$gene %in% mira_genes$V1 & ISKS_ASRB_sub_old$pval_SKATbin < 0.1,c(3,13:15)]

SKAT_aug31 = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/isks_PC1234_ver4_clinrect_Aug31.rds")
##filtered with positive weight difference
#SKAT_aug31_wt_diff = SKAT_aug31[SKAT_aug31$wt_diff > 0 & SKAT_aug31$MGRB <= 14 & SKAT_aug31$ISKS > 1,]
SKAT_aug31_wt_diff = SKAT_aug31[SKAT_aug31$wt_diff > 0 & SKAT_aug31$ISKS > 1,]
SKAT_aug31_wt_neg = SKAT_aug31[SKAT_aug31$wt_diff < 0 & SKAT_aug31$ISKS > 1,]

dim(SKAT_aug31_wt_diff)
dim(SKAT_aug31_wt_diff[SKAT_aug31_wt_diff$gene %in% mira_genes$V1,])
par(mfrow=c(1,2))
hist(SKAT_aug31[SKAT_aug31$gene %in% mira_genes$V1 & SKAT_aug31$ISKS > 1 & SKAT_aug31$MGRB > 1,]$pval_SKATbin, breaks = 20, 
     ylab = "Number of Mirabello genes",
     xlab = "P-value Distribution")

hist(SKAT_aug31_wt_diff[SKAT_aug31_wt_diff$gene %in% mira_genes$V1 & SKAT_aug31_wt_diff$ISKS > 1 & SKAT_aug31_wt_diff$MGRB > 1,]$pval_SKATbin, breaks = 10, 
     ylab = "Number of Mirabello genes",
     xlab = "P-value Distribution")

hist(SKAT_aug31_wt_neg[SKAT_aug31_wt_neg$gene %in% mira_genes$V1 & SKAT_aug31_wt_neg$ISKS > 1 & SKAT_aug31_wt_neg$MGRB > 1,]$pval_SKATbin, breaks = 10, 
     ylab = "Number of Mirabello genes",
     xlab = "P-value Distribution")

hist(SKAT_aug31[SKAT_aug31$gene %in% mira_genes$V1 & SKAT_aug31$ISKS > 1 & SKAT_aug31$MGRB > 1,]$pval_SKATbin, breaks = 10, 
     ylab = "Number of Mirabello genes",
     xlab = "P-value Distribution")

library(ggplot2)
gg <- ggplot()
df_05 = SKAT_aug31_wt_diff[SKAT_aug31_wt_diff$gene %in% mira_genes$V1 & SKAT_aug31_wt_diff$ISKS > 1 ,c(3,14)]
df_05$group = "A"
df_1 = SKAT_aug31_wt_diff[SKAT_aug31_wt_diff$gene %in% mira_genes$V1 & SKAT_aug31_wt_diff$ISKS > 1 ,c(3,14)]
df_1$group = "B"
df = rbind.data.frame(df_05, df_1)
##Alternate representation
ggplot(df_05, aes(pval_SKATbin)) + geom_histogram(binwidth=0.5)
ggplot(df_05, aes(pval_SKATbin)) +  geom_bar() + scale_x_binned()
source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")
pp_05 = ggplot(df_05, aes(pval_SKATbin)) +  geom_bar() + scale_x_binned(n.breaks = 20, limits = c(0,1)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 10, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10, face = "bold"))
pp_05 = pp_05 + ylab("Gene counts") + xlab("P-value") + theme()
pp_1 = ggplot(df_1, aes(pval_SKATbin)) +  geom_bar() + scale_x_binned(n.breaks = 10, limits = c(0,1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pp_1 = pp_1 + xlab("P-value")
multiplot(pp_05, pp_1, cols = 2)
pp_05

mira_05 = ISKS_ASRB_sub_old[ISKS_ASRB_sub_old$gene %in% mira_genes$V1 & ISKS_ASRB_sub_old$pval_SKATbin < 0.05 ,c(3,13:15)]
SKAT_05 = SKAT_aug31_wt_diff[SKAT_aug31_wt_diff$pval_SKATbin < 0.05 ,]
mira_05_1 = ISKS_ASRB_sub_old[ISKS_ASRB_sub_old$gene %in% mira_genes$V1 & ISKS_ASRB_sub_old$pval_SKATbin < 0.1 
                              & ISKS_ASRB_sub_old$pval_SKATbin > 0.05 ,c(3,13:15)]
SKAT_05_1 = SKAT_aug31_wt_diff[SKAT_aug31_wt_diff$pval_SKATbin < 0.1 & SKAT_aug31_wt_diff$pval_SKATbin > 0.05 ,]

mira_0_1 = SKAT_aug31_wt_diff[SKAT_aug31_wt_diff$gene %in% mira_genes$V1 & SKAT_aug31_wt_diff$pval_SKATbin < 0.1,]
SKAT_0_1 = SKAT_aug31_wt_diff[SKAT_aug31_wt_diff$pval_SKATbin < 0.1,]

