##WKS test
##weighted test for the entire SKAT output

##SKAT output
Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_skat_para_result_isks_combset2020_uni_MAF_PC1234_ver4_clinrect_Aug31.rds")
df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})
Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_para_pc1234_SKAT_Enriched_ISKS_2020_Aug31.rds")
##Add p-value
Exome_pc123_srt_SKAT_case_enr_nCH_pval <- Map(cbind.data.frame, Exome_pc123_srt_SKAT_case_enr_nCH, Exome_pc123_srt_SKAT)
##save as rds file
SKAT_all_PC4_gender <- Exome_pc123_srt_SKAT_case_enr_nCH_pval[[4]]
saveRDS(SKAT_all_PC4_gender, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/SKAT_all_PC4_gender.rds", compress = T)
#SKAT_all_PC4_gender <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/SKAT_all_PC4_gender.rds")
#SKAT_all_PC4_gender <- SKAT_all_PC4_gender[,c(1:6,14)]
SKAT_all_PC4_gender_ranked <- SKAT_all_PC4_gender[order(SKAT_all_PC4_gender$wt_diff, decreasing = T),]
SKAT_all_PC4_gender_ranked <- SKAT_all_PC4_gender_ranked[,c(1:6,14)]

##run Weighted Kolmogorov Smirnov test
source("~/WKS/wks/wks.r")
##geneset = C2 from MsigDB
C2 <- readRDS("~/WKS/wks/C2.rds")
##weighted gene list (weighted using wt_diff)
gen_wt_diff <- SKAT_all_PC4_gender_ranked$wt_diff
names(gen_wt_diff) <- SKAT_all_PC4_gender_ranked$gene
pv_wt_diff_wks <- WKS.test(gen_wt_diff,C2, alternative = "greater")

gen_logSKATbin <- -log10(SKAT_all_PC4_gender_ranked$pval_SKATbin)
names(gen_logSKATbin) <- SKAT_all_PC4_gender_ranked$gene
pv_skatbin_wks <- WKS.test(gen_logSKATbin,C2, alternative = "greater")

gen_rank_score <- -log10(SKAT_all_PC4_gender_ranked$pval_SKATbin) * SKAT_all_PC4_gender_ranked$wt_diff
names(gen_rank_score) <- SKAT_all_PC4_gender_ranked$gene
pv_rank_score_wks <- WKS.test(gen_rank_score,C2, alternative = "greater")

##cpx_list from paper_figs_tab.Rmd
pv_rank_score_wks <- WKS.test(gen_rank_score,cpx_list, alternative = "greater")


##iDEA test data input generation
#devtools::install_github('xzhoulab/iDEA')
#https://xzhoulab.github.io/iDEA/
library(iDEA)
data(humanGeneSets)
SKAT_all_PC4_gender <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/SKAT_all_PC4_gender.rds")

zscore <- qnorm(SKAT_all_PC4_gender$pval_SKATbin/2.0, lower.tail=FALSE)
beta <- (SKAT_all_PC4_gender$wt_diff - mean(SKAT_all_PC4_gender$wt_diff))/sd(SKAT_all_PC4_gender$wt_diff) ## Cohen's effect size
se_beta <- abs(beta/zscore) ## to approximate the standard error of beta
beta_var = se_beta^2  ### square 
summary = data.frame(beta = beta,beta_var = beta_var)
## add the gene names as the rownames of summary
rownames(summary) = (SKAT_all_PC4_gender$gene) ### or the gene id column in the res_DE results

##Geneset file
Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")

Glutamate_NMDA <- c("DLG1", "GRIN2A", "GRIA1", "CAM2KB", "CAM2KD")

ZBTB16_complex <- c("ZBTB16", "ATG7", "ASB10", "ASB8", "FBXO7")

CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")

CENP_complex <- c("CENPC", "CENPF", "BUB1", "BUB1B", "AURKB", "RANGAP1", "RACGAP1", "MAD2L2")

Sarcoma_genes <- c("TP53", "SDHA", "SDHB", "SDHD")

TP53_control <- c("TP53") #positive control

MPNST_pos <- c("NF1", "LZTR1", "SDHA", "SDHB", "SDHD")

Breast_cancer_genes <- c("BRCA1", "BRCA2", "PALB2")

cpx_list <- list(Shelterin, Glutamate_NMDA, ZBTB16_complex, CEP_HAUS_core, CENP_complex, Breast_cancer_genes, Sarcoma_genes,
                 TP53_control, MPNST_pos)
names(cpx_list) <- c("Shelterin geneset", "Glutamate_NMDA", "ZBTB16_complex", "Centrosome geneset", "CENP_complex", 
                     "BRCA geneset", "Sarcoma geneset", "TP53", "NF1 geneset")

gene_annot_list <- lapply(cpx_list, function(x)ifelse(rownames(summary) %in% x, 1, 0))
gene_annot_list_df <- do.call("cbind.data.frame", gene_annot_list)
rownames(gene_annot_list_df) <- rownames(summary)

##create idea object
idea <- CreateiDEAObject(summary, gene_annot_list_df, max_var_beta = 15000, min_precent_annot = 0.0000025, num_core=10)

idea_h_gset <- CreateiDEAObject(summary, humanGeneSets, max_var_beta = 100, min_precent_annot = 0.0025, num_core=15)
##Fit the model EM-MCMC
idea <- iDEA.fit(idea,
                 fit_noGS=FALSE,
                 init_beta=NULL, 
                 init_tau=c(-2,0.5),
                 min_degene=1,
                 em_iter=15,
                 mcmc_iter=1000, 
                 fit.tol=1e-5,
                 modelVariant = F,
                 verbose=TRUE)

idea_h_gset <- iDEA.fit(idea_h_gset,
                 fit_noGS=FALSE,
                 init_beta=NULL, 
                 init_tau=c(-2,0.5),
                 min_degene=5,
                 em_iter=15,
                 mcmc_iter=1000, 
                 fit.tol=1e-5,
                 modelVariant = F,
                 verbose=TRUE)
##Correct p-values

idea <- iDEA.louis(idea)
head(idea@gsea)
##Bayesian model averaging (BMA)
idea <- iDEA.BMA(idea)
head(idea@BMA_pip)
