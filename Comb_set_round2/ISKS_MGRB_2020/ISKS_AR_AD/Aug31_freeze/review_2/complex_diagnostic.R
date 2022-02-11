##For review 2
##Contains all complex related analysis (Fig.3a, Figure S7, Table S13)
##Note, table S13 needs rectification


Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")

CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
#CEP_HAUS_core_mod <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1")

MPNST_pos <- c("NF1", "LZTR1", "SDHA", "SDHB", "SDHC", "SDHD")

##ISKS vs MGRB
isks_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_para_tab_ISKS_MGRB_C345_Aug31.rds")
isks_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Exome_para_tab_ISKS_MGRB_minusC3_Aug31.rds")

##ISKS vs MGRB (minus_13)
isks_mgrb_genes_m13 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Exome_para_tab_ISKS_MGRB_C345_Aug31_sub_PC1.rds")
isks_mgrb_genes_noC3_m13 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Exome_para_tab_ISKS_MGRB_minusC3_Aug31_sub_PC1.rds")

##ISKS vs MGRB (minus_13_cc)
isks_mgrb_genes_cc <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_2/Exome_para_tab_ISKS_MGRB_C345_Jan10_ccrect.rds")
isks_mgrb_genes_noC3_cc <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_2/Exome_para_tab_ISKS_MGRB_minusC3_Jan10_ccrect.rds")


##ISKS_MPNST vs MGRB
isks_mpnst_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/MPNST/Exome_para_ISKS_MPNST_sep30.rds")
isks_mpnst_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/MPNST/Exome_para_ISKS_MPNST_sep30_noC3.rds")



cpx_OR_fisher_one <- function(ppi_res,case_coh_size, cont_coh_size, cpx_list, sub_case, sub_cont, sub_ISKS, sub_MGRB, coh){
  
  ppi_res_tab <- ppi_res[ppi_res$gene %in% cpx_list,]
  # ppi_res_tab[,2] <- ifelse(ppi_res_tab[,2] == 0, 1, ppi_res_tab[,2])
 # inp <- c(sum(ppi_res_tab[,1]) - sub_case , case_coh_size - sum(ppi_res_tab[,1]) - sub_ISKS , 
 #          sum(ppi_res_tab[,2]) - sub_cont, cont_coh_size - sum(ppi_res_tab[,2]) - sub_MGRB)
  inp <- c(sum(ppi_res_tab[,1]) - sub_case , case_coh_size - sub_ISKS - (sum(ppi_res_tab[,1]) - sub_case), 
           sum(ppi_res_tab[,2]) - sub_cont, cont_coh_size - sub_MGRB - (sum(ppi_res_tab[,2]) - sub_cont))
  sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
  colnames(sim_mat) <- c("case", "cont")
  rownames(sim_mat) <- c("hits", "no_hits")
  #ft <- fisher.test(sim_mat, alternative = "greater")
  #ft <- fisher.test(sim_mat)
  cpx_name = deparse(substitute(cpx_list))
  ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.95)
  ft_df <- cbind.data.frame("gene" = cpx_name ,"Cases" = sim_mat[1,1],
                                 "Controls" = sim_mat[1,2],
                                 "Fish_pval" = ft$p.value,"CI_lower" = ft$conf.int[1],
                                 "CI_upper" = ft$conf.int[2],
                                 "OR_Fish" = ft$estimate, "case_coh_size" = sum(sim_mat[,1]),
                                "cont_coh_size" = sum(sim_mat[,2]),
                                 "Coh" = coh)
  return(ft_df)
}

cpx_res_list = list()
list_names = c("original", "PCs_min13", "PCs_EMMAX_min13", "kingGRM_EMMAX_min13", "kingPCs_EMMAX_min13")

######CEP HAUS complex
##with C345 variants
# cpx_res_list[[1]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 0 , 0 , 0, 0, "original")
# ##totex + rm_13AFR (SSNA1 p-val < 0.1, 3 ISKS lost)
# cpx_res_list[[2]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 3 , 0 , 13, 0, "PCs_min13")
# ##totex + rm_13AFR + ancestry PCs + ancestry adjusted kinship matrix (CEP72 p-val < 0.1, 7 ISKS + 5 MGRB)
# cpx_res_list[[3]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 7 , 5 , 13, 0, "PCs_EMMAX_min13")
# ##totex + rm_13AFR + GRM (HAUS4 & SSNA1 p-val < 0.1, 6 ISKS)
# cpx_res_list[[4]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 6 , 0 , 13, 0, "kingGRM_EMMAX_min13")
# ##totex + rm_13AFR + kinship  PCs + ancestry adjusted kinship matrix (HAUS4 & SSNA1 p-val < 0.1, 6 ISKS)
# cpx_res_list[[5]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 6 , 0 , 13, 0, "kingPCs_EMMAX_min13")
# 
# tab_C345_CEP_HAUS = do.call("rbind.data.frame", cpx_res_list)
# path = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/review_2/Results/"
# write.table(tab_C345_CEP_HAUS, paste(path, "tab_fisher_C345_CEP_HAUS.tsv", sep = ""), sep = "\t",
#             row.names = F, quote = F)

##noC3 variants
path = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/review_2/Results/"
cpx_res_list_noC3 = list()
cpx_res_list_noC3[[1]] = cpx_OR_fisher_one(isks_mgrb_genes_noC3, 1644, 3205, CEP_HAUS_core, 0 , 0 , 0, 0, "original")
##totex + rm_13AFR (SSNA1 p-val < 0.1, 3 ISKS lost)
cpx_res_list_noC3[[2]] = cpx_OR_fisher_one(isks_mgrb_genes_noC3, 1644, 3205, CEP_HAUS_core, 0 , 0 , 13, 0, "PCs_min13")
##totex + rm_13AFR + ancestry PCs + ancestry adjusted kinship matrix (CEP72 p-val < 0.1, 6 ISKS + 2 MGRB)
cpx_res_list_noC3[[3]] = cpx_OR_fisher_one(isks_mgrb_genes_noC3, 1644, 3205, CEP_HAUS_core, 6 , 2 , 13, 0, "PCs_EMMAX_min13")
##totex + rm_13AFR + GRM (HAUS4 p-val < 0.1, 2 ISKS)
cpx_res_list_noC3[[4]] = cpx_OR_fisher_one(isks_mgrb_genes_noC3, 1644, 3205, CEP_HAUS_core, 2 , 0 , 13, 0, "kingGRM_EMMAX_min13")
##totex + rm_13AFR + kinship  PCs + ancestry adjusted kinship matrix (HAUS4 < 0.1, 2 ISKS)
cpx_res_list_noC3[[5]] = cpx_OR_fisher_one(isks_mgrb_genes_noC3, 1644, 3205, CEP_HAUS_core, 2 , 0 , 13, 0, "kingPCs_EMMAX_min13")

tab_noC3_CEP_HAUS = do.call("rbind.data.frame", cpx_res_list_noC3)
write.table(tab_noC3_CEP_HAUS, paste(path, "tab_noC3_CEP_HAUS.tsv", sep = ""), sep = "\t",
            row.names = F, quote = F)

##Shelterin Complex
# cpx_res_list_Shel = list()
# ##with C345 variants
# cpx_res_list_Shel[[1]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 0 , 0 , 0, 0, "original")
# ##totex + rm_13AFR (SSNA1 p-val < 0.1, 3 ISKS lost)
# cpx_res_list_Shel[[2]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 0 , 0 , 13, 0, "PCs_min13")
# ##totex + rm_13AFR + ancestry PCs + ancestry adjusted kinship matrix (CEP72 p-val < 0.1, 7 ISKS + 5 MGRB)
# cpx_res_list_Shel[[3]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 0 , 5 , 13, 0, "PCs_EMMAX_min13")
# ##totex + rm_13AFR + GRM (HAUS4 & SSNA1 p-val < 0.1, 6 ISKS)
# cpx_res_list_Shel[[4]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 0 , 0 , 13, 0, "kingGRM_EMMAX_min13")
# ##totex + rm_13AFR + kinship  PCs + ancestry adjusted kinship matrix (HAUS4 & SSNA1 p-val < 0.1, 6 ISKS)
# cpx_res_list_Shel[[5]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 0 , 0 , 13, 0, "kingPCs_EMMAX_min13")

# tab_C345_SHEL = do.call("rbind.data.frame", cpx_res_list_Shel)
# write.table(tab_C345_SHEL, paste(path, "tab_C345_SHEL.tsv", sep = ""), sep = "\t",
#             row.names = F, quote = F)

##noC3 variants
cpx_res_list_Shel_noC3 = list()
cpx_res_list_Shel_noC3[[1]] = cpx_OR_fisher_one(isks_mgrb_genes_noC3, 1644, 3205, Shelterin, 0 , 0 , 0, 0, "original")
##totex + rm_13AFR (SSNA1 p-val < 0.1, 3 ISKS lost)
cpx_res_list_Shel_noC3[[2]] = cpx_OR_fisher_one(isks_mgrb_genes_noC3, 1644, 3205, Shelterin, 0 , 0 , 13, 0, "PCs_min13")
##totex + rm_13AFR + ancestry PCs + ancestry adjusted kinship matrix (CEP72 p-val < 0.1, 6 ISKS + 2 MGRB)
cpx_res_list_Shel_noC3[[3]] = cpx_OR_fisher_one(isks_mgrb_genes_noC3, 1644, 3205, Shelterin, 0 , 0 , 13, 0, "PCs_EMMAX_min13")
##totex + rm_13AFR + GRM (HAUS4 p-val < 0.1, 2 ISKS)
cpx_res_list_Shel_noC3[[4]] = cpx_OR_fisher_one(isks_mgrb_genes_noC3, 1644, 3205, Shelterin, 0 , 0 , 13, 0, "kingGRM_EMMAX_min13")
##totex + rm_13AFR + kinship  PCs + ancestry adjusted kinship matrix (HAUS4 < 0.1, 2 ISKS)
cpx_res_list_Shel_noC3[[5]] = cpx_OR_fisher_one(isks_mgrb_genes_noC3, 1644, 3205, Shelterin, 0 , 0 , 13, 0, "kingPCs_EMMAX_min13")

tab_noC3_SHEL = do.call("rbind.data.frame", cpx_res_list_Shel_noC3)
write.table(tab_noC3_SHEL, paste(path, "tab_noC3_SHEL.tsv", sep = ""), sep = "\t",
            row.names = F, quote = F)

##CEP HAUS in MPNST subset
##get MPNST samples and check if there are any 
##noC3
# cpx_res_list_noC3_mpnst = list()
# cpx_res_list_noC3_mpnst[[1]] = cpx_OR_fisher_one(isks_mpnst_mgrb_genes_noC3, 157, 3205, CEP_HAUS_core, 0 , 0 , 0, 0, "original")
# ##totex + rm_13AFR (SSNA1 p-val < 0.1, 3 ISKS lost)
# cpx_res_list_noC3_mpnst[[2]] = cpx_OR_fisher_one(isks_mpnst_mgrb_genes_noC3, 157, 3205, CEP_HAUS_core, 0 , 0 , 2, 0, "PCs_min13")
# ##totex + rm_13AFR + ancestry PCs + ancestry adjusted kinship matrix (CEP72 p-val < 0.1, 6 ISKS + 2 MGRB)
# cpx_res_list_noC3_mpnst[[3]] = cpx_OR_fisher_one(isks_mpnst_mgrb_genes_noC3, 157, 3205, CEP_HAUS_core, 2 , 2 , 2, 0, "PCs_EMMAX_min13")
# ##totex + rm_13AFR + GRM (HAUS4 p-val < 0.1, 2 ISKS)
# cpx_res_list_noC3_mpnst[[4]] = cpx_OR_fisher_one(isks_mpnst_mgrb_genes_noC3, 157, 3205, CEP_HAUS_core, 1 , 0 , 2, 0, "kingGRM_EMMAX_min13")
# ##totex + rm_13AFR + kinship  PCs + ancestry adjusted kinship matrix (HAUS4 < 0.1, 2 ISKS)
# cpx_res_list_noC3_mpnst[[5]] = cpx_OR_fisher_one(isks_mpnst_mgrb_genes_noC3, 157, 3205, CEP_HAUS_core, 1 , 0 , 2, 0, "kingPCs_EMMAX_min13")
# 
# tab_noC3_CEP_HAUS_mpnst = do.call("rbind.data.frame", cpx_res_list_noC3_mpnst)
# write.table(tab_noC3_CEP_HAUS_mpnst, paste(path, "tab_noC3_CEP_HAUS_mpnst.tsv", sep = ""), sep = "\t",
#             row.names = F, quote = F)


##Shelterin no subtraction
##see pca_pop_plots_country_ISKS.R; lines 104 - 107
 # cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 0, 0, 0, 0, "ISKSvsMGRB")
 # cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 2 , 0 , 13, 0, "ISKSvsMGRB")
 # cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 0 , 0 , 0, 0, "ISKSvsMGRB")
 # cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core_mod, 0 , 0 , 0, 0, "ISKSvsMGRB")
 # cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 0 , 0 , 13, 0, "ISKSvsMGRB")
 # cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core_mod, 0 , 0 , 13, 0, "ISKSvsMGRB")

##EMMAX output comparison
# std_PC1234 = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_skat_para_result_isks_combset2020_uni_MAF_PC1234_ver4_clinrect_Aug31.rds")
# std_PC1234_min13 = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/Exome_skat_para_result_isks_totex_count_PC1234_ver4_clinrect_Aug31_subPC1.rds")
# std_PC1234_EMMAX_min13 = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/Exome_skat_para_result_isks_totex_count_PC1234_subPC1_EMMAX.rds")
# std_PCair1234_EMMAX_min13 = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/Exome_skat_para_result_isks_totex_count_pcair1234_subPC1_kinEMMAX.rds")
# #std_grm_EMMAX_min13 = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/Exome_skat_para_result_isks_totex_count_noPC1234_ver4_clinrect_Aug31_subPC1_grmEMMAX.rds")
# std_grm_GCTA_min13 = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/Exome_skat_para_result_isks_totex_count_noPC1234_ver4_clinrect_Aug31_subPC1_grmGCTA_cpxgenes.rds")


std_PC1234 = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_skat_para_result_isks_combset2020_uni_MAF_PC1234_ver4_clinrect_Aug31.rds")
std_PCs = std_PC1234[[4]]
#std_PC1234_min13 = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/Exome_skat_para_result_isks_totex_count_PC1234_ver4_clinrect_Aug31_subPC1.rds")
#std_PCs_min13 = std_PC1234_min13[[4]]
std_PC1234_min13_cc = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_2/SKAT/Exome_skat_para_result_isks_totex_count_PC1234_ver4_clinrect_Jan10_ccrect_PCnew.rds")
std_PC1234_EMMAX_min13_cc = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_2/SKAT/Exome_skat_para_result_isks_totex_count_PC1234_subPC1_EMMAX_Jan10_ccrect_PCnew.rds")
#std_PCair1234_EMMAX_min13_cc = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_2/SKAT/Exome_skat_para_result_isks_totex_count_pcair1234_kinEMMAX_Jan10_ccrect.rds")
std_grm_GCTA_min13_cc = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_2/SKAT/Exome_skat_para_result_isks_totex_count_noPC1234_ver4_clinrect_grmGCTA_Jan10_ccrect_PCnew.rds")
std_PCs_EUR = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_2/SKAT/Exome_skat_para_result_isks_totex_count_PC1234_ver4_clinrect_Jan10_EUR.rds")

#Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")
##for paper : Nov-19-2021
Shelterin <- c("POT1", "TINF2", "TERF1", "SMARCAL1", "STAG3")

CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
#CEP_HAUS_core_mod <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1")

#MPNST_pos <- c("NF1", "LZTR1", "SDHA", "SDHB", "SDHC", "SDHD")
##for paper : Nov-19-2021
MPNST_pos <- c("NF1")

all_genes = unique(c(Shelterin, CEP_HAUS_core, MPNST_pos, "TP53", "BRCA2", "BRCA1", "PALB2"))


get_gene_pval <- function(skatlist_inp)
{
  skat_df <- do.call("rbind.data.frame", skatlist_inp)
  skat_df_sub =  skat_df[skat_df[,1] %in% all_genes,]
  obj <- deparse(substitute(skatlist_inp))
  skat_df_sub$Class = obj
  cnames = colnames(skat_df_sub)
  cnames = gsub("pval_SKATbin", "pval_SKATburden", cnames)
  colnames(skat_df_sub) = cnames
  return(skat_df_sub)
}



std_PC1234_cpx = get_gene_pval(std_PCs)
std_PC1234_cpx = std_PC1234_cpx[,c(1,3,6)]
colnames(std_PC1234_cpx)[2] = "pval_SKATburden"
#std_PCs_min13 = std_PC1234_min13[[4]]
#std_PC1234_min13_cpx = get_gene_pval(std_PCs_min13)
#std_PC1234_min13_cpx = std_PC1234_min13_cpx[,c(1,3,6)]
#colnames(std_PC1234_min13_cpx)[2] = "pval_SKATburden"
std_PC1234_min13_cc_cpx = get_gene_pval(std_PC1234_min13_cc)
std_PC1234_EMMAX_min13_cc_cpx = get_gene_pval(std_PC1234_EMMAX_min13_cc)
#std_PCair1234_EMMAX_min13_cc_cpx = get_gene_pval(std_PCair1234_EMMAX_min13_cc)
#std_grm_EMMAX_min13_cpx = get_gene_pval(std_grm_EMMAX_min13)
std_grm_GCTA_min13_cc_cpx = get_gene_pval(std_grm_GCTA_min13_cc)
#std_grm_GCTA_min13_cpx = get_gene_pval(std_grm_GCTA_min13)
std_PC1234_EUR_cpx = get_gene_pval(std_PCs_EUR)

##signal for centrosome complex was lost upon exclusion of non-EUR samples in std_PC1234_EUR_cpx.
################
##make separate combined_df_pval dataframe for different cohort sizes

#comb_df_pval = rbind.data.frame(std_PC1234_cpx, std_PC1234_min13_cpx, std_PC1234_EMMAX_min13_cpx, std_PCair1234_EMMAX_min13_cpx, std_grm_EMMAX_min13_cpx)
##for fig S7
comb_df_pval = rbind.data.frame(std_PC1234_cpx, std_PC1234_min13_cc_cpx, std_PC1234_EMMAX_min13_cc_cpx,
                                std_grm_GCTA_min13_cc_cpx)
##for table S13
#comb_df_pval = rbind.data.frame(std_PC1234_min13_cc_cpx, std_PC1234_EMMAX_min13_cc_cpx,
#                                std_PCair1234_EMMAX_min13_cc_cpx, std_grm_GCTA_min13_cc_cpx)

#comb_df_pval = rbind.data.frame(std_PC1234_cpx, std_PC1234_min13_cpx, std_PC1234_EMMAX_min13_cpx, std_grm_GCTA_min13_cpx)

comb_df_pval$log10pval = -log10(comb_df_pval$pval_SKATburden)

library(ggplot2)
library(cowplot)
samp_ord = c("TP53", MPNST_pos, Shelterin, CEP_HAUS_core,  "BRCA1", "BRCA2", "PALB2")
comb_df_pval$eg_ID <- factor(comb_df_pval$eg_ID,levels=samp_ord,ordered=TRUE)

# ggplot(comb_df_pval, aes(x = eg_ID, y = log10pval, color = Class)) + 
#   geom_point() + 
#   geom_hline(yintercept = 1, lty=2, colour = "black") + theme_cowplot(font_size = 12) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   theme(legend.position = "bottom")
# 
# ggplot(comb_df_pval, aes(x = eg_ID, y = log10pval)) + 
#   geom_boxplot() + 
#   geom_hline(yintercept = 1, lty=2, colour = "black") + theme_cowplot(font_size = 12) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   theme(legend.position = "bottom")
# 
`%nin%` = Negate(`%in%`)
#comb_df_pval = comb_df_pval[comb_df_pval$Class %nin% c("std_PCs", "std_PCair1234_EMMAX_min13"),]
#comb_df_pval = comb_df_pval[comb_df_pval$Class %nin% c("std_PCair1234_EMMAX_min13_cc"),]
comb_df_pval$Class = gsub("std_PCs", "PCs_base",comb_df_pval$Class)
comb_df_pval$Class = gsub("std_PC1234_min13_cc", "PCs_min61",comb_df_pval$Class)
comb_df_pval$Class = gsub("std_PC1234_EMMAX_min13_cc", "PCs_nullEMMAX_min61",comb_df_pval$Class)
comb_df_pval$Class = gsub("std_grm_GCTA_min13_cc", "stdGRM_nullEMMAX_min61",comb_df_pval$Class)
# comb_df_pval$complex = ifelse(comb_df_pval$eg_ID %in% Shelterin, "Shelterin", 
#                               ifelse(comb_df_pval$eg_ID %in% CEP_HAUS_core, "Centrosome", 
#                                      ifelse(comb_df_pval$eg_ID %in% MPNST_pos , "MPNST",
#                                             ifelse(comb_df_pval$eg_ID %in% c("BRCA1", "BRCA2", "PALB2"), "HBOC", "TP53"))))

##for paper Nov-19-2021
comb_df_pval$complex = ifelse(comb_df_pval$eg_ID %in% Shelterin, "Shelterin", 
                              ifelse(comb_df_pval$eg_ID %in% CEP_HAUS_core, "Centrosome", 
                                     ifelse(comb_df_pval$eg_ID %in% MPNST_pos , "NF1",
                                            ifelse(comb_df_pval$eg_ID %in% c("BRCA1", "BRCA2", "PALB2"), "HBOC", "TP53"))))

###Final plot
ggplot(comb_df_pval, aes(x=eg_ID, y=log10pval)) +  
  geom_boxplot(show.legend = F) +
  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.3), 
             aes(color=factor(Class)), show.legend = T) +
  geom_hline(yintercept = 1, lty=2, colour = "black") +
  theme_cowplot(font_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  xlab("") + ylab("-log10(P.value)") 
#theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))
#+ guides(fill = guide_legend(reverse = TRUE))
p1 = ggplot(comb_df_pval, aes(x=eg_ID, y=log10pval)) +  
  geom_boxplot(show.legend = F) +
  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.3), 
             aes(color=factor(Class)), show.legend = T) +
  geom_hline(yintercept = 1, lty=2, colour = "black") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "top", legend.title = element_blank()) +
  xlab("") + ylab("-log10(P.value)") 

#p1 + theme(legend.text.align = 3)


#####Table S13 redo

##make table for hypergeometric test
tab_fisher = function(comb_model_cpx_genes, count_tab){
  classes = unique(comb_model_cpx_genes$Class)
  if(length(classes) > 1) {
  sig_genes_all = list()
  for(i in 1:length(classes)){
  sig_genes = comb_model_cpx_genes[comb_model_cpx_genes$Class %in% classes[i] & 
                                     comb_model_cpx_genes$pval_SKATburden < 0.1, ]
  sig_genes[,6:8] = count_tab[match(sig_genes[,1], count_tab$gene),1:3 ]
  colnames(sig_genes[6:8]) = c("case", "control", "gene")
  sig_genes_all[[i]] = sig_genes
  }
  sig_genes_all_df = do.call("rbind.data.frame", sig_genes_all)
  return(sig_genes_all_df)
  }
  else{
    sig_genes = comb_model_cpx_genes[comb_model_cpx_genes$Class %in% classes & 
                                       comb_model_cpx_genes$pval_SKATburden < 0.1, ]
    print(dim(sig_genes))
    sig_genes[,4:6] = count_tab[match(sig_genes[,1], count_tab$gene),1:3 ]
    colnames(sig_genes[4:6]) = c("case", "control", "gene")
    return(sig_genes)
  }
  #  sig_genes = sig_genes[!is.na(sig_genes$ISKS),]
}

##make separate combined_df_pval dataframe for different cohort sizes
comb_df_pval_counts_noC3_cc = tab_fisher(comb_df_pval, isks_mgrb_genes_noC3_cc) ##with case control optimization
##minus13 (sub_PC1)
#comb_df_pval_m13 = rbind.data.frame(std_PC1234_min13_cpx, std_grm_GCTA_min13_cpx)
PC1234_min13_pval_counts_noC3_m13 = tab_fisher(std_PC1234_min13_cpx, isks_mgrb_genes_noC3_m13) ##with case control optimization
GRM_min13_pval_counts_noC3_m13 = tab_fisher(std_grm_GCTA_min13_cpx, isks_mgrb_genes_noC3_m13) ##with case control optimization
comb_df_pval_counts_noC3_m13 = rbind.data.frame(PC1234_min13_pval_counts_noC3_m13, GRM_min13_pval_counts_noC3_m13)
comb_df_pval_counts_noC3_m13$complex = comb_df_pval_counts_noC3_cc[match(comb_df_pval_counts_noC3_m13$eg_ID, comb_df_pval_counts_noC3_cc$eg_ID),5]


##original data
comb_df_pval_counts_noC3 = tab_fisher(std_PC1234_cpx, isks_mgrb_genes_noC3) ##with case control optimization

comb_df_pval_counts_noC3$complex = comb_df_pval_counts_noC3_cc[match(comb_df_pval_counts_noC3$eg_ID, comb_df_pval_counts_noC3_cc$eg_ID),5]


##Run Fisher test for all models

cpx_OR_fisher_count <- function(comb_df_inp, mod_name, case_coh_size, cont_coh_size, cpx_name, sub_case, sub_cont, sub_ISKS, sub_MGRB){
  print(cpx_name)
  print(mod_name)
  comb_df_inp_tab <- comb_df_inp[comb_df_inp$complex %in% cpx_name & comb_df_inp$Class %in% mod_name,]
  print((comb_df_inp_tab))
  # ppi_res_tab[,2] <- ifelse(ppi_res_tab[,2] == 0, 1, ppi_res_tab[,2])
  # inp <- c(sum(ppi_res_tab[,1]) - sub_case , case_coh_size - sum(ppi_res_tab[,1]) - sub_ISKS , 
  #          sum(ppi_res_tab[,2]) - sub_cont, cont_coh_size - sum(ppi_res_tab[,2]) - sub_MGRB)
  inp <- c(sum(comb_df_inp_tab[,6], na.rm = T) - sub_case , case_coh_size - sub_ISKS - (sum(comb_df_inp_tab[,6], na.rm = T) - sub_case), 
           sum(comb_df_inp_tab[,7],na.rm = T) - sub_cont, cont_coh_size - sub_MGRB - (sum(comb_df_inp_tab[,7], na.rm = T) - sub_cont))
  sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
  colnames(sim_mat) <- c("case", "cont")
  rownames(sim_mat) <- c("hits", "no_hits")
  print(sim_mat)
  #ft <- fisher.test(sim_mat, alternative = "greater")
  #ft <- fisher.test(sim_mat)
  cpx_lab = deparse(substitute(cpx_name))
  ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.95)
  ft_df <- cbind.data.frame("gene" = cpx_lab ,"Cases" = sim_mat[1,1],
                            "Controls" = sim_mat[1,2],
                            "Fish_pval" = ft$p.value,"CI_lower" = ft$conf.int[1],
                            "CI_upper" = ft$conf.int[2],
                            "OR_Fish" = ft$estimate, "case_coh_size" = sum(sim_mat[,1]),
                            "cont_coh_size" = sum(sim_mat[,2]), "model" = mod_name)
  return(ft_df)
}

##reports only genes that are significant in SKATburden test and are present in the complex
##CC
cpx_OR_fisher_count(comb_df_pval_counts_noC3_cc, "std_PC1234_min13_cc", 1585, 3203, "Shelterin", 0, 0, 0, 0)
cpx_OR_fisher_count(comb_df_pval_counts_noC3_cc, "std_PC1234_min13_cc", 1585, 3203, "Centrosome", 0, 0, 0, 0)
cpx_OR_fisher_count(comb_df_pval_counts_noC3_cc, "std_PC1234_EMMAX_min13_cc", 1585, 3203, "Shelterin", 0, 0, 0, 0)
cpx_OR_fisher_count(comb_df_pval_counts_noC3_cc, "std_PC1234_EMMAX_min13_cc", 1585, 3203, "Centrosome", 0, 0, 0, 0)
cpx_OR_fisher_count(comb_df_pval_counts_noC3_cc, "std_grm_GCTA_min13_cc", 1585, 3203, "Shelterin", 0, 0, 0, 0)
cpx_OR_fisher_count(comb_df_pval_counts_noC3_cc, "std_grm_GCTA_min13_cc", 1585, 3203, "Centrosome", 0, 0, 0, 0)

##minus13
comb_df_pval_counts_noC3_m13 = comb_df_pval_counts_noC3_m13[,c(1:3,6:7,4:5)]
cpx_OR_fisher_count(comb_df_pval_counts_noC3_m13, "std_PCs_min13", 1631, 3205, "Shelterin", 0, 0, 0, 0)
cpx_OR_fisher_count(comb_df_pval_counts_noC3_m13, "std_PCs_min13", 1631, 3205, "Centrosome", 0, 0, 0, 0)
cpx_OR_fisher_count(comb_df_pval_counts_noC3_m13, "std_grm_GCTA_min13", 1631, 3205, "Shelterin", 0, 0, 0, 0)
cpx_OR_fisher_count(comb_df_pval_counts_noC3_m13, "std_grm_GCTA_min13", 1631, 3205, "Centrosome", 0, 0, 0, 0)


##base model
comb_df_pval_counts_noC3 = comb_df_pval_counts_noC3[,c(1:3,6:7,4:5)]
cpx_OR_fisher_count(comb_df_pval_counts_noC3, "std_PCs", 1644, 3205, "Shelterin", 0, 0, 0, 0)
cpx_OR_fisher_count(comb_df_pval_counts_noC3, "std_PCs", 1644, 3205, "Centrosome", 0, 0, 0, 0)



########################################################
##SKAT for complexes
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)
library(SKAT)


##use fil_tab and PCs after removal of samples using centroid method
fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze_rev_PC1_sub.tsv",
                      sep = "\t", stringsAsFactors = F, header = T) ##generated by R2Q6 script pca_pop*ISKS.R
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)

##Additive model
Ex_samp_id <- unique(fil_tab$SAMPLE)

######get phenotype data to control for age and sex and PC's; 

QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv", 
                      header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %in% "TRUE",]
QC2_dat_pass$isFemale <- ifelse(QC2_dat_pass$f_stat < 0.2, 1, 
                                ifelse(QC2_dat_pass$f_stat > 0.8, 0, 2))
##for future QC with prob.NFE correction use the following file and adjust the column for gender from 39 to 53

#p_Data_noCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_2/rect_case_cont_ratio_PCA_clustered_jan2022.rds")
##p_Data_noCH below is generated by pca_tsne_diagnostics.Rmd script: Added country and cluster
p_Data_noCH <- readRDS("~/RVAS/comb_set_2020/pop_PCA/review_2/MGRB_ISKS_1000G_rect_PCnew.rds") #change column names in Null model back to 54:55
p_Data_noCH$rect_sam = p_Data_noCH$sample # the samples name are rectified at HAIL stage
p_Data_noCH <- p_Data_noCH[as.character(p_Data_noCH$rect_sam) %in% Ex_samp_id,]
##drop MGRB sample BAAUD; not in latest call
#samp_ID_match <- samp_ID_match[grep("BAAUD", samp_ID_match, invert = T)]
#p_Data_noCH <- p_Data_noCH[match(samp_ID_match, p_Data_noCH$sample),]
Ex_samp_id <- Ex_samp_id[match(p_Data_noCH$rect_sam, Ex_samp_id)]
##filter out QC fail cases
fil_tab <- fil_tab[fil_tab$SAMPLE %in% Ex_samp_id,]

p_Data_noCH$gender <- QC2_dat_pass[match(p_Data_noCH$rect_sam, QC2_dat_pass$new_sampleid), 21]
#write.table(p_Data_noCH, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/pheno_data_SKAT_final_freeze.tsv",
#            sep = "\t", row.names = F, quote = F)

##Add exome count
##ref: Exome sequencing in amyotrophic lateral sclerosis implicates a novel gene, DNAJC7, encoding a heat-shock protein
##total exome count (summation of synonymous variants, benign missense variants, damaging missense variants and PTVs) 
tot_exome_count <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/all_isks_mgrb_total_exome_count_clingene_rnd3_freeze.tsv",
                              sep = "\t", header = T, stringsAsFactors = F)
p_Data_noCH$totex_count <- tot_exome_count[match(p_Data_noCH$rect_sam, tot_exome_count$Var1), 2]

#binary phenotype vector 
p_vec <- ifelse(!is.na(as.numeric(as.character(p_Data_noCH$rect_sam))) | grepl("^CR|^LK",as.character(p_Data_noCH$rect_sam)), 1, 0)


##########################################
##for complexes after removal of POPMAX variant > 0.001, with new PCs
##parallel process complexes
para_SKAT_cpx_var <- function(genes, cpx_name, rm_var_file){
  skat_pvals_sex_pc12 = list()
  
  ftemp_tab <- fil_tab[fil_tab$gene_symbol %in% genes,]
  ftemp_tab <- ftemp_tab[ftemp_tab$comb_score >= 5.6 & as.numeric(ftemp_tab$VAF) >= 0.35, ]
  ftemp_tab <- ftemp_tab[ftemp_tab$VARIANT %nin% rm_var_file,]
  
  print(max(ftemp_tab$comb_score))
  if(dim(ftemp_tab) == 0 ){
    next
  }
  else{
    ftemp_tab_var_id <- unique(ftemp_tab$VARIANT)
    samp_vec <- list()
    for(m in 1:length(ftemp_tab_var_id)){
      sam_gene_gt <- ftemp_tab[ftemp_tab$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,127:128)]
      sam_gene_gt <- unique(sam_gene_gt)
      if(dim(sam_gene_gt)[1] > 1 & length(unique(sam_gene_gt$vep_consequence)) > 1){
        # if(dim(sam_gene_gt)[1] > 1 & length(unique(sam_gene_gt$vep_consequence)) > 1 & length(unique(sam_gene_gt$comb_score)) > 1){
        vep_con <- unique(sam_gene_gt$vep_consequence)
        samp_vec_con <- list()
        for(k in 1:length(vep_con)){
          sam_gene_gt_con <- sam_gene_gt[sam_gene_gt$vep_consequence %in% vep_con[k],]
          sam_gene_gt_con <- unique(sam_gene_gt_con)
          maf_vec_cont <- sum(dim(sam_gene_gt_con[grepl("^[ABZ]",sam_gene_gt_con$SAMPLE),])[1])/(2*length(p_vec))
          maf_vec_case <- sum(dim(sam_gene_gt_con[!grepl("^[ABZ]",sam_gene_gt_con$SAMPLE),])[1])/(2*length(p_vec))
          ##genotype matrix  
          sam_gene_gt_con$add_mod <- as.numeric(sam_gene_gt_con$GT)
          sam10 <- ifelse(Ex_samp_id %in% sam_gene_gt_con$SAMPLE, 1, 0)
          sam10[which(sam10 != 0)] <- sam_gene_gt_con$add_mod ##additive model
          samp_vec_con[[k]] <- c(unique(sam_gene_gt_con$VARIANT),
                                 unique(sam_gene_gt_con$gene_symbol), 
                                 unique(sam_gene_gt_con$vep_consequence), 
                                 unique(as.character(sam_gene_gt_con$auto_call)),
                                 as.numeric(unique(sam_gene_gt_con$comb_score)), as.numeric(maf_vec_case),
                                 as.numeric(maf_vec_cont), sam10)
        }
        #  vep_cod <- sapply(str_extract_all(sam_gene_gt$vep_Codons, "[A-Z]+"), paste, collapse= '_')
        #  ref_alt <- unlist(lapply(strsplit(sam_gene_gt$VARIANT, split = ":"), function(x) paste(x[3], x[4], sep = "_")))
        #  sam_gene_gt <- sam_gene_gt[which(vep_cod %in% ref_alt),]
        
      }
      else{
        ##compute cohort specific MAF
        maf_vec_cont <- sum(ifelse(is.na(as.numeric(sam_gene_gt$SAMPLE)), 1, 0))/(2*length(p_vec))
        maf_vec_case <- sum(ifelse(!is.na(as.numeric(sam_gene_gt$SAMPLE)) | 
                                     grepl("^CR|^LK", as.character(sam_gene_gt$SAMPLE)), 1, 0))/(2*length(p_vec))
        #maf_vec <- (maf_vec_cont + maf_vec_case)/(2*(1572 + 1110))
        ##genotype matrix  
        sam_gene_gt$add_mod <- as.numeric(sam_gene_gt$GT)
        sam10 <- ifelse(Ex_samp_id %in% sam_gene_gt$SAMPLE, 1, 0)
        sam10[which(sam10 != 0)] <- sam_gene_gt$add_mod ##additive model
        ##account for discrepancy in vep_consequence
        #vep_con <- names(sort(table(unlist(strsplit(sam_gene_gt$vep_consequence, split = "&"))),decreasing=TRUE)[1])
        # sam_gene_gt$vep_consequence <- vep_con
        
        
        samp_vec[[m]] <- c(ftemp_tab_var_id[m],
                           unique(sam_gene_gt$gene_symbol), 
                           unique(sam_gene_gt$vep_consequence), 
                           unique(as.character(sam_gene_gt$auto_call)),
                           as.numeric(unique(sam_gene_gt$comb_score)), as.numeric(maf_vec_case),
                           as.numeric(maf_vec_cont), sam10)
      }
    }
    samp_vec_mat_uni <- do.call("rbind.data.frame", samp_vec)
    colnames(samp_vec_mat_uni) <- c("VARIANT", "cpx_name", "vep_consequence", "auto_call", "comb_score", "coh_MAF_case",
                                    "coh_MAF_cont", Ex_samp_id)
    if(exists("samp_vec_con")){
      samp_vec_mat_con <- do.call("rbind.data.frame", samp_vec_con)
      colnames(samp_vec_mat_con) <- c("VARIANT", "cpx_name", "vep_consequence", "auto_call", "comb_score", "coh_MAF_case",
                                      "coh_MAF_cont", Ex_samp_id)
      samp_vec_mat <- rbind.data.frame(samp_vec_mat_uni, samp_vec_mat_con)
    }else{
      samp_vec_mat <- samp_vec_mat_uni
    }
    print(dim(samp_vec_mat))
    ##Intracohort filter  : equivalent to ~ 5/1661*2 for case ; ~ 5/3209*2 for control
    ##changed to uniform intra_cohort_MAF filter
    ## change to entire cohort MAF filter to 3/(4849*2)##latest Aug.31
    #     samp_vec_mat <- samp_vec_mat[as.numeric(as.character(samp_vec_mat$coh_MAF_case)) <= 0.0015 | 
    #                                    as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) <= 0.001,]
    samp_vec_mat <- samp_vec_mat[!(as.numeric(as.character(samp_vec_mat$coh_MAF_case)) >= 0.00035 | 
                                     as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) >= 0.00035),]
    #  samp_vec_mat <- samp_vec_mat[as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) <= 0.001,]
    
    if(is.null(dim(samp_vec_mat)) | dim(samp_vec_mat)[1] == 0) {
      
      skat_pvals_sex_pc12 <- NULL
      
      next
    }
    ##compute maf SKAT takes care of the AFs internally by assigning higher weights to rare variants 
    
    else {
      gene_mat_comb <- as.numeric(as.character(samp_vec_mat[,5]))
      gene_mat <- as.matrix(samp_vec_mat[,-c(1:7)])
      class(gene_mat) <- "numeric"
      
      ##gender + totex + pc1:4
      skat_pvals_sex_pc12 <- SKAT_run(geno_mat = gene_mat, gene = cpx_name, x=p_Data_noCH[,c(55:56,19:22)], p_vec, cust_weight = gene_mat_comb, rho = 1)
      
    }
    
  } ##end of nested else loop
  return(skat_pvals_sex_pc12)
}##end of gene loop ; i loop   



Shelterin <- c("POT1", "TINF2", "TERF1", "SMARCAL1", "STAG3")
Sarcoma_genes <- c("TP53", "NF1", "SDHA", "SDHB", "SDHD")
CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")

Breast_cancer_genes = c("BRCA2", "BRCA1", "PALB2")

MPNST_pos <- c("NF1", "LZTR1", "SDHA", "SDHB", "SDHC", "SDHD")

cpx_list <- list(Shelterin, CEP_HAUS_core, Breast_cancer_genes, Sarcoma_genes, MPNST_pos)
names(cpx_list) <- c("Shelterin", "Centrosome", "BRCA_genes", "SARC_genes", "MPNST")

##Parallelised SKAT  
library(doParallel)
library(doMC)
registerDoMC(30)

##POPMAX filtered variants
dbnsfp_cpx = read.delim("~/RVAS/dbNSFP_4.1a/isks_cpx_var_dbnsfp_out_fin.tsv", sep = "\t", header = T, stringsAsFactors = F)
for(i in 6:9){
  dbnsfp_cpx[,i] <- as.numeric(dbnsfp_cpx[,i])
}
dbnsfp_cpx[is.na(dbnsfp_cpx)] = 0
dbnsfp_cpx =  unique(dbnsfp_cpx)
rownames(dbnsfp_cpx) = gsub(" ", "", apply(dbnsfp_cpx[ ,1:4] , 1 , paste, collapse = ":" ))

dbnsfp_cpx_ind_01 = apply(dbnsfp_cpx[,6:9], 1, function(x)ifelse(max(x) >= 0.001, 1, 0))
dbnsfp_cpx_samp_01 = names(dbnsfp_cpx_ind_01[dbnsfp_cpx_ind_01 == 1])


res20_cpx_var_pc <- list()
#genes <- unique(fil_tab$gene_symbol)
#genes <- names(cpx_list)
system.time(res20_cpx_var_pc <- foreach(i=1:length(cpx_list), .errorhandling = 'remove') %dopar% 
{para_SKAT_cpx_var(cpx_list[[i]], names(cpx_list)[i], dbnsfp_cpx_samp_01)}) ##from POPMAX : line 454

res20_cpx_var_pc_df = do.call("rbind.data.frame", res20_cpx_var_pc)
res20_cpx_var_pc_df$p.adj = p.adjust(res20_cpx_var_pc_df$pval_SKATbin)
