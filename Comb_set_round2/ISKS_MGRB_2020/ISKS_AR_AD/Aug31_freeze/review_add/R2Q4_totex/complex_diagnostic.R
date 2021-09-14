##R2Q4; remove SSNA1 from CEP_HAUS_core
Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")

CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
CEP_HAUS_core_mod <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1")

MPNST_pos <- c("NF1", "LZTR1", "SDHA", "SDHB", "SDHC", "SDHD")

##ISKS vs MGRB
isks_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_para_tab_ISKS_MGRB_C345_Aug31.rds")
isks_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Exome_para_tab_ISKS_MGRB_minusC3_Aug31.rds")

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
cpx_res_list[[1]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 0 , 0 , 0, 0, "original")
##totex + rm_13AFR (SSNA1 p-val < 0.1, 3 ISKS lost)
cpx_res_list[[2]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 3 , 0 , 13, 0, "PCs_min13")
##totex + rm_13AFR + ancestry PCs + ancestry adjusted kinship matrix (CEP72 p-val < 0.1, 7 ISKS + 5 MGRB)
cpx_res_list[[3]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 7 , 5 , 13, 0, "PCs_EMMAX_min13")
##totex + rm_13AFR + GRM (HAUS4 & SSNA1 p-val < 0.1, 6 ISKS)
cpx_res_list[[4]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 6 , 0 , 13, 0, "kingGRM_EMMAX_min13")
##totex + rm_13AFR + kinship  PCs + ancestry adjusted kinship matrix (HAUS4 & SSNA1 p-val < 0.1, 6 ISKS)
cpx_res_list[[5]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 6 , 0 , 13, 0, "kingPCs_EMMAX_min13")

tab_C345_CEP_HAUS = do.call("rbind.data.frame", cpx_res_list)
path = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/review_add/R2Q6_popstruc/Results/"
write.table(tab_C345_CEP_HAUS, paste(path, "tab_fisher_C345_CEP_HAUS.tsv", sep = ""), sep = "\t",
            row.names = F, quote = F)
##noC3 variants
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
cpx_res_list_Shel = list()
##with C345 variants
cpx_res_list_Shel[[1]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 0 , 0 , 0, 0, "original")
##totex + rm_13AFR (SSNA1 p-val < 0.1, 3 ISKS lost)
cpx_res_list_Shel[[2]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 0 , 0 , 13, 0, "PCs_min13")
##totex + rm_13AFR + ancestry PCs + ancestry adjusted kinship matrix (CEP72 p-val < 0.1, 7 ISKS + 5 MGRB)
cpx_res_list_Shel[[3]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 0 , 5 , 13, 0, "PCs_EMMAX_min13")
##totex + rm_13AFR + GRM (HAUS4 & SSNA1 p-val < 0.1, 6 ISKS)
cpx_res_list_Shel[[4]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 0 , 0 , 13, 0, "kingGRM_EMMAX_min13")
##totex + rm_13AFR + kinship  PCs + ancestry adjusted kinship matrix (HAUS4 & SSNA1 p-val < 0.1, 6 ISKS)
cpx_res_list_Shel[[5]] = cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 0 , 0 , 13, 0, "kingPCs_EMMAX_min13")

tab_C345_SHEL = do.call("rbind.data.frame", cpx_res_list_Shel)
write.table(tab_C345_SHEL, paste(path, "tab_C345_SHEL.tsv", sep = ""), sep = "\t",
            row.names = F, quote = F)

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
cpx_res_list_noC3_mpnst = list()
cpx_res_list_noC3_mpnst[[1]] = cpx_OR_fisher_one(isks_mpnst_mgrb_genes_noC3, 157, 3205, CEP_HAUS_core, 0 , 0 , 0, 0, "original")
##totex + rm_13AFR (SSNA1 p-val < 0.1, 3 ISKS lost)
cpx_res_list_noC3_mpnst[[2]] = cpx_OR_fisher_one(isks_mpnst_mgrb_genes_noC3, 157, 3205, CEP_HAUS_core, 0 , 0 , 2, 0, "PCs_min13")
##totex + rm_13AFR + ancestry PCs + ancestry adjusted kinship matrix (CEP72 p-val < 0.1, 6 ISKS + 2 MGRB)
cpx_res_list_noC3_mpnst[[3]] = cpx_OR_fisher_one(isks_mpnst_mgrb_genes_noC3, 157, 3205, CEP_HAUS_core, 2 , 2 , 2, 0, "PCs_EMMAX_min13")
##totex + rm_13AFR + GRM (HAUS4 p-val < 0.1, 2 ISKS)
cpx_res_list_noC3_mpnst[[4]] = cpx_OR_fisher_one(isks_mpnst_mgrb_genes_noC3, 157, 3205, CEP_HAUS_core, 1 , 0 , 2, 0, "kingGRM_EMMAX_min13")
##totex + rm_13AFR + kinship  PCs + ancestry adjusted kinship matrix (HAUS4 < 0.1, 2 ISKS)
cpx_res_list_noC3_mpnst[[5]] = cpx_OR_fisher_one(isks_mpnst_mgrb_genes_noC3, 157, 3205, CEP_HAUS_core, 1 , 0 , 2, 0, "kingPCs_EMMAX_min13")

tab_noC3_CEP_HAUS_mpnst = do.call("rbind.data.frame", cpx_res_list_noC3_mpnst)
write.table(tab_noC3_CEP_HAUS_mpnst, paste(path, "tab_noC3_CEP_HAUS_mpnst.tsv", sep = ""), sep = "\t",
            row.names = F, quote = F)


##Shelterin no subtraction
##see pca_pop_plots_country_ISKS.R; lines 104 - 107
 # cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 0, 0, 0, 0, "ISKSvsMGRB")
 # cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 2 , 0 , 13, 0, "ISKSvsMGRB")
 # cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 0 , 0 , 0, 0, "ISKSvsMGRB")
 # cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core_mod, 0 , 0 , 0, 0, "ISKSvsMGRB")
 # cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 0 , 0 , 13, 0, "ISKSvsMGRB")
 # cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core_mod, 0 , 0 , 13, 0, "ISKSvsMGRB")

##EMMAX output comparison
std_PC1234 = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_skat_para_result_isks_combset2020_uni_MAF_PC1234_ver4_clinrect_Aug31.rds")
std_PC1234_min13 = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/Exome_skat_para_result_isks_totex_count_PC1234_ver4_clinrect_Aug31_subPC1.rds")
std_PC1234_EMMAX_min13 = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/Exome_skat_para_result_isks_totex_count_PC1234_subPC1_EMMAX.rds")
std_PCair1234_EMMAX_min13 = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/Exome_skat_para_result_isks_totex_count_pcair1234_subPC1_kinEMMAX.rds")
std_grm_EMMAX_min13 = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/Exome_skat_para_result_isks_totex_count_noPC1234_ver4_clinrect_Aug31_subPC1_grmEMMAX.rds")

Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")

CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
CEP_HAUS_core_mod <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1")

MPNST_pos <- c("NF1", "LZTR1", "SDHA", "SDHB", "SDHC", "SDHD")

all_genes = unique(c(Shelterin, CEP_HAUS_core, MPNST_pos, "TP53", "BRCA2", "BRCA1", "PALB2"))


get_gene_pval <- function(skatlist_inp)
{
  skat_df <- do.call("rbind.data.frame", skatlist_inp)
  skat_df_sub =  skat_df[skat_df[,1] %in% all_genes,]
  obj <- deparse(substitute(skatlist_inp))
  skat_df_sub$Class = obj
  return(skat_df_sub)
}


std_PCs = std_PC1234[[4]]
std_PC1234_cpx = get_gene_pval(std_PCs)
std_PC1234_cpx = std_PC1234_cpx[,c(1,3,6)]
colnames(std_PC1234_cpx)[2] = "pval_SKATburden"
std_PCs_min13 = std_PC1234_min13[[4]]
std_PC1234_min13_cpx = get_gene_pval(std_PCs_min13)
std_PC1234_min13_cpx = std_PC1234_min13_cpx[,c(1,3,6)]
colnames(std_PC1234_min13_cpx)[2] = "pval_SKATburden"
std_PC1234_EMMAX_min13_cpx = get_gene_pval(std_PC1234_EMMAX_min13)
std_PCair1234_EMMAX_min13_cpx = get_gene_pval(std_PCair1234_EMMAX_min13)
std_grm_EMMAX_min13_cpx = get_gene_pval(std_grm_EMMAX_min13)

comb_df_pval = rbind.data.frame(std_PC1234_cpx, std_PC1234_min13_cpx, std_PC1234_EMMAX_min13_cpx, std_PCair1234_EMMAX_min13_cpx, std_grm_EMMAX_min13_cpx)
comb_df_pval$log10pval = -log10(comb_df_pval$pval_SKATburden)

library(ggplot2)
library(cowplot)
samp_ord = c("TP53", Shelterin, CEP_HAUS_core, MPNST_pos, "BRCA1", "BRCA2", "PALB2")
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
comb_df_pval = comb_df_pval[comb_df_pval$Class %nin% "std_PCs",]
comb_df_pval$Class = gsub("std_PCs_min13", "PCs_min13",comb_df_pval$Class)
comb_df_pval$Class = gsub("std_PC1234_EMMAX_min13", "PCs_EMMAX_min13",comb_df_pval$Class)
comb_df_pval$Class = gsub("std_PCair1234_EMMAX_min13", "kingPCs_EMMAX_min13",comb_df_pval$Class)
comb_df_pval$Class = gsub("std_grm_EMMAX_min13", "kingGRM_EMMAX_min13",comb_df_pval$Class)
comb_df_pval$complex = ifelse(comb_df_pval$eg_ID %in% Shelterin, "Shelterin", 
                              ifelse(comb_df_pval$eg_ID %in% CEP_HAUS_core, "Centrosome", 
                                     ifelse(comb_df_pval$eg_ID %in% MPNST_pos , "MPNST",
                                            ifelse(comb_df_pval$eg_ID %in% c("BRCA1", "BRCA2", "PALB2"), "HBOC", "TP53"))))

###Final plot
ggplot(comb_df_pval, aes(x=eg_ID, y=log10pval)) +  
  geom_boxplot(show.legend = F) +
  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.3), 
             aes(color=factor(Class)), show.legend = T) +
  geom_hline(yintercept = 1, lty=2, colour = "black") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  xlab("") + ylab("-log10(P.value)") 
#theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))

ggplot(comb_df_pval, aes(x=eg_ID, y=log10pval)) +  
  geom_boxplot(show.legend = F) +
  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.3), 
             aes(color=factor(Class)), show.legend = T) +
  geom_hline(yintercept = 1, lty=2, colour = "black") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = as.factor(complex))) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  xlab("") + ylab("-log10(P.value)") 