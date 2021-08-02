
Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")
CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
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
  inp <- c(sum(ppi_res_tab[,1]) - sub_case , case_coh_size - sum(ppi_res_tab[,1]) - sub_ISKS , 
           sum(ppi_res_tab[,2]) - sub_cont, cont_coh_size - sum(ppi_res_tab[,2]) - sub_MGRB)
  sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
  colnames(sim_mat) <- c("case", "cont")
  rownames(sim_mat) <- c("hits", "no_hits")
  #ft <- fisher.test(sim_mat, alternative = "greater")
  #ft <- fisher.test(sim_mat)
  cpx_name = deparse(substitute(cpx_list))
  ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.99)
  ft_df <- cbind.data.frame("gene" = cpx_name ,"Cases" = sim_mat[1,1],
                                 "Controls" = sim_mat[1,2],
                                 "Fish_pval" = ft$p.value,"CI_lower" = ft$conf.int[1],
                                 "CI_upper" = ft$conf.int[2],
                                 "OR_Fish" = ft$estimate, "case_coh_size" = sum(sim_mat[,1]),
                                 "Coh" = coh)
  return(ft_df)
}

##Shelterin no subtraction
##see pca_pop_plots_country_ISKS.R; lines 104 - 107
# cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 0, 0, 0, 0)
# cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, Shelterin, 2 , 0 , 13, 0)
# cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 0 , 0 , 0, 0)
# cpx_OR_fisher_one(isks_mgrb_genes, 1644, 3205, CEP_HAUS_core, 0 , 0 , 13, 0)
