##Extended Sheltrin_complex 
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Exome_skat_para_result_isks_combset2020_uni_MAF_PC1234_ver4_rect_test.rds")
df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})
Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Exome_para_pc123_SKAT_Enriched_ISKS_2020.rds")
##Add p-value
Exome_pc123_srt_SKAT_case_enr_nCH_pval <- Map(cbind.data.frame, Exome_pc123_srt_SKAT_case_enr_nCH, Exome_pc123_srt_SKAT)
Sheltrin_comp_extn = c("ACD", "POT1", "TERF1", "TERF2", "TERF2IP", "TINF2", "ATM",
"BAG3", "BLM", "BRCA1", "CALD1", "CLK3", "DCLRE1B", "FANCD2",
"FBXO4", "HSPA4", "KIAA1191", "MRE11A", "NBN", "PINX1", "PRKDC",
"RAD50", "SLX4", "STUB1", "TNKS", "TNKS2", "U2AF2", "UCHL1",
"WRN", "XRCC5", "XRCC6")
sheltrin_extn_skatbin <- Exome_pc123_srt_SKAT_case_enr_nCH_pval[[4]][Exome_pc123_srt_SKAT_case_enr_nCH_pval[[4]]$gene %in% Sheltrin_comp_extn,]

para_fisher_weighted_new <- function(gene_sym){
  
  gene_sym <- gene_sym[3]
  
  ppi_res <- sheltrin_extn_skatbin[sheltrin_extn_skatbin$gene %in% gene_sym,]
  
  ##simple Fisher test
  #inp <- c(deg_union, deg_biog, tot_union - deg_union, tot_biog - deg_biog)
  inp <- c(ppi_res[,1], 1646 - ppi_res[,1] , ppi_res[,2], 3205 - ppi_res[,2])
  sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
  colnames(sim_mat) <- c("isks_enr", "mgrb_enr")
  rownames(sim_mat) <- c("hits", "no_hits")
  #ft <- fisher.test(mgrb_tsg, alternative = "greater")
  ft <- fisher.test(sim_mat, alternative = "greater")
  ##Weighted Fisher
  inp1 <- c(ppi_res[,4], 1646 - ppi_res[,1] , ppi_res[,5], 3205 - ppi_res[,2])
  sim_mat1 <- matrix(inp1 ,nrow = 2, ncol = 2)
  colnames(sim_mat1) <- c("isks_wt", "mgrb_wt")
  rownames(sim_mat1) <- c("hits", "no_hits")
  ft1 <- fisher.test(sim_mat1, alternative = "greater")
  
  ft_df <- cbind.data.frame("gene" =  gene_sym, "Fish_pval" = ft$p.value,
                            "OR_Fish" = ft$estimate, "Wt_Fish_pval" = ft1$p.value,
                            "Wt_Fish_OR" = ft1$estimate)
  
  return(ft_df)
}



library(parallel)
cl <- makeCluster(25)
clusterExport(cl, c("sheltrin_extn_skatbin", "para_fisher_weighted_new"))
system.time(para_fish_wt <- parApply(cl, sheltrin_extn_skatbin, 1, para_fisher_weighted_new))
stopCluster(cl)
wt_fisher_res <- do.call("rbind.data.frame",para_fish_wt)
sheltrin_extn_skatbin <- cbind.data.frame(sheltrin_extn_skatbin, wt_fisher_res)

write.table(sheltrin_extn_skatbin, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/sheltrin_extn_SKATbin_wt_fisher.tsv", sep = "\t", row.names = F, quote = F)
