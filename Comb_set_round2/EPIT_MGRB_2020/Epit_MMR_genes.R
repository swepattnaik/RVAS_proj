##MMR genes check: EPIT dataset 
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/round2/Exome_skat_para_result_EPIT_combset2020_uni_MAF_PC1234_ver4_rect_test.rds")

df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})
Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/round2/Exome_para_pc123_SKAT_Enriched_EPIT_2020.rds")
##Add p-value
Exome_pc123_srt_SKAT_case_enr_nCH_pval <- Map(cbind.data.frame, Exome_pc123_srt_SKAT_case_enr_nCH, Exome_pc123_srt_SKAT)

MMR_genes = c("MLH1","MSH2", "MSH6", "PMS2", "POLE", "BRCA1", "BRCA2", "TP53")
MMR_genes_skatbin <- Exome_pc123_srt_SKAT_case_enr_nCH_pval[[4]][Exome_pc123_srt_SKAT_case_enr_nCH_pval[[4]]$gene %in% MMR_genes,]

para_fisher_weighted_new <- function(gene_sym){
  
  gene_sym <- gene_sym[3]
  
  ppi_res <- MMR_genes_skatbin[MMR_genes_skatbin$gene %in% gene_sym,]
  
  ##simple Fisher test
  #inp <- c(deg_union, deg_biog, tot_union - deg_union, tot_biog - deg_biog)
  inp <- c(ppi_res[,1], 920 - ppi_res[,1] , ppi_res[,2], 3209 - ppi_res[,2])
  sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
  colnames(sim_mat) <- c("isks_enr", "mgrb_enr")
  rownames(sim_mat) <- c("hits", "no_hits")
  #ft <- fisher.test(mgrb_tsg, alternative = "greater")
  ft <- fisher.test(sim_mat, alternative = "greater")
  ##Weighted Fisher
  inp1 <- c(ppi_res[,4], 920 - ppi_res[,1] , ppi_res[,5], 3209 - ppi_res[,2])
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
clusterExport(cl, c("MMR_genes_skatbin", "para_fisher_weighted_new"))
system.time(para_fish_wt <- parApply(cl, MMR_genes_skatbin, 1, para_fisher_weighted_new))
stopCluster(cl)
wt_fisher_res <- do.call("rbind.data.frame",para_fish_wt)
MMR_genes_skatbin <- cbind.data.frame(MMR_genes_skatbin, wt_fisher_res)

write.table(MMR_genes_skatbin, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/round2/MMR_genes_SKATbin_wt_fisher.tsv", sep = "\t", row.names = F, quote = F)


##MMR gene variants
##variants with corresponding SKAT p-values
fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/all_epit_mgrb_combset2020_variants_filt_all_fields_risk_factor.tsv",
                      sep = "\t", header = T, stringsAsFactors = F)
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
QC2_dat_pass_mgrb <- QC2_dat_pass[grepl("^[ABZ]", QC2_dat_pass$new_sampleid),]
QC2_dat_pass_mgrb <- QC2_dat_pass_mgrb[,c(1,20:21)]

QC2_bat1 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set/EPIT_MGRB_2020/olga_batch1.final_qc_output.tsv",
                       sep = "\t", header = T, stringsAsFactors = F)
QC2_bat1$isFemale <- ifelse(QC2_bat1$f_stat < 0.2, 1, 
                            ifelse(QC2_bat1$f_stat > 0.8, 0, 2))
QC2_bat1 <- QC2_bat1[,c(1,27:28)]

QC2_bat2 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set/EPIT_MGRB_2020/olga_batch2.final_qc_output.tsv",
                       sep = "\t", header = T, stringsAsFactors = F)
QC2_bat2$isFemale <- ifelse(QC2_bat2$f_stat < 0.2, 1, 
                            ifelse(QC2_bat2$f_stat > 0.8, 0, 2))
QC2_bat2 <- QC2_bat2[,c(1,17:18)]

QC2_epit <- rbind.data.frame(QC2_bat1, QC2_bat2)
QC2_epit_pass <- QC2_epit[QC2_epit$passes_qc2 %in% "TRUE",]

QC2_epit_mgrb_pass <- rbind.data.frame(QC2_epit_pass, QC2_dat_pass_mgrb)
##note p_Data_PC_comb does not change
#p_Data <- read.table("~/RVAS/ISKS_MGRB_gender_age_3665.tsv", header = T, sep = "\t")
p_Data <- read.delim("~/RVAS/Epi_set_2020/pop_PCA/MGRB_EPIT_1000G_combset_pca.scores_clustered.tsv", header = T, sep = "\t", stringsAsFactors = F)
p_Data <- p_Data[p_Data$sample %in% QC2_epit_mgrb_pass$new_sampleid,]

p_Data_noCH <- p_Data[as.character(p_Data$sample) %in% Ex_samp_id,]
##drop MGRB sample BAAUD; not in latest call
#samp_ID_match <- samp_ID_match[grep("BAAUD", samp_ID_match, invert = T)]
#p_Data_noCH <- p_Data_noCH[match(samp_ID_match, p_Data_noCH$sample),]
Ex_samp_id <- Ex_samp_id[match(p_Data_noCH$sample, Ex_samp_id)]
##filter out QC fail cases
fil_tab <- fil_tab[fil_tab$SAMPLE %in% Ex_samp_id,]

fil_tab_MMR_gene_var <- fil_tab[fil_tab$gene_symbol %in% MMR_genes,]
fil_tab_MMR_gene_var <- fil_tab_MMR_gene_var[fil_tab_MMR_gene_var$comb_score >= 5.6 & fil_tab_MMR_gene_var$VAF >= 0.35,]
write.table(fil_tab_MMR_gene_var, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/round2/MMR_genes_filt_VARIANTS.tsv", sep = "\t", row.names = F, quote = F)
