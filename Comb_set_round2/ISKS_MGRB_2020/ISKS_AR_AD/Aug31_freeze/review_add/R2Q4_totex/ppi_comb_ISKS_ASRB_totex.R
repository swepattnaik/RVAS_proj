##copied from ~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_EPIT/isks_epit_cmp.R
##ISKSvMGRB, EPIvsMGRB(no familial breast), CAIRNSvsMGRB(QC1 pass) combine
#.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)

Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/Exome_skat_para_result_isks_totex_count_PC1234_ver4_clinrect_Aug31.rds")
df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})
Exome_pc123_srt_SKAT_df <- Exome_pc123_srt_SKAT[[4]]
Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_para_pc1234_SKAT_Enriched_ISKS_2020_Aug31.rds")
Exome_pc123_srt_SKAT_case_enr_nCH_df <- Exome_pc123_srt_SKAT_case_enr_nCH[[4]]

##Add p-value
#Exome_pc123_srt_SKAT_case_enr_nCH_pval <- Map(cbind.data.frame, Exome_pc123_srt_SKAT_case_enr_nCH, Exome_pc123_srt_SKAT)
Exome_pc123_srt_SKAT_case_enr_nCH_df_srt <- Exome_pc123_srt_SKAT_case_enr_nCH_df[match(as.character(Exome_pc123_srt_SKAT_df$symbol),as.character(Exome_pc123_srt_SKAT_case_enr_nCH_df$gene)),]

Exome_pc123_srt_SKAT_case_enr_nCH_pval <- cbind.data.frame(Exome_pc123_srt_SKAT_case_enr_nCH_df_srt,
                                                           Exome_pc123_srt_SKAT_df)
isks_mgrb_genes <- Exome_pc123_srt_SKAT_case_enr_nCH_pval

#Add PPI column for ISKS top genes (degree, rank_PPI, PPI_p_val_wt, OR)
#isks_mgrb_genes[,17:20] <- NA
#PPI_col <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/ppi_res_fil_final_SKATbin.tsv", 
#                      sep = "\t", header = T, stringsAsFactors = F)
#isks_mgrb_genes[,17:20] <- PPI_col[match(isks_mgrb_genes$gene, PPI_col$gene), c(18,20,22,25)]
##composite score added
isks_mgrb_genes[,17:36] <- NA

PPI_col <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/ppi_res_fil_final_SKATbin_comb_str_biog_totex.tsv", 
                      sep = "\t", header = T, stringsAsFactors = F)
#PPI_col$PPI_p_val_wt <- ifelse(PPI_col$degree == 0, 1, PPI_col$PPI_p_val_wt)
#isks_mgrb_genes[,17:23] <- PPI_col[match(isks_mgrb_genes$gene, PPI_col$gene), c(15,17:22)]
isks_mgrb_genes[,17:36] <- PPI_col[match(as.character(isks_mgrb_genes$gene), PPI_col$gene), c(17,19:37)]
#colnames(isks_mgrb_genes)[17:23] <- c("degree", "rank_PPI", "PPI_p_val_wt", "new_comp_score","t_stat", "pval_t.test",
#                                      "comb_fisher")
colnames(isks_mgrb_genes)[17:36] <- colnames(PPI_col)[c(17,19:37)]
isks_mgrb_genes_top <- isks_mgrb_genes[isks_mgrb_genes$wt_diff > 0 & isks_mgrb_genes$pval_SKATbin < 0.1,]

##Epithelial data
# Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/round2/Exome_skat_para_result_EPIT_combset2020_uni_MAF_PC1234_ver4_no_fbrca.rds")
# 
# df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
# Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
# Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})
# Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/round2/Exome_para_pc123_SKAT_Enriched_EPIT_2020_no_fbrca.rds")
# ##Add p-value
# Exome_pc123_srt_SKAT_case_enr_nCH_pval <- Map(cbind.data.frame, Exome_pc123_srt_SKAT_case_enr_nCH, Exome_pc123_srt_SKAT)
# 
# epit_mgrb_genes <- Exome_pc123_srt_SKAT_case_enr_nCH_pval[[4]]
# epit_mgrb_genes_top <- epit_mgrb_genes[epit_mgrb_genes$wt_diff > 0 & epit_mgrb_genes$pval_SKATbin < 0.1,]

##Cairns data
Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/Exome_skat_para_result_asrb_totex_count_PC1234_ver4_clinrect_Aug31.rds")

df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})
Exome_pc123_srt_SKAT_df <- Exome_pc123_srt_SKAT[[4]]
Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CAIRNS/round2/ASRB_AR_AD/SKAT/Exome_para_pc123_SKAT_Enriched_CAIRNS_2020.rds")
##Add p-value
Exome_pc123_srt_SKAT_case_enr_nCH_df <- Exome_pc123_srt_SKAT_case_enr_nCH[[4]]

##Add p-value
#Exome_pc123_srt_SKAT_case_enr_nCH_pval <- Map(cbind.data.frame, Exome_pc123_srt_SKAT_case_enr_nCH, Exome_pc123_srt_SKAT)
Exome_pc123_srt_SKAT_case_enr_nCH_df_srt <- Exome_pc123_srt_SKAT_case_enr_nCH_df[match(as.character(Exome_pc123_srt_SKAT_df$symbol),as.character(Exome_pc123_srt_SKAT_case_enr_nCH_df$gene)),]

Exome_pc123_srt_SKAT_case_enr_nCH_pval <- cbind.data.frame(Exome_pc123_srt_SKAT_case_enr_nCH_df_srt,
                                                           Exome_pc123_srt_SKAT_df)

cairns_mgrb_genes <- Exome_pc123_srt_SKAT_case_enr_nCH_pval
cairns_mgrb_genes_top <- cairns_mgrb_genes[cairns_mgrb_genes$wt_diff > 0 & cairns_mgrb_genes$pval_SKATbin < 0.1,]


##combine
#isks_epit_mgrb_genes_top <- isks_mgrb_genes[isks_mgrb_genes$gene %in% epit_mgrb_genes_top$gene,]
isks_cairns_mgrb_genes <- isks_mgrb_genes[isks_mgrb_genes$gene %in% cairns_mgrb_genes_top$gene,]
isks_cairns_mgrb_genes_top <- rbind.data.frame(isks_mgrb_genes_top, isks_cairns_mgrb_genes)
isks_cairns_mgrb_genes_top <- unique(isks_cairns_mgrb_genes_top)

#isks_cairns_mgrb_genes_top[,21:36] <- cairns_mgrb_genes[match(isks_cairns_mgrb_genes_top$gene, cairns_mgrb_genes$gene),1:16]
#isks_cairns_mgrb_genes_top[,24:39] <- cairns_mgrb_genes[match(isks_cairns_mgrb_genes_top$gene, cairns_mgrb_genes$gene),1:16]

isks_cairns_mgrb_genes_top[,37:52] <- cairns_mgrb_genes[match(isks_cairns_mgrb_genes_top$gene, cairns_mgrb_genes$gene),1:16]


##remove all genes with only one variant in ISKS
isks_cairns_mgrb_genes_top <- isks_cairns_mgrb_genes_top[isks_cairns_mgrb_genes_top$ISKS > 1,]
#isks_cairns_mgrb_genes_top <- isks_cairns_mgrb_genes_top[,-c(13,15:16,36,38:39)]
#isks_cairns_mgrb_genes_top <- isks_cairns_mgrb_genes_top[,-c(13,15:16,57:58,60,62:63)]


write.table(isks_cairns_mgrb_genes_top, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/SKAT/isks_cyto_ppifisher_ASRB_comb_genes_totex_2021.tsv",
            sep = "\t", row.names = F, quote = F)
#write.table(isks_cairns_mgrb_genes_top, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/isks_cyto_ppifisher_ASRB_comb_genes_top_new_rand_comp_score_Aug31.tsv",
#            sep = "\t", row.names = F, quote = F)
#write.table(isks_cairns_mgrb_genes_top, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/isks_ppifisher_ASRB_comb_genes_top_new_rand_comp_score_Aug31.tsv",
#            sep = "\t", row.names = F, quote = F)