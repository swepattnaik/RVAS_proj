#fil_tab_noCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_noNAs_auto_4Aug.rds")
#fil_tab_noCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_noNAs_auto_6Aug.rds")
#fil_tab_noCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_auto_7Aug.rds")
#fil_tab_noCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_12Aug.rds")
#fil_tab_noCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_05Sept_rect.rds")
#fil_tab_noCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_05Sept_splice.rds")
#fil_tab_noCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_05Sept_rect_ASP.rds")
fil_tab_noCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_05Sept_rect_ASP.rds")
fil_tab_noCH <- fil_tab_noCH[!is.na(fil_tab_noCH$SAMPLE),]

get_coh_dist <- function(gene_sym, fil_db){
  res <- list()
  for(i in 1:length(gene_sym)){
    fil_db_genes <- fil_db[grep(paste("^", as.character(gene_sym[i]), "$", sep = ""), fil_db$gene_symbol),]
    fil_db_genes <- fil_db_genes[fil_db_genes$comb_score >= 5,]
    tot_var <- dim(fil_db_genes)[1]
    isks <- sum(ifelse(fil_db_genes$is_MGRB == 0, 1, 0))
    control <- tot_var - isks
    case_wt <- sum(fil_db_genes[fil_db_genes$is_MGRB == 0,]$comb_score)
    control_wt <- sum(fil_db_genes[fil_db_genes$is_MGRB == 1,]$comb_score)
    wt_diff <- case_wt - control_wt
case_call_auto <- paste(apply(as.data.frame(table(as.character(fil_db_genes[fil_db_genes$is_MGRB == 0,]$auto_call))) , 1 , paste , collapse = ":" ), collapse = ",")
    control_call_auto <- paste(apply(as.data.frame(table(as.character(fil_db_genes[fil_db_genes$is_MGRB == 1,]$auto_call))) , 1 , paste , collapse = ":" ), collapse = ",")
    case_vep_var <- paste(apply(as.data.frame(table(as.character(fil_db_genes[fil_db_genes$is_MGRB == 0,]$vep_consequence))) , 1 , paste , collapse = ":" ), collapse = ",")
    control_vep_var <- paste(apply(as.data.frame(table(as.character(fil_db_genes[fil_db_genes$is_MGRB == 1,]$vep_consequence))) , 1 , paste , collapse = ":" ), collapse = ",")
    res[[i]] <- cbind.data.frame("ISKS" = isks, "Control" = control, "gene" = gene_sym[i], "case_wt" = case_wt, "control_wt" = control_wt, "wt_diff" = wt_diff, "case_call_auto" = case_call_auto, "control_call_auto" = control_call_auto, "case_vep_var" = case_vep_var, "control_vep_var" = control_vep_var)
  }
  return(do.call("rbind.data.frame",res))
}


##Additive model
#Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_skat_wsing_load123_noCH_C5eqC4_noNA_gt_4Aug.rds")
#Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_skat_wsing_load123_noCH_C5eqC4_noNA_gt_6Aug.rds")
#Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_skat_wsing_load123_noCH_C5eqC4_nonmds_gt_isksrisc_12Aug.rds")
#Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_skat_wsing_load123_noCH_C5eqC4_nonmds_gt_isksrisc_12Aug_rerun.rds")
#Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_skat_wsing_load123_noCH_C5eqC4_nonmds_gt_isksrisc_sept05_rect.rds")
#Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_skat_wsing_load123_noCH_C5eqC4_nonmds_gt_isksrisc_sept05_splice.rds")
#Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_skat_wsing_load123_noCH_C5eqC4_nonmds_gt_isksrisc_sept05_rect_ASP.rds")
Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Exome_skat_wsing_load123_noCH_C5eqC4_nonmds_gt_isksrisc_sept05_rect_ASP_cmaf.rds")
df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
#Exome_pc123_srt_SKATO <- lapply(df_skat, function(x) x[order(x$pval_SKATO, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})

Exome_pc123_srt_SKAT_case_enr_nCH <- lapply(Exome_pc123_srt_SKAT, 
                                         function(x)get_coh_dist(x[,1], fil_db = fil_tab_noCH))
#saveRDS(Exome_pc123_srt_SKAT_case_enr_nCH, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_pc123_srt_SKAT_case_enr_nCH.rds", compress = T)
#saveRDS(Exome_pc123_srt_SKAT_case_enr_nCH, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_pc123_srt_SKAT_case_enr_nCH_Aug6.rds", compress = T)
#saveRDS(Exome_pc123_srt_SKAT_case_enr_nCH, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_pc123_srt_SKAT_case_enr_nCH_Aug7.rds", compress = T)
#saveRDS(Exome_pc123_srt_SKAT_case_enr_nCH, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_pc123_srt_SKAT_case_enr_nCH_iskrisc_Aug12.rds", compress = T)
#saveRDS(Exome_pc123_srt_SKAT_case_enr_nCH, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_pc123_srt_SKAT_case_enr_nCH_iskrisc_Aug12_rerun.rds", compress = T)
#saveRDS(Exome_pc123_srt_SKAT_case_enr_nCH, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_pc123_srt_SKAT_case_enr_nCH_iskrisc_05Sept_rect.rds", compress = T)
#saveRDS(Exome_pc123_srt_SKAT_case_enr_nCH, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_pc123_srt_SKAT_case_enr_nCH_iskrisc_05Sept_splice.rds", compress = T)
#saveRDS(Exome_pc123_srt_SKAT_case_enr_nCH, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_pc123_srt_SKAT_case_enr_nCH_iskrisc_05Sept_rect_ASP.rds", compress = T)
saveRDS(Exome_pc123_srt_SKAT_case_enr_nCH, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Exome_pc123_srt_SKAT_case_enr_nCH_iskrisc_05Sept_rect_ASP_cmaf.rds", compress = T)
