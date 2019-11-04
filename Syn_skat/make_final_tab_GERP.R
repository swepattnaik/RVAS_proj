##This script generates a summary of SKAT analysis by combining various attributes to better visualise case control
##variants with corresponding SKAT p-values
fil_tab_noCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Syn/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_05Sept_rect_ASP_syn.rds")
fil_tab_noCH <- fil_tab_noCH[!is.na(fil_tab_noCH$SAMPLE),]
fil_tab_noCH <- fil_tab_noCH[!is.na(fil_tab_noCH$GERP_SCORE),]
fil_tab_noCH <- fil_tab_noCH[fil_tab_noCH$GERP_SCORE > 1,]

get_coh_dist <- function(gene_sym, fil_db){
  res <- list()
  for(i in 1:length(gene_sym)){
    print(i)
    print(gene_sym[i])
    fil_db_genes <- fil_db[grep(paste("^", as.character(gene_sym[i]), "$", sep = ""), fil_db$gene_symbol),]
    if(dim(fil_db_genes)[1] != 0){
    #fil_db_genes <- fil_db_genes[fil_db_genes$comb_score >= 5,]
    ##intra cohort filter
    ftemp_tab_var_id <- unique(fil_db_genes$VARIANT)
    var_vec <- list()
    for(m in 1:length(ftemp_tab_var_id)){
      sam_gene_gt <- fil_db_genes[fil_db_genes$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,58)]
      sam_gene_gt <- unique(sam_gene_gt)
      maf_vec_cont <- sum(ifelse(is.na(as.numeric(sam_gene_gt$SAMPLE)), 1, 0))/(2*1572)
      
      maf_vec_case <- sum(ifelse(!is.na(as.numeric(sam_gene_gt$SAMPLE)) | 
                                   grepl("^CR", as.character(sam_gene_gt$SAMPLE)), 1, 0))/(2*1110)
      
      
      if(!(as.numeric(as.character(maf_vec_cont)) >= 0.001 | 
           as.numeric(as.character(maf_vec_case)) >= 0.0015)){
        
        var_vec[[m]] <- ftemp_tab_var_id[m]
      }
      else{ 
        next 
      }
    }
    }
    else{
      next
    }
    ##report generation  
   # var_vec_fin <- do.call("rbind.data.frame",var_vec)
    fil_db_genes <- fil_db_genes[fil_db_genes$VARIANT %in% unlist(var_vec),]
   # fil_db_genes$wt_case <- var_vec_fin[match(fil_db_genes$VARIANT, var_vec_fin$VAR), 2]
   # fil_db_genes$wt_cont <- var_vec_fin[match(fil_db_genes$VARIANT, var_vec_fin$VAR), 3]
    tot_var <- dim(fil_db_genes)[1]
    isks <- sum(ifelse(fil_db_genes$is_MGRB == 0, 1, 0))
    control <- tot_var - isks
    case_wt <- sum(fil_db_genes[fil_db_genes$is_MGRB == 0,]$GERP_SCORE)
    control_wt <- sum(fil_db_genes[fil_db_genes$is_MGRB == 1,]$GERP_SCORE)
    wt_diff <- case_wt - control_wt
    res[[i]] <- cbind.data.frame("ISKS" = isks, "Control" = control, "gene" = gene_sym[i], 
                                 "case_wt" = case_wt, "control_wt" = control_wt, "wt_diff" = wt_diff)
  }
  return(do.call("rbind.data.frame",res))
}

##SKAT output is used to add p-values to all enriched genes
Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Syn/Exome_skat_filt3_nCH_C5eqC4_nonmds_iskrisc_05Sept_rect_ASP_syn_GERP.rds")
df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
#Exome_pc123_srt_SKATO <- lapply(df_skat, function(x) x[order(x$pval_SKATO, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})

Exome_pc123_srt_SKAT_case_enr_nCH <- lapply(Exome_pc123_srt_SKAT, 
                                            function(x)get_coh_dist(x[,1], fil_db = fil_tab_noCH))
saveRDS(Exome_pc123_srt_SKAT_case_enr_nCH, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Syn/Exome_pc123_srt_SKAT_case_enr_nCH_iskrisc_05Sept_rect_ASP_cmaf_syn_GERP.rds", compress = T)
