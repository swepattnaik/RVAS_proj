.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)

repset_pois_563_PID <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_POIS_filt_C345_tcga_nosarc_berg_mgrb_PID.tsv",
                         sep = "\t", header = T, stringsAsFactors = F)
repset_pois_563 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_POIS_filt_C345_tcga_nosarc_berg_mgrb.tsv",
                                   sep = "\t", header = T, stringsAsFactors = F)

get_coh_dist_para <- function(gene_sym, fil_db, tot_coh_size){
  #res <- list()
  # for(i in 1:length(gene_sym)){
  #   print(i)
  print(as.character(gene_sym))
  # fil_db_genes <- fil_db[grep(paste("^", as.character(gene_sym), "$", sep = ""), fil_db$gene_symbol),]
  fil_db_genes <- fil_db[fil_db$gene_symbol %in% as.character(gene_sym),]
  fil_db_genes <- fil_db_genes[fil_db_genes$auto_call %nin% "C3",]
 # fil_db_genes <- fil_db_genes[fil_db_genes$comb_score >= 5.6 & fil_db_genes$VAF >= 0.35,]
      if(is.null(dim(fil_db_genes)) | dim(fil_db_genes)[1] < 1 ){
        res <- NULL
        return(res)
      }
      else{
  ##intra cohort filter; uniform MAF filter
  ftemp_tab_var_id <- unique(fil_db_genes$VARIANT)
  if(length(ftemp_tab_var_id) == 0 ){ ##added on apr.22
    res <- NULL
    return(res)
  }
  else {
    var_vec <- list()
    for(m in 1:length(ftemp_tab_var_id)){
      sam_gene_gt <- fil_db_genes[fil_db_genes$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,127:128)]
      sam_gene_gt <- unique(sam_gene_gt)
      #maf_vec_cont <- sum(ifelse(is.na(as.numeric(sam_gene_gt$SAMPLE)), 1, 0))/(tot_coh_size*2)
      ##it should be for MGRB and to not penalise CR's and LK's:Oct18-2020
      maf_vec_cont <- length(grep("^[ABZ]", sam_gene_gt$SAMPLE))/(2*tot_coh_size)
      maf_vec_case <- length(grep("^[ABZ]", sam_gene_gt$SAMPLE, invert = T))/(2*tot_coh_size)
      #maf_vec_case <- sum(ifelse(!is.na(as.numeric(sam_gene_gt$SAMPLE)) | 
      #                             grepl("^CR|^LK", as.character(sam_gene_gt$SAMPLE)), 1, 0))/(tot_coh_size*2)
      ##MAF filter = 5/(1661*2); change to 3/(1661*2)
      ##changed to cohort MAF threshold : 3/((563 + 3205)*2) 
      maf_cutoff = 3.5/(tot_coh_size*2)
      if(!(as.numeric(as.character(maf_vec_cont)) >= maf_cutoff | 
           as.numeric(as.character(maf_vec_case)) >= maf_cutoff)){
        var_vec[[m]] <- ftemp_tab_var_id[m]
      }
      else{ 
        res <- NULL
        return(res) 
      }
    }
  }
  ##report generation  
  fil_db_genes <- fil_db_genes[fil_db_genes$VARIANT %in% unlist(var_vec),]
  fil_db_genes$is_case <- ifelse(grepl("^MGRB", fil_db_genes$set), 0, 1)
  tot_var <- dim(fil_db_genes)[1]
  isks <- sum(fil_db_genes$is_case)
  control <- tot_var - isks
  case_wt <- sum(fil_db_genes[fil_db_genes$is_case == 1,]$comb_score)
  control_wt <- sum(fil_db_genes[fil_db_genes$is_case == 0,]$comb_score)
  wt_diff <- case_wt - control_wt
  case_call_auto <- paste(apply(as.data.frame(table(as.character(fil_db_genes[fil_db_genes$is_case == 1,]$auto_call))) , 1 , paste , collapse = ":" ), collapse = ",")
  control_call_auto <- paste(apply(as.data.frame(table(as.character(fil_db_genes[fil_db_genes$is_case == 0,]$auto_call))) , 1 , paste , collapse = ":" ), collapse = ",")
  case_vep_var <- paste(apply(as.data.frame(table(as.character(fil_db_genes[fil_db_genes$is_case == 1,]$vep_consequence))) , 1 , paste , collapse = ":" ), collapse = ",")
  control_vep_var <- paste(apply(as.data.frame(table(as.character(fil_db_genes[fil_db_genes$is_case == 0,]$vep_consequence))) , 1 , paste , collapse = ":" ), collapse = ",")
  res <- cbind.data.frame("Sarc_comb" = isks, "MGRB" = control, "gene" = gene_sym, 
                          "case_wt" = case_wt, "control_wt" = control_wt, "wt_diff" = wt_diff, 
                          "case_call_auto" = case_call_auto, "control_call_auto" = control_call_auto, 
                          "case_vep_var" = case_vep_var, "control_vep_var" = control_vep_var,
                          "maf_filt_var" = paste(unlist(var_vec), collapse = ","))
  
  return(res)
}
}

##generate table in parallel
library(doParallel)
library(doMC)
registerDoMC(30)

para_make_tab <- function(df_inp, fname , tot_coh_size){
  res20 <- list()
  res_list <- list()
  
  genes <- unique(as.character(df_inp[,9]))
  system.time(res_list <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
  {get_coh_dist_para(genes[i], df_inp, tot_coh_size)})
  res20 <- do.call("rbind.data.frame", res_list)
  saveRDS(res20, paste0("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/",
                        fname,".rds"), compress = T)
}

para_make_tab(repset_pois_563_PID, "repset_563_POIS_PID_noC3", tot_coh_size = 3768) #563 + 3205
para_make_tab(repset_pois_563, "repset_563_POIS_noC3", tot_coh_size = 3768) #563 + 3205