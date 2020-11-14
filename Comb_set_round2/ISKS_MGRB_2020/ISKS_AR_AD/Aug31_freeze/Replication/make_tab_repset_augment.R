##This script generates a summary of SKAT analysis by combining various attributes to better visualise case control
##variants with corresponding SKAT p-values
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)
##PID repset
#repset_563 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/all_PID_genes_tcga_nosarc_berg_mgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_all_fields_rnd3_Aug31.tsv",
#                           sep = "\t", header = T, stringsAsFactors = F)
#Ex_samp_id <- unique(fil_tab_noCH$SAMPLE)
##Augmented repsets PID
repset_663 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_tcga_nosarc_berg_isks100_mgrb_PID.tsv",
                         sep = "\t", header = T, stringsAsFactors = F)
repset_563 <- repset_663[repset_663$set %nin% "ISKS_AR_AD",]
isks_min100_663 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_isks_mgrb_minus100.tsv",
                              sep = "\t", header = T, stringsAsFactors = F)
repset_763 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_tcga_nosarc_berg_isks200_mgrb_PID.tsv",
                         sep = "\t", header = T, stringsAsFactors = F)
isks_min200_763 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_isks_mgrb_minus200.tsv",
                              sep = "\t", header = T, stringsAsFactors = F)
repset_963 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_tcga_nosarc_berg_isks400_mgrb_PID.tsv",
                         sep = "\t", header = T, stringsAsFactors = F)
isks_min400_963 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_isks_mgrb_minus400.tsv",
                              sep = "\t", header = T, stringsAsFactors = F)
isks_mgrb_all <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_isks_mgrb_maf_filt.tsv",
                            sep = "\t", header = T, stringsAsFactors = F)


get_coh_dist_para <- function(gene_sym, fil_db, tot_coh_size){
  #res <- list()
  # for(i in 1:length(gene_sym)){
  #   print(i)
  print(as.character(gene_sym))
  # fil_db_genes <- fil_db[grep(paste("^", as.character(gene_sym), "$", sep = ""), fil_db$gene_symbol),]
  fil_db_genes <- fil_db[fil_db$gene_symbol %in% as.character(gene_sym),]
  fil_db_genes <- fil_db_genes[fil_db_genes$comb_score >= 5.6 & fil_db_genes$VAF >= 0.35,]
  #    if(is.null(dim(fil_db_genes)) | dim(fil_db_genes)[1] < 1 ){
  #      next
  #    }
  #    else{
  ##intra cohort filter; uniform MAF filter
  ftemp_tab_var_id <- unique(fil_db_genes$VARIANT)
  if(length(ftemp_tab_var_id) == 0 ){ ##added on apr.22
    next
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
        next 
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

para_make_tab(repset_563, "repset_563", tot_coh_size = 3768) #563 + 3205
para_make_tab(repset_663, "repset_663", tot_coh_size = 3868) #663 + 3205
para_make_tab(repset_763, "repset_763", tot_coh_size = 3968) #763 + 3205
para_make_tab(repset_963, "repset_963", tot_coh_size = 4168) #963 + 3205

##isks complements
para_make_tab(isks_min100_663, "isks_min100", tot_coh_size = 4749) # 1544 + 3205
para_make_tab(isks_min200_763, "isks_min200", tot_coh_size = 4649) # 1444 + 3205
para_make_tab(isks_min400_963, "isks_min400", tot_coh_size = 4449) # 1244 + 3205
para_make_tab(isks_mgrb_all, "isks_mgrb_all", tot_coh_size = 4849) # 1644 + 3205