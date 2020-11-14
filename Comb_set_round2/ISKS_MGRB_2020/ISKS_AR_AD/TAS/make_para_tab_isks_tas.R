##This script generates a summary of SKAT analysis by combining various attributes to better visualise case control
##variants with corresponding SKAT p-values
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
fil_tab_noCH <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_all_fields_rnd3.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
Ex_samp_id <- unique(fil_tab_noCH$SAMPLE)


dup_samp <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/fin_samp_dup_drop.txt", sep = "", header = T,
                       stringsAsFactors = F)
`%nin%` = Negate(`%in%`)
Ex_samp_id <- Ex_samp_id[Ex_samp_id %nin% dup_samp$x]


######get phenotype data to control for age and sex and PC's; 

QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv", 
                      header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %in% "TRUE",]
QC2_dat_pass$isFemale <- ifelse(QC2_dat_pass$f_stat < 0.2, 1, 
                                ifelse(QC2_dat_pass$f_stat > 0.8, 0, 2))
p_Data <- read.table("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_combset_pca.scores_clustered.tsv", 
                     header = T, sep = "\t", stringsAsFactors = F)
p_Data <- p_Data[p_Data$sample %in% QC2_dat_pass$new_sampleid,]
p_Data_noCH <- p_Data[as.character(p_Data$sample) %in% Ex_samp_id,]
Ex_samp_id <- Ex_samp_id[match(p_Data_noCH$sample, Ex_samp_id)]
##filter out QC fail cases
fil_tab_noCH <- fil_tab_noCH[fil_tab_noCH$SAMPLE %in% Ex_samp_id,]

library(readxl)
comb_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/PID_Master_file_290420.xlsx", sheet = 1, col_types = c("list"))
comb_pheno <- as.data.frame(comb_pheno)
comb_pheno1 <- sapply(comb_pheno, unlist)
colnames(comb_pheno1) <- colnames(comb_pheno)
comb_pheno <- comb_pheno1
comb_pheno <- as.data.frame(comb_pheno, stringsAsFactors = F)
comb_pheno <- unique(comb_pheno)
comb_pheno <- comb_pheno[!is.na(comb_pheno$pid),]
comb_pheno$`age at dateExtracted` <- as.numeric(comb_pheno$`age at dateExtracted`)
comb_pheno$AgeatSarcoma <- as.numeric(comb_pheno$AgeatSarcoma)
comb_pheno$SubjectAgeCancer <- as.numeric(comb_pheno$SubjectAgeCancer)
#genomicclass
comb_pheno_no_comp <- as.character(comb_pheno[comb_pheno$genomicclass %nin% "TAS",]$pmn)
##all TAS cases are annotated in Mandy's data
not_annot <- Ex_samp_id[Ex_samp_id %nin% comb_pheno[comb_pheno$genomicclass %in% "TAS",]$pmn & !grepl("^[ABZ]", Ex_samp_id)]

##filter out QC fail cases
fil_tab_noCH <- fil_tab_noCH[fil_tab_noCH$SAMPLE %nin% comb_pheno_no_comp,]
fil_tab_noCH <- fil_tab_noCH[fil_tab_noCH$SAMPLE %nin% not_annot,]
fil_tab_noCH <- fil_tab_noCH[fil_tab_noCH$SAMPLE %in% Ex_samp_id,]
length(unique(fil_tab_noCH$SAMPLE))
length(unique(fil_tab_noCH[fil_tab_noCH$is_case == 1,]$SAMPLE))

get_coh_dist_para <- function(gene_sym, fil_db){
  #res <- list()
  # for(i in 1:length(gene_sym)){
  #   print(i)
  print(as.character(gene_sym))
  fil_db_genes <- fil_db[grep(paste("^", as.character(gene_sym), "$", sep = ""), fil_db$gene_symbol),]
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
      #maf_vec_cont <- sum(ifelse(is.na(as.numeric(sam_gene_gt$SAMPLE)), 1, 0))/((354 + 3205)*2)
      maf_vec_cont <- length(grep("^[ABZ]", sam_gene_gt$SAMPLE))/((354 + 3205)*2)
      maf_vec_case <- sum(ifelse(!is.na(as.numeric(sam_gene_gt$SAMPLE)) | 
                                   grepl("^CR|^LK", as.character(sam_gene_gt$SAMPLE)), 1, 0))/((354 + 3205)*2)
      ##MAF filter = 5/(1661*2); change to 3/(1661*2)
      ##changed to cohort MAF threshold : 3/((354 + 3209)*2) 
      if(!(as.numeric(as.character(maf_vec_cont)) >= 0.00043 | 
           as.numeric(as.character(maf_vec_case)) >= 0.00043)){
        var_vec[[m]] <- ftemp_tab_var_id[m]
      }
      else{ 
        next 
      }
    }
  }
  ##report generation  
  fil_db_genes <- fil_db_genes[fil_db_genes$VARIANT %in% unlist(var_vec),]
  fil_db_genes$is_case <- ifelse(grepl("^ISKS", fil_db_genes$set), 1, 0)
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
  res <- cbind.data.frame("ISKS_tas" = isks, "MGRB" = control, "gene" = gene_sym, 
                          "case_wt" = case_wt, "control_wt" = control_wt, "wt_diff" = wt_diff, 
                          "case_call_auto" = case_call_auto, "control_call_auto" = control_call_auto, 
                          "case_vep_var" = case_vep_var, "control_vep_var" = control_vep_var,
                          "maf_filt_var" = paste(unlist(var_vec), collapse = ","))
  
  return(res)
}


##SKAT output is used to add p-values to all enriched genes
# Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Exome_skat_para_result_isks_combset2020_uni_MAF_PC1234_ver4_complex.rds")
# df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
# Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
# Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})

##generate table in parallel
library(doParallel)
library(doMC)
registerDoMC(30)

res20 <- list()


#for(k in 1:length(unique(fil_tab_noCH$gene_symbol))){
res_list <- list()
#  genes <-  as.character(Exome_pc123_srt_SKAT[[k]][,1])
genes <- unique(fil_tab_noCH$gene_symbol)
system.time(res_list <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
{get_coh_dist_para(genes[i], fil_tab_noCH)})
res20 <- do.call("rbind.data.frame", res_list)
#}

saveRDS(res20, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/TAS/Exome_para_ISKS_TAS.rds", compress = T)
