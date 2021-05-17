##This script generates a summary of SKAT analysis by combining various attributes to better visualise case control
##variants with corresponding SKAT p-values
fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CVD_DIANE/round2/CVD_AR_AD/all_CVD_mgrb_skatinp_combset2021_clin_C3C4C5_NFE0002_AD_all_fields_rnd3.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)

##Additive model
Ex_samp_id <- unique(fil_tab$SAMPLE)

#remove duplicates(not needed)
#dup_samp <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/fin_samp_dup_drop.txt", sep = "", header = T,
#                       stringsAsFactors = F)
#`%nin%` = Negate(`%in%`)
#Ex_samp_id <- Ex_samp_id[Ex_samp_id %nin% dup_samp$x]
######get phenotype data to control for age and sex and PC's; 

QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv",
                      header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %in% "TRUE",]
QC2_dat_pass$isFemale <- ifelse(QC2_dat_pass$f_stat < 0.2, 1,
                                ifelse(QC2_dat_pass$f_stat > 0.8, 0, 2))
QC2_dat_pass_mgrb <- QC2_dat_pass[grepl("^[ABZ]", QC2_dat_pass$new_sampleid),]
QC2_dat_pass_mgrb <- QC2_dat_pass_mgrb[,c(1,21)]

##Pheno data
p_Data <- read.delim("~/RVAS/CVD_Diane/pop_PCA/MGRB_CVD_1000G_combset_pca.scores_clustered.tsv", header = T, sep = "\t", stringsAsFactors = F)
p_Data_noCH <- p_Data[as.character(p_Data$sample) %in% Ex_samp_id,]
##drop MGRB sample BAAUD; not in latest call
#samp_ID_match <- samp_ID_match[grep("BAAUD", samp_ID_match, invert = T)]
#p_Data_noCH <- p_Data_noCH[match(samp_ID_match, p_Data_noCH$sample),]

##remove from gene matrix
#col_rm <- which(colnames(D_tab) ==  "BAAUD")
#D_tab <- D_tab[,-c(col_rm)]
##Add gender information
Ex_samp_id_df <- as.data.frame(Ex_samp_id)
Ex_samp_id_df$Ex_samp_id <- as.character(Ex_samp_id_df$Ex_samp_id)
Ex_samp_id_df$gender <- QC2_dat_pass_mgrb[match(Ex_samp_id_df$Ex_samp_id, QC2_dat_pass_mgrb$new_sampleid), 2]
#p_Data_noCH$gender <- QC2_dat_pass_mgrb[match(p_Data_noCH$sample, QC2_dat_pass_mgrb$new_sampleid), 2]

##Add CVD phenotype
CVD_pheno <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/CVD_MGRB_2021/CVD_AR_AD/CVD_pheno.tsv",
                        sep = "\t", header = T, stringsAsFactors = F)
CVD_pheno$Gender <- ifelse(CVD_pheno$Gender %in% "F", 1, 0)
Ex_samp_id_df$gender_cvd <- CVD_pheno[match(Ex_samp_id_df$Ex_samp_id,CVD_pheno$Sample.ID), 2]
Ex_samp_id_df$gender <- ifelse(is.na(Ex_samp_id_df$gender), Ex_samp_id_df$gender_cvd, Ex_samp_id_df$gender)

p_Data_noCH$gender <- Ex_samp_id_df[match(p_Data_noCH$sample, Ex_samp_id_df$Ex_samp_id),2]

##remove cases with no gender (24 Cairns lost)
p_Data_noCH <- p_Data_noCH[!is.na(p_Data_noCH$gender),]
Ex_samp_id <- Ex_samp_id[match(p_Data_noCH$sample, Ex_samp_id)]
##filter out QC fail cases
fil_tab <- fil_tab[fil_tab$SAMPLE %in% Ex_samp_id,]

#write.table(fil_tab, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CAIRNS/all_cairns_mgrb_combset2020_variants_filt_all_fields_rmdup_fin.tsv",
#            sep = "\t", row.names = F, quote = F)

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
        toMatch <- c("-", "^F")
        maf_vec_cont <- dim(sam_gene_gt[!grepl(paste(toMatch,collapse="|"),sam_gene_gt$SAMPLE),])[1]/((117 + 3205)*2)
        maf_vec_case <- dim(sam_gene_gt[grepl(paste(toMatch,collapse="|"),sam_gene_gt$SAMPLE),])[1]/((117 + 3205)*2)
        ##MAF filter = 5/(1661*2); change to 3/(1662*2)
        ##changed to cohort MAF threshold : 3/((117 + 3205)*2) 
        if(!(as.numeric(as.character(maf_vec_cont)) >= 0.00046 | 
             as.numeric(as.character(maf_vec_case)) >= 0.00046)){
          var_vec[[m]] <- ftemp_tab_var_id[m]
        }
        else{ 
          next 
        }
      }
    }
    ##report generation  
    fil_db_genes <- fil_db_genes[fil_db_genes$VARIANT %in% unlist(var_vec),]
    #this was not performed in the DT_exome_clean script for EPIT
    toMatch <- c("-", "^F")
    fil_db_genes$is_case <- ifelse(grepl(paste(toMatch,collapse="|"), fil_db_genes$SAMPLE), 1, 0)
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
    res_df <- cbind.data.frame("CVD" = isks, "MGRB" = control, "gene" = gene_sym, 
                                 "case_wt" = case_wt, "control_wt" = control_wt, "wt_diff" = wt_diff, 
                                 "case_call_auto" = case_call_auto, "control_call_auto" = control_call_auto, 
                                 "case_vep_var" = case_vep_var, "control_vep_var" = control_vep_var,
                                 "maf_filt_var" = paste(unlist(var_vec), collapse = ","))
    
  return(res_df)
}


##SKAT output is used to add p-values to all enriched genes
Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CVD_DIANE/round2/CVD_AR_AD/SKAT/Exome_skat_para_result_CVD_combset2020_uni_MAF_PC1234_ver4.rds")
df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
#Exome_pc123_srt_SKATO <- lapply(df_skat, function(x) x[order(x$pval_SKATO, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})

#Exome_pc123_srt_SKAT_case_enr_nCH <- lapply(Exome_pc123_srt_SKAT, 
 #                                           function(x)get_coh_dist(x[,1], fil_db = fil_tab_noCH))

##Parallelised SKAT  
library(doParallel)
library(doMC)
registerDoMC(30)

res20 <- list()

#genes <- unique(fil_tab_noCH$gene_symbol)
for(k in 1:length(Exome_pc123_srt_SKAT)){
  res_list <- list()
  genes <-  as.character(Exome_pc123_srt_SKAT[[k]][,1])
system.time(res_list <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
{get_coh_dist_para(genes[i], fil_tab)})
res20[[k]] <- do.call("rbind.data.frame", res_list)
}


saveRDS(res20, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CVD_DIANE/round2/CVD_AR_AD/SKAT/Exome_para_pc123_SKAT_Enriched_CVD_2021.rds", compress = T)
