##This script generates a summary of SKAT analysis by combining various attributes to better visualise case control
##variants with corresponding SKAT p-values
fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/all_epit_mgrb_combset2020_variants_filt_all_fields_risk_factor.tsv",
                      sep = "\t", header = T, stringsAsFactors = F)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]

dim(fil_tab)

##Additive model
Ex_samp_id <- unique(fil_tab$SAMPLE)
Ex_var_id <- unique(fil_tab$VARIANT)

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
fil_tab_noCH <- fil_tab[fil_tab$SAMPLE %in% Ex_samp_id,]

get_coh_dist_para <- function(gene_sym, fil_db){
  #res <- list()
 # for(i in 1:length(gene_sym)){
 #   print(i)
    print(as.character(gene_sym))
    fil_db_genes <- fil_db[grep(paste("^", as.character(gene_sym), "$", sep = ""), fil_db$gene_symbol),]
    fil_db_genes <- fil_db_genes[fil_db_genes$comb_score >= 3 & fil_db_genes$VAF >= 0.35,]
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
        sam_gene_gt <- fil_db_genes[fil_db_genes$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,125:126)]
        sam_gene_gt <- unique(sam_gene_gt)
        maf_vec_cont <- dim(sam_gene_gt[grepl("^[ABZ]",sam_gene_gt$SAMPLE),])[1]/((920 + 3209)*2)
        maf_vec_case <- dim(sam_gene_gt[!grepl("^[ABZ]",sam_gene_gt$SAMPLE),])[1]/((920 + 3209)*2)
        ##MAF filter = 5/(1661*2); change to 3/(1662*2)
        ##changed to cohort MAF threshold : 3/((920 + 3209)*2) 
        if(!(as.numeric(as.character(maf_vec_cont)) >= 0.00037 | 
             as.numeric(as.character(maf_vec_case)) >= 0.00037)){
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
    fil_db_genes$is_case <- ifelse(grepl("^[ABZ]", fil_db_genes$SAMPLE), 0, 1)
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
    res_df <- cbind.data.frame("EPIT" = isks, "Control" = control, "gene" = gene_sym, 
                                 "case_wt" = case_wt, "control_wt" = control_wt, "wt_diff" = wt_diff, 
                                 "case_call_auto" = case_call_auto, "control_call_auto" = control_call_auto, 
                                 "case_vep_var" = case_vep_var, "control_vep_var" = control_vep_var,
                                 "maf_filt_var" = paste(unlist(var_vec), collapse = ","))
    
  return(res_df)
}


##SKAT output is used to add p-values to all enriched genes
Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/Exome_skat_para_result_EPIT_combset2020_uni_MAF_PC1234_ver4_rect_test.rds")
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
{get_coh_dist_para(genes[i], fil_tab_noCH)})
res20[[k]] <- do.call("rbind.data.frame", res_list)
}


#saveRDS(Exome_pc123_srt_SKAT_case_enr_nCH, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/Exome_pc123_SKAT_Enriched_combset2020.rds", compress = T)
saveRDS(res20, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/Exome_para_pc123_SKAT_Enriched_EPIT_2020.rds", compress = T)
