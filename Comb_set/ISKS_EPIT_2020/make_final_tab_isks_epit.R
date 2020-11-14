##This script generates a summary of SKAT analysis by combining various attributes to better visualise case control
##variants with corresponding SKAT p-values
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

fil_tab_isks <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/all_isksmgrb_combset2020_variants_filt_all_fields.tsv",
                           sep = "\t", stringsAsFactors = F, header = T)
##fil_tab_isks has risk_factor assigned as C5
fil_tab_isks <- fil_tab_isks[!is.na(fil_tab_isks$SAMPLE),]
fil_tab_isks_set <- fil_tab_isks[fil_tab_isks$is_case == 1,]
#dim(fil_tab_isks_set)

fil_tab_epit <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/all_epit_mgrb_combset2020_variants_filt_all_fields_risk_factor.tsv",
                           sep = "\t", stringsAsFactors = F, header = T)
fil_tab_epit <- fil_tab_epit[!is.na(fil_tab_epit$SAMPLE),]
fil_tab_epit_set <- fil_tab_epit[fil_tab_epit$is_case == 1,]
#dim(fil_tab_epit_set)
fil_tab_epit_set$is_case <- 0

fil_tab <- rbind.data.frame(fil_tab_isks_set, fil_tab_epit_set)

##Additive model
Ex_samp_id <- unique(fil_tab$SAMPLE)
Ex_var_id <- unique(fil_tab$VARIANT)

######get phenotype data to control for age and sex and PC's; 

QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv", 
                      header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %in% "TRUE",]
QC2_dat_pass$isFemale <- ifelse(QC2_dat_pass$f_stat < 0.2, 1, 
                                ifelse(QC2_dat_pass$f_stat > 0.8, 0, 2))
QC2_dat_pass_isks <- QC2_dat_pass[!grepl("^[ABZ]", QC2_dat_pass$new_sampleid),]
QC2_dat_pass_isks <- QC2_dat_pass_isks[,c(1,20:21)]

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

QC2_epit_isks_pass <- rbind.data.frame(QC2_epit_pass, QC2_dat_pass_isks)

p_Data <- read.delim("~/RVAS/comb_set_2020/epit_isks_mgrb_pop_PCA/MGRB_ISKS_EPIT_1000G_combset_pca.scores_clustered.tsv", 
                     header = T, sep = "\t", stringsAsFactors = F)

p_Data_noCH <- p_Data[as.character(p_Data$sample) %in% Ex_samp_id,]
##retain QC pass samples## 66 EPIT samples failed QC
p_Data_noCH <- p_Data_noCH[p_Data_noCH$sample %in% QC2_epit_isks_pass$new_sampleid,]
Ex_samp_id <- Ex_samp_id[match(p_Data_noCH$sample, Ex_samp_id)]
##filter out QC fail cases
fil_tab_noCH <- fil_tab[fil_tab$SAMPLE %in% Ex_samp_id,]

get_coh_dist <- function(gene_sym, fil_db){
  res <- list()
  for(i in 1:length(gene_sym)){
    print(i)
    print(as.character(gene_sym[i]))
    fil_db_genes <- fil_db[grep(paste("^", as.character(gene_sym[i]), "$", sep = ""), fil_db$gene_symbol),]
    fil_db_genes <- fil_db_genes[fil_db_genes$comb_score >= 3 & fil_db_genes$VAF >= 0.35,]
#    if(is.null(dim(fil_db_genes)) | dim(fil_db_genes)[1] < 1 ){
#      next
#    }
#    else{
    ##intra cohort filter; uniform MAF filter
    ftemp_tab_var_id <- unique(fil_db_genes$VARIANT)
    var_vec <- list()
    for(m in 1:length(ftemp_tab_var_id)){
      sam_gene_gt <- fil_db_genes[fil_db_genes$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,11,125:126)]
      sam_gene_gt <- unique(sam_gene_gt)
      maf_vec_cont <- sum(ifelse(is.na(as.numeric(sam_gene_gt$SAMPLE)), 1, 0))/(2*920)
      maf_vec_case <- sum(ifelse(!is.na(as.numeric(sam_gene_gt$SAMPLE)) | 
                                   grepl("^CR|^LK", as.character(sam_gene_gt$SAMPLE)), 1, 0))/(2*1661)
      print(ftemp_tab_var_id[m])
      print(as.numeric(as.character(maf_vec_cont)))
      print(as.numeric(as.character(maf_vec_case)))
      ##MAF filter = 5/(1661*2); change to 3/(1662*2); change to entire cohort MAF = 3/(2581)
      if(!(as.numeric(as.character(maf_vec_cont)) >= 0.00116234 | 
           as.numeric(as.character(maf_vec_case)) >= 0.00116234)){
           var_vec[[m]] <- ftemp_tab_var_id[m]
                                                            }
      else{ 
          next 
          }
                                             }
  ##report generation  
    fil_db_genes <- fil_db_genes[fil_db_genes$VARIANT %in% unlist(var_vec),]
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
    res[[i]] <- cbind.data.frame("ISKS" = isks, "Control" = control, "gene" = gene_sym[i], 
                                 "case_wt" = case_wt, "control_wt" = control_wt, "wt_diff" = wt_diff, 
                                 "case_call_auto" = case_call_auto, "control_call_auto" = control_call_auto, 
                                 "case_vep_var" = case_vep_var, "control_vep_var" = control_vep_var,
                                 "maf_filt_var" = paste(unlist(var_vec), collapse = ","))
  }
  return(do.call("rbind.data.frame",res))
}
#}

##SKAT output is used to add p-values to all enriched genes
#Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/Exome_skat_para_result_isks_combset2020_uni_MAF_PC1234_ver4_rect.rds")
Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/Exome_skat_para_result_isks_epit_combset2020_PC1234_ver4_rect_test.rds")
df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
#Exome_pc123_srt_SKATO <- lapply(df_skat, function(x) x[order(x$pval_SKATO, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})

Exome_pc123_srt_SKAT_case_enr_nCH <- lapply(Exome_pc123_srt_SKAT, 
                                         function(x)get_coh_dist(x[,1], fil_db = fil_tab_noCH))

#saveRDS(Exome_pc123_srt_SKAT_case_enr_nCH, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/Exome_pc123_SKAT_Enriched_combset2020.rds", compress = T)
saveRDS(Exome_pc123_srt_SKAT_case_enr_nCH, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/Exome_isks_epit_pc1234_SKAT_Enriched_combset2020_ver4.rds", compress = T)
