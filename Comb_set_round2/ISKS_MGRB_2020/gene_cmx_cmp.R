###############SKATbin for complexes identified by PPI approach


library(data.table)
library(VariantAnnotation)
library(stringr)
library(SKAT)
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

`%nin%` = Negate(`%in%`)

##SKAT null function with customised covariate
SKAT_fun_null <- function(x=NULL, pheno_file, p_vec){
  if(is.null(x)){
    obj_N <- SKAT::SKAT_Null_Model(p_vec ~ 1,out_type="D")
  }
  # else if(is.integer(x) || is.numeric()){
  else if(is.null(dim(x))){
    obj_N <- SKAT::SKAT_Null_Model(p_vec ~ x,out_type="D")
  }
  else if(dim(x)[2] > 1){
    nul_for <- as.formula(paste("p_vec", paste(colnames(x), collapse = " + "), sep = " ~ "))
    obj_N <- SKAT::SKAT_Null_Model(nul_for, data = pheno_file, out_type="D")
    #  obj_N <- SKAT_Null_Model(nul_for, data = x, out_type="D")
  }
  return(obj_N)
}

#####make genotype matrix function
make_geno_mat <- function(ftemp_file, maf_thresh, p_vec, samp_id){
  ftemp_tab_var_id <- unique(ftemp_file$VARIANT)
  samp_vec <- list()
  for(m in 1:length(ftemp_tab_var_id)){
    sam_gene_gt <- ftemp_file[ftemp_file$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,82,125:126)]
    sam_gene_gt <- unique(sam_gene_gt)
   # print(dim(sam_gene_gt)[1])
    if(dim(sam_gene_gt)[1] > 1 & length(unique(sam_gene_gt$vep_consequence)) > 1){
      # if(dim(sam_gene_gt)[1] > 1 & length(unique(sam_gene_gt$vep_consequence)) > 1 & length(unique(sam_gene_gt$comb_score)) > 1){
      vep_con <- unique(sam_gene_gt$vep_consequence)
      samp_vec_con <- list()
      for(k in 1:length(vep_con)){
        sam_gene_gt_con <- sam_gene_gt[sam_gene_gt$vep_consequence %in% vep_con[k],]
        sam_gene_gt_con <- unique(sam_gene_gt_con)
        cont_wt <- sum(sam_gene_gt_con[grepl("^[ABZ]",sam_gene_gt_con$SAMPLE),]$comb_score)
        case_wt <- sum(sam_gene_gt_con[!grepl("^[ABZ]",sam_gene_gt_con$SAMPLE),]$comb_score)
        maf_vec_cont <- sum(dim(sam_gene_gt_con[grepl("^[ABZ]",sam_gene_gt_con$SAMPLE),])[1])/(2*length(p_vec))
        maf_vec_case <- sum(dim(sam_gene_gt_con[!grepl("^[ABZ]",sam_gene_gt_con$SAMPLE),])[1])/(2*length(p_vec))
        
        ##genotype matrix  
        sam_gene_gt_con$add_mod <- as.numeric(sam_gene_gt_con$GT)
        sam10 <- ifelse(samp_id %in% sam_gene_gt_con$SAMPLE, 1, 0)
        sam10[which(sam10 != 0)] <- sam_gene_gt_con$add_mod ##additive model
        samp_vec_con[[k]] <- c(unique(sam_gene_gt_con$VARIANT),
                               unique(sam_gene_gt_con$gene_symbol), 
                               unique(sam_gene_gt_con$vep_consequence), 
                               unique(as.character(sam_gene_gt_con$auto_call)),
                               as.numeric(unique(sam_gene_gt_con$comb_score)), 
                               as.numeric(maf_vec_case),
                               as.numeric(maf_vec_cont), 
                               as.numeric(case_wt),
                               as.numeric(cont_wt),
                               sam10)
      }
      
    }
    else{
      ##compute cohort specific MAF
      cont_wt <- sum(sam_gene_gt[grepl("^[ABZ]",sam_gene_gt$SAMPLE),]$comb_score)
      case_wt <- sum(sam_gene_gt[!grepl("^[ABZ]",sam_gene_gt$SAMPLE),]$comb_score)
      maf_vec_cont <- sum(ifelse(is.na(as.numeric(sam_gene_gt$SAMPLE)), 1, 0))/(2*length(p_vec))
      maf_vec_case <- sum(ifelse(!is.na(as.numeric(sam_gene_gt$SAMPLE)) | 
                                   grepl("^CR|^LK", as.character(sam_gene_gt$SAMPLE)), 1, 0))/(2*length(p_vec))
      #maf_vec <- (maf_vec_cont + maf_vec_case)/(2*(1572 + 1110))
      ##genotype matrix  
      sam_gene_gt$add_mod <- as.numeric(sam_gene_gt$GT)
      sam10 <- ifelse(samp_id %in% sam_gene_gt$SAMPLE, 1, 0)
      sam10[which(sam10 != 0)] <- sam_gene_gt$add_mod ##additive model
      ##account for discrepancy in vep_consequence
      #vep_con <- names(sort(table(unlist(strsplit(sam_gene_gt$vep_consequence, split = "&"))),decreasing=TRUE)[1])
      # sam_gene_gt$vep_consequence <- vep_con
      
      
      samp_vec[[m]] <- c(ftemp_tab_var_id[m],
                         unique(sam_gene_gt$gene_symbol), 
                         unique(sam_gene_gt$vep_consequence), 
                         unique(as.character(sam_gene_gt$auto_call)),
                         as.numeric(unique(sam_gene_gt$comb_score)), 
                         as.numeric(maf_vec_case),
                         as.numeric(maf_vec_cont), 
                         as.numeric(case_wt),
                         as.numeric(cont_wt),
                         sam10)
    }
  }
  samp_vec_mat_uni <- do.call("rbind.data.frame", samp_vec)
  colnames(samp_vec_mat_uni) <- c("VARIANT", "gene_symbol", "vep_consequence", "auto_call", "comb_score", "coh_MAF_case",
                                  "coh_MAF_cont", "case_wt", "cont_wt", samp_id)
  if(exists("samp_vec_con")){
    samp_vec_mat_con <- do.call("rbind.data.frame", samp_vec_con)
    colnames(samp_vec_mat_con) <- c("VARIANT", "gene_symbol", "vep_consequence", "auto_call", "comb_score", "coh_MAF_case",
                                    "coh_MAF_cont", "case_wt", "cont_wt", samp_id)
    samp_vec_mat <- rbind.data.frame(samp_vec_mat_uni, samp_vec_mat_con)
  }else{
    samp_vec_mat <- samp_vec_mat_uni
  }
  
  ##Intracohort filter  : equivalent to ~ 5/1661*2 for case ; ~ 5/3209*2 for control
  ##changed to uniform intra_cohort_MAF filter
  ## change to MAF filter to 3/1661*2
  #     samp_vec_mat <- samp_vec_mat[as.numeric(as.character(samp_vec_mat$coh_MAF_case)) <= 0.0015 | 
  #                                    as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) <= 0.001,]
  samp_vec_mat <- samp_vec_mat[!(as.numeric(as.character(samp_vec_mat$coh_MAF_case)) > maf_thresh | 
                                   as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) > maf_thresh),]
  print(dim(samp_vec_mat))
  
  return(samp_vec_mat)
  
}

##function to return table with number of variants in case control
make_table <- function(gene_sym, var_file){
  ftemp_tab <- var_file[var_file$gene_symbol %in% gene_sym,] 
  ftemp_tab <- ftemp_tab[ftemp_tab$VAF >= 0.35 & ftemp_tab$comb_score >= 5.6,]
  ##make table
  
  cs <- c("C3", "C4", "C5")
  cnt_df <- as.data.frame(table(ftemp_tab[ftemp_tab$is_case == 0,]$auto_call))
  if(length(cs[cs %nin% cnt_df$Var1] > 0)){
    cnt_df_add <- cbind.data.frame("Var1" = cs[cs %nin% cnt_df$Var1], "Freq" = 0)
    cnt_df <- rbind.data.frame(cnt_df, cnt_df_add)
    cnt_df <- cnt_df[cnt_df$Var1 %in% cs,]
  }else {cnt_df <- cnt_df }
  case_df <- as.data.frame(table(ftemp_tab[ftemp_tab$is_case == 1,]$auto_call))
  if(length(cs[cs %nin% case_df$Var1] > 0)){
    case_df_add <- cbind.data.frame("Var1" = cs[cs %nin% case_df$Var1], "Freq" = 0)
    case_df <- rbind.data.frame(case_df, case_df_add)
    case_df <- case_df[case_df$Var1 %in% cs,]
  }else {case_df <- case_df }
  mytab <- cbind.data.frame(cnt_df, case_df)
  rownames(mytab) <- mytab$Var1
  mytab <- mytab[,-c(1,3)]
  mytab <- t(mytab)
  rownames(mytab) <- c("control", "case")
  return(mytab)
}

##Burden test

eigen_gene_complex <- function(gene_sym, var_file, pheno_file, p_vec, samp_id){
  
  x <- pheno_file[,c(53,19:22)]
  #null_d <- SKAT_fun_null(x = pheno_file[,c(53,19:22)], p_vec = p_vec)
  nul_for <- as.formula(paste("p_vec", paste(colnames(x), collapse = " + "), sep = " ~ "))
  obj_N <- SKAT::SKAT_Null_Model(nul_for, data = pheno_file, out_type="D")
  #cut_off <- seq(0,10)
  #  score_scale <- c((0.001:50)^2/100,0.5,1000)
#  p_val_cutoff <- list()
#  for(i in 1:length(cut_off))  {
  maf_cut <- 3.5/(2*length(p_vec))
    ftemp_tab <- var_file[var_file$gene_symbol %in% gene_sym,] 
    ftemp_tab <- ftemp_tab[ftemp_tab$VAF >= 0.35 & ftemp_tab$comb_score >= 5.6,]
    if(dim(ftemp_tab)[1] > 1){
      samp_vec_mat <- make_geno_mat(ftemp_tab, maf_cut, p_vec, samp_id)
      ftemp_tab_mod <- samp_vec_mat[,c(1:9)]
      ftemp_tab_mod$comb_score <- as.numeric(as.character(ftemp_tab_mod$comb_score))
      ftemp_tab_mod$case_wt <- as.numeric(as.character(ftemp_tab_mod$case_wt))
      ftemp_tab_mod$cont_wt <- as.numeric(as.character(ftemp_tab_mod$cont_wt))
      ftemp_tab_mod$auto_call <- as.character(ftemp_tab_mod$auto_call)
      gene_mat <- as.matrix(samp_vec_mat[,-c(1:9)])
      class(gene_mat) <- "numeric"
      
      ##SKAT
      p_val_complex <- tryCatch(SKAT::SKATBinary(t(gene_mat), obj_N, method = "Burden", weights = ftemp_tab_mod$comb_score)$p.value, error=function(x) 0.99)
     # p_val_complex <- tryCatch(SKAT::SKATBinary_Robust(t(gene_mat), obj_N, method = "Burden", weights = ftemp_tab_mod$comb_score)$p.value, error=function(x) 0.99)
      wt_pval <-  cbind.data.frame("case_wt" = sum(ftemp_tab_mod$case_wt), 
                                   "cont_wt" = sum(ftemp_tab_mod$cont_wt), p_val_complex)
      }
    else {
      p_val_complex <- 0.99
      wt_pval <-  cbind.data.frame("case_wt" = sum(ftemp_tab_mod$case_wt), 
                                   "cont_wt" = sum(ftemp_tab_mod$cont_wt), p_val_complex)
    }
 
  #return(p_val_complex)
    return(wt_pval)
}



########ISKS
##Input file:
fil_tab_isks <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/all_isksmgrb_combset2020_variants_filt_all_fields_rm_dup.tsv",
                                 sep = "\t", stringsAsFactors = F, header = T)
isks_pheno <- read.table("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_combset_pca.scores_clustered.tsv", 
                                   header = T, sep = "\t", stringsAsFactors = F)
QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv", 
                      header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %in% "TRUE",]
QC2_dat_pass$isFemale <- ifelse(QC2_dat_pass$f_stat < 0.2, 1, 
                                ifelse(QC2_dat_pass$f_stat > 0.8, 0, 2))
Ex_samp_id_isks <- unique(fil_tab_isks$SAMPLE)
isks_pheno <- isks_pheno[isks_pheno$sample %in% Ex_samp_id_isks,]
isks_pheno$gender <- QC2_dat_pass[match(isks_pheno$sample, QC2_dat_pass$new_sampleid), 21]
Ex_samp_id_isks <- Ex_samp_id_isks[match(isks_pheno$sample, Ex_samp_id_isks)]
p_vec <- ifelse(!is.na(as.numeric(as.character(Ex_samp_id_isks))) | grepl("^CR|^LK",as.character(Ex_samp_id_isks)), 1, 0)


eigen_gene_complex("TP53", fil_tab_isks,isks_pheno,p_vec, Ex_samp_id_isks) 
eigen_gene_complex("TINF2", fil_tab_isks,isks_pheno,p_vec, Ex_samp_id_isks)
eigen_gene_complex("TERF1", fil_tab_isks,isks_pheno,p_vec, Ex_samp_id_isks)
eigen_gene_complex("POT1", fil_tab_isks,isks_pheno,p_vec, Ex_samp_id_isks)

##all genes combined
gene_all <- c("TP53", "NF1", "EXT1", "EXT2", "BRCA2", "ERCC2", "SDHA", "SDHB", "SDHD")
eigen_gene_complex(gene_all, fil_tab_isks,isks_pheno,p_vec, Ex_samp_id_isks)
##test on Sheltrin complex
##all genes combined
gene_sheltrin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "ACD")
eigen_gene_complex(gene_sheltrin, fil_tab_isks,isks_pheno,p_vec, Ex_samp_id_isks)


##Extended Sheltrin complex
Sheltrin_comp_extn = c("ACD", "POT1", "TERF1", "TERF2", "TERF2IP", "TINF2", "ATM", 
                       "BAG3", "BLM", "BRCA1", "CALD1", "CLK3", "DCLRE1B", "FANCD2", 
                       "FBXO4", "HSPA4", "KIAA1191", "MRE11A", "NBN", "PINX1", "PRKDC", 
                       "RAD50", "SLX4", "STUB1", "TNKS", "TNKS2", "U2AF2", "UCHL1", 
                       "WRN", "XRCC5", "XRCC6")

eigen_gene_complex(Sheltrin_comp_extn, fil_tab_isks,isks_pheno,p_vec, Ex_samp_id_isks)


##Repeat the above for EPIT and CAIRNS

##EPIT
fil_tab_epit <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/all_epit_mgrb_combset2020_variants_filt_all_fields_nof_brca.tsv", 
           sep = "\t", header = T, stringsAsFactors = F)
Ex_samp_id_epit <- unique(fil_tab_epit$SAMPLE)
p_Data <- read.delim("~/RVAS/Epi_set_2020/pop_PCA/MGRB_EPIT_1000G_combset_pca.scores_clustered.tsv", header = T, sep = "\t", stringsAsFactors = F)
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
##
p_Data_noCH <- p_Data[as.character(p_Data$sample) %in% Ex_samp_id_epit,]
p_Data_noCH$gender <- QC2_epit_mgrb_pass[match(p_Data_noCH$sample, QC2_epit_mgrb_pass$new_sampleid), 3]
##f_brca already removed
#p_Data_noCH$study <- epit_pheno[match(p_Data_noCH$sample, epit_pheno$Sample_Id),18]
#p_Data_noCH$study <- ifelse(is.na(p_Data_noCH$study), "MGRB", p_Data_noCH$study)
##remove 78 Familial breast cases
#`%nin%` = Negate(`%in%`)
#p_Data_noCH <- p_Data_noCH[p_Data_noCH$study %nin% "Familial breast",]

Ex_samp_id_epit <- Ex_samp_id_epit[match(p_Data_noCH$sample, Ex_samp_id_epit)]
p_vec_epit <- ifelse(grepl("^[ABZ]",as.character(p_Data_noCH$sample)), 0, 1)

eigen_gene_complex(gene_all, fil_tab_epit,p_Data_noCH,p_vec_epit,Ex_samp_id_epit)

eigen_gene_complex(gene_sheltrin, fil_tab_epit,p_Data_noCH,p_vec_epit,Ex_samp_id_epit)

eigen_gene_complex(Sheltrin_comp_extn, fil_tab_epit,p_Data_noCH,p_vec_epit,Ex_samp_id_epit)


##Odds ratio vs phenotype

##weighted odd ratio using logistic and firth logistic regression
##logistic

library("MASS")
logistic_gene_complex <- function(gene_sym, var_file, pheno_file, p_vec, samp_id){
  
  #x <- pheno_file[,c(53,19:22)]
  #null_d <- SKAT_fun_null(x = pheno_file[,c(53,19:22)], p_vec = p_vec)
  #nul_for <- as.formula(paste("p_vec", paste(colnames(x), collapse = " + "), sep = " ~ "))
  #obj_N <- SKAT::SKAT_Null_Model(nul_for, data = pheno_file, out_type="D")
  
  maf_cut <- 3.5/(2*length(p_vec))
  ftemp_tab <- var_file[var_file$gene_symbol %in% gene_sym,] 
  ftemp_tab <- ftemp_tab[ftemp_tab$VAF >= 0.35 & ftemp_tab$comb_score >= 5.6,]
  if(dim(ftemp_tab)[1] > 1){
    samp_vec_mat <- make_geno_mat(ftemp_tab, maf_cut, p_vec, samp_id)
    ftemp_tab_mod <- samp_vec_mat[,c(1:9)]
    ftemp_tab_mod$comb_score <- as.numeric(as.character(ftemp_tab_mod$comb_score))
    ftemp_tab_mod$case_wt <- as.numeric(as.character(ftemp_tab_mod$case_wt))
    ftemp_tab_mod$cont_wt <- as.numeric(as.character(ftemp_tab_mod$cont_wt))
    ftemp_tab_mod$auto_call <- as.character(ftemp_tab_mod$auto_call)
    gene_mat <- as.matrix(samp_vec_mat[,-c(1:9)])
    class(gene_mat) <- "numeric"
    
    ##logistic regression
    Weighted_mat <- t(gene_mat) * ftemp_tab_mod$comb_score
    Weighted_mat_rv_sum <- apply(Weighted_mat, 1, sum)
    Weighted_mat_df <- cbind.data.frame("p_vec" = p_vec, "weight" = Weighted_mat_rv_sum, pheno_file[,c(53,19:22)])
    Weighted_mat_df$p_vec <- p_vec
    output <- summary(glm(p_vec ~ weight + gender + PC1 + PC2 + PC3 + PC4, data=Weighted_mat_df, family=binomial(link="logit")))
    ##Adjust for PC's and gender
    OR <- exp(output$coefficients[2,1])
    z_score <- output$coefficients[2,3]
    p_value <- output$coefficients[2,4]
  #  p_val_complex <- tryCatch(SKAT::SKATBinary(t(gene_mat), obj_N, method = "Burden", weights = ftemp_tab_mod$comb_score)$p.value, error=function(x) 0.99)
    wt_pval <-  cbind.data.frame("case_wt" = sum(ftemp_tab_mod$case_wt), 
                                 "cont_wt" = sum(ftemp_tab_mod$cont_wt), OR, z_score, p_value)
  }
  else {
  #  cont_wt = as.numeric(ftemp_tab[grepl("^[ABZ]",as.character(ftemp_tab$SAMPLE)),]$comb_score)
  #  case_wt = as.numeric(ftemp_tab[!grepl("^[ABZ]",as.character(ftemp_tab$SAMPLE)),]$comb_score)
    OR <- 1
    z_score <- 0
    p_value <- 0.99
    wt_pval <-  cbind.data.frame(OR, z_score, p_value)
  }
  
  #return(p_val_complex)
  return(wt_pval)
}

##ISKS

logistic_gene_complex("TP53", fil_tab_isks,isks_pheno,p_vec,Ex_samp_id_isks)

logistic_gene_complex(gene_all, fil_tab_isks,isks_pheno,p_vec,Ex_samp_id_isks)

logistic_gene_complex(gene_sheltrin, fil_tab_isks,isks_pheno,p_vec,Ex_samp_id_isks)

logistic_gene_complex(Sheltrin_comp_extn, fil_tab_isks,isks_pheno,p_vec,Ex_samp_id_isks)

##EPIT

logistic_gene_complex("TP53", fil_tab_epit,p_Data_noCH,p_vec_epit,Ex_samp_id_epit)

logistic_gene_complex(gene_all, fil_tab_epit,p_Data_noCH,p_vec_epit,Ex_samp_id_epit)

logistic_gene_complex(gene_sheltrin, fil_tab_epit,p_Data_noCH,p_vec_epit,Ex_samp_id_epit)

logistic_gene_complex(Sheltrin_comp_extn, fil_tab_epit,p_Data_noCH,p_vec_epit,Ex_samp_id_epit)


##firth logistic 
library(logistf)

logistic_firth_gene_complex <- function(gene_sym, var_file, pheno_file, p_vec, samp_id){
  
  #x <- pheno_file[,c(53,19:22)]
  #null_d <- SKAT_fun_null(x = pheno_file[,c(53,19:22)], p_vec = p_vec)
  #nul_for <- as.formula(paste("p_vec", paste(colnames(x), collapse = " + "), sep = " ~ "))
  #obj_N <- SKAT::SKAT_Null_Model(nul_for, data = pheno_file, out_type="D")
  
  maf_cut <- 3.5/(2*length(p_vec))
  ftemp_tab <- var_file[var_file$gene_symbol %in% gene_sym,] 
  ftemp_tab <- ftemp_tab[ftemp_tab$VAF >= 0.35 & ftemp_tab$comb_score >= 5.6,]
  if(dim(ftemp_tab)[1] > 1){
    samp_vec_mat <- make_geno_mat(ftemp_tab, maf_cut, p_vec, samp_id)
    ftemp_tab_mod <- samp_vec_mat[,c(1:9)]
    ftemp_tab_mod$comb_score <- as.numeric(as.character(ftemp_tab_mod$comb_score))
    ftemp_tab_mod$case_wt <- as.numeric(as.character(ftemp_tab_mod$case_wt))
    ftemp_tab_mod$cont_wt <- as.numeric(as.character(ftemp_tab_mod$cont_wt))
    ftemp_tab_mod$auto_call <- as.character(ftemp_tab_mod$auto_call)
    gene_mat <- as.matrix(samp_vec_mat[,-c(1:9)])
    class(gene_mat) <- "numeric"
    
    ##logistic regression
    Weighted_mat <- t(gene_mat) * ftemp_tab_mod$comb_score
    Weighted_mat_rv_sum <- apply(Weighted_mat, 1, sum)
    Weighted_mat_df <- cbind.data.frame("p_vec" = p_vec, "weight" = Weighted_mat_rv_sum, pheno_file[,c(53,19:22)])
    Weighted_mat_df$p_vec <- p_vec
    output <- logistf(p_vec ~ weight + gender + PC1 + PC2 + PC3 + PC4, data=Weighted_mat_df)
    ##Adjust for PC's and gender
    OR <- exp(output$coefficients[2])
   # z_score <- output$coefficients[23]
    p_value <- output$prob[2]
    #  p_val_complex <- tryCatch(SKAT::SKATBinary(t(gene_mat), obj_N, method = "Burden", weights = ftemp_tab_mod$comb_score)$p.value, error=function(x) 0.99)
    wt_pval <-  cbind.data.frame("case_wt" = sum(ftemp_tab_mod$case_wt), 
                                 "cont_wt" = sum(ftemp_tab_mod$cont_wt), OR, p_value)
  }
  else {
    #  cont_wt = as.numeric(ftemp_tab[grepl("^[ABZ]",as.character(ftemp_tab$SAMPLE)),]$comb_score)
    #  case_wt = as.numeric(ftemp_tab[!grepl("^[ABZ]",as.character(ftemp_tab$SAMPLE)),]$comb_score)
    OR <- 1
 #   z_score <- 0
    p_value <- 0.99
    wt_pval <-  cbind.data.frame(OR, p_value)
  }
  
  #return(p_val_complex)
  return(wt_pval)
}

logistic_firth_gene_complex("TP53", fil_tab_isks,isks_pheno,p_vec,Ex_samp_id_isks)
logistic_firth_gene_complex("TINF2", fil_tab_isks,isks_pheno,p_vec,Ex_samp_id_isks)
logistic_firth_gene_complex(gene_all, fil_tab_isks,isks_pheno,p_vec,Ex_samp_id_isks)
logistic_firth_gene_complex(gene_sheltrin, fil_tab_isks,isks_pheno,p_vec,Ex_samp_id_isks)
logistic_firth_gene_complex(Sheltrin_comp_extn, fil_tab_isks,isks_pheno,p_vec,Ex_samp_id_isks)

##Cairns data
fil_tab_cairns <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CAIRNS/all_cairns_mgrb_combset2020_variants_filt_all_fields_rmdup_fin.tsv", 
                             sep = "\t",header = T, stringsAsFactors = F)
Ex_samp_id_cairns <- unique(fil_tab_cairns$SAMPLE)
QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv",
                      header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %in% "TRUE",]
QC2_dat_pass$isFemale <- ifelse(QC2_dat_pass$f_stat < 0.2, 1,
                                ifelse(QC2_dat_pass$f_stat > 0.8, 0, 2))
QC2_dat_pass_mgrb <- QC2_dat_pass[grepl("^[ABZ]", QC2_dat_pass$new_sampleid),]
QC2_dat_pass_mgrb <- QC2_dat_pass_mgrb[,c(1,21)]

p_Data_cairn <- read.delim("~/RVAS/Cairns_set/pop_PCA/MGRB_CAIRNS_1000G_combset_pca.scores_clustered.tsv", header = T, sep = "\t", stringsAsFactors = F)
p_Data_noCH_cairn <- p_Data_cairn[as.character(p_Data_cairn$sample) %in% Ex_samp_id_cairns,]
p_Data_noCH_cairn$gender <- QC2_dat_pass_mgrb[match(p_Data_noCH_cairn$sample, QC2_dat_pass_mgrb$new_sampleid), 2]
##Add Cairns phenotype
Cairns_pheno <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/CAIRNS_MGRB_2020/cairns_ages_gender_telolength.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
Cairns_pheno$Gender <- ifelse(Cairns_pheno$Gender %in% "F", 1, 0)
p_Data_noCH_cairn$gender <- ifelse(is.na(p_Data_noCH_cairn$gender), Cairns_pheno[match(p_Data_noCH_cairn$sample, Cairns_pheno$Sample), 4],p_Data_noCH_cairn$gender)
Ex_samp_id_cairns <- Ex_samp_id_cairns[match(p_Data_noCH_cairn$sample, Ex_samp_id_cairns)]

#binary phenotype vector 
p_vec_cairn <- ifelse(grepl("^[ABZ]",as.character(p_Data_noCH_cairn$sample)), 0, 1)

##ISKS
cpx_list <- list(gene_all,gene_sheltrin, Sheltrin_comp_extn)
cpx_list_OR_ISKS <- list()
cpx_list_OR_EPIT <- list()
cpx_list_OR_CAIRNS <- list()
for(i in 1:3){
  cpx_list_OR_ISKS[[i]] <- logistic_firth_gene_complex(cpx_list[[i]], fil_tab_isks,isks_pheno,p_vec,Ex_samp_id_isks)
  cpx_list_OR_EPIT[[i]] <- logistic_firth_gene_complex(cpx_list[[i]], fil_tab_epit,p_Data_noCH,p_vec_epit,Ex_samp_id_epit)
  cpx_list_OR_CAIRNS[[i]] <- logistic_firth_gene_complex(cpx_list[[i]], fil_tab_cairns,p_Data_noCH_cairn,p_vec_cairn,Ex_samp_id_cairns)
}
cpx_ISKS_OR_df <- do.call("rbind.data.frame", cpx_list_OR_ISKS)
rownames(cpx_ISKS_OR_df) <- c("pos_cont", "Sheltrin", "Sheltrin_extn")
cpx_EPIT_OR_df <- do.call("rbind.data.frame", cpx_list_OR_EPIT)
rownames(cpx_EPIT_OR_df) <- c("pos_cont", "Sheltrin", "Sheltrin_extn")
cpx_CAIRNS_OR_df <- do.call("rbind.data.frame", cpx_list_OR_CAIRNS)
rownames(cpx_CAIRNS_OR_df) <- c("pos_cont", "Sheltrin", "Sheltrin_extn")


##################################Acceptable for paper
##gene_complexes
##all genes combined:positive control
gene_all <- c("TP53", "NF1", "EXT1", "EXT2", "BRCA2", "ERCC2", "SDHA", "SDHB", "SDHD")
gene_all <- c("TP53", "NF1")
##Sheltrin complex
gene_sheltrin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "ACD")
##Extended Sheltrin complex
Sheltrin_comp_extn = c("ACD", "POT1", "TERF1", "TERF2", "TERF2IP", "TINF2", "ATM", 
                       "BAG3", "BLM", "BRCA1", "CALD1", "CLK3", "DCLRE1B", "FANCD2", 
                       "FBXO4", "HSPA4", "KIAA1191", "MRE11A", "NBN", "PINX1", "PRKDC", 
                       "RAD50", "SLX4", "STUB1", "TNKS", "TNKS2", "U2AF2", "UCHL1", 
                       "WRN", "XRCC5", "XRCC6")

##Odds ratio using fisher's exact test(use ppi_final_res... file for fisher's test)
Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Exome_skat_para_result_isks_combset2020_uni_MAF_PC1234_ver4_rect_test.rds")

df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})
Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Exome_para_pc123_SKAT_Enriched_ISKS_2020.rds")
##Add p-value
Exome_pc123_srt_SKAT_case_enr_nCH_pval <- Map(cbind.data.frame, Exome_pc123_srt_SKAT_case_enr_nCH, Exome_pc123_srt_SKAT)

isks_mgrb_genes <- Exome_pc123_srt_SKAT_case_enr_nCH_pval[[4]]

##Epithelial data
Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/round2/Exome_skat_para_result_EPIT_combset2020_uni_MAF_PC1234_ver4_no_fbrca.rds")

df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})
Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/round2/Exome_para_pc123_SKAT_Enriched_EPIT_2020_no_fbrca.rds")
##Add p-value
Exome_pc123_srt_SKAT_case_enr_nCH_pval <- Map(cbind.data.frame, Exome_pc123_srt_SKAT_case_enr_nCH, Exome_pc123_srt_SKAT)

epit_mgrb_genes <- Exome_pc123_srt_SKAT_case_enr_nCH_pval[[4]]

##Cairns data
Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CAIRNS/Exome_skat_para_result_CAIRNS_combset2020_uni_MAF_PC1234_ver4.rds")

df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})
Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CAIRNS/Exome_para_pc123_SKAT_Enriched_CAIRNS_2020.rds")
##Add p-value
Exome_pc123_srt_SKAT_case_enr_nCH_pval <- Map(cbind.data.frame, Exome_pc123_srt_SKAT_case_enr_nCH, Exome_pc123_srt_SKAT)

cairns_mgrb_genes <- Exome_pc123_srt_SKAT_case_enr_nCH_pval[[4]]

########

isks_epit_cairn <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_EPIT_cmp/round2/isks_ppi_epit_no_fbrca_cairns_comb_genes_top.tsv", 
                              sep = "\t", header = T, stringsAsFactors = F)
##use only enriched(SKATBin <0.1) genes for sheltrin extn complex, for sarc genes use TP53, NF1.
Sheltrin_comp_extn_enr <- isks_epit_cairn[isks_epit_cairn$gene %in% Sheltrin_comp_extn & isks_epit_cairn$pval_SKATbin < 0.1,]$gene

##Exact test OR and Pvalues
gene_all <- c("TP53", "NF1")
cpx_list <- list(gene_all,gene_sheltrin, Sheltrin_comp_extn, Sheltrin_comp_extn_enr)
names(cpx_list) <- c("Sarc_genes", "Sheltrin", "Sheltrin_extn", "Sheltrin_enr")
cpx_OR_fisher <- function(ppi_res,case_coh_size, cont_coh_size, coh){
  ft_df <- list()
for(i in 1:length(cpx_list)){
  ppi_res_tab <- ppi_res[ppi_res$gene %in% cpx_list[[i]],]
  inp <- c(sum(ppi_res_tab[,1]), case_coh_size - sum(ppi_res_tab[,1]) , 
           sum(ppi_res_tab[,2]), cont_coh_size - sum(ppi_res_tab[,2]))
  sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
  colnames(sim_mat) <- c("case", "cont")
  rownames(sim_mat) <- c("hits", "no_hits")
  #ft <- fisher.test(sim_mat, alternative = "greater")
  ft <- fisher.test(sim_mat)
  ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.99)
  ft_df[[i]] <- cbind.data.frame("gene" = names(cpx_list)[i] ,"Cases" = sum(ppi_res_tab[,1]),
                                 "Controls" = sum(ppi_res_tab[,2]),
                                 "Fish_pval" = ft$p.value,"CI_lower" = ft$conf.int[1],
                                 "CI_upper" = ft$conf.int[2],
                            "OR_Fish" = ft$estimate, "Coh" = coh)
}
  return(ft_df)
}

cpx_ISKS_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_mgrb_genes, 1646, 3205, "ISKSvsMGRB"))
cpx_EPIT_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(epit_mgrb_genes, 842, 3205, "EPITvsMGRB"))
cpx_CAIRN_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(cairns_mgrb_genes, 413, 3205, "CAIRNvsMGRB"))
df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_OR_df, cpx_EPIT_OR_df, cpx_CAIRN_OR_df)
df_comb_rnd2_filt <- df_comb_rnd2[df_comb_rnd2$Fish_pval < 0.1, ]

##make_tab for ISKS_MGRB, EPIT_MGRB and CAIRNS_MGRB without C3.
##CAIRNS
Cairns_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CAIRNS/Exome_para_pc123_SKAT_Enriched_CAIRNS_2020_minusC3.rds")
Cairns_noC3_fin <- Cairns_noC3[[4]]
isks_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Exome_para_pc123_SKAT_Enriched_ISKS_2020_minusC3.rds")
isks_noC3_fin <- isks_noC3[[4]]
epit_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/round2/Exome_para_pc123_SKAT_Enriched_EPIT_2020_no_fbrca_minusC3.rds")
epit_noC3_fin <- epit_noC3[[4]]
##no_C3
cpx_ISKS_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_noC3_fin, 1646, 3205, "ISKSvsMGRB"))
cpx_EPIT_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(epit_noC3_fin, 842, 3205, "EPITvsMGRB"))
cpx_CAIRN_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(Cairns_noC3_fin, 413, 3205, "CAIRNvsMGRB"))
noC3_df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_OR_df_noC3, cpx_EPIT_OR_df_noC3, cpx_CAIRN_OR_df_noC3)

noC3_df_comb_rnd2$OR_Fish <- ifelse(noC3_df_comb_rnd2$OR_Fish == 0, 1, noC3_df_comb_rnd2$OR_Fish)
noC3_df_comb_rnd2_filt <- noC3_df_comb_rnd2[noC3_df_comb_rnd2$Fish_pval < 0.1,]
##Forest plot
#define colours for dots and bars
library(ggplot2)
forest_custplot <- function(df){
dotCOLS = c("#a6d8f0","#f9b282", "#78f542")
barCOLS = c("#008fd5","#de6b35", "#7a9406")

p <- ggplot(df, aes(x=gene, y=OR_Fish, ymin=CI_lower, ymax=CI_upper,col=Coh,fill=Coh)) + 
  #specify position here
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  #specify position here too
  geom_point(size=5, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS)+
  scale_x_discrete(name="Protein modules") +
  scale_y_continuous(name="Odds ratio") +
  coord_flip() +
  theme_minimal()
return(p)
}

`%nin%` = Negate(`%in%`)
#t1_df <- noC3_df_comb_rnd2[noC3_df_comb_rnd2$gene %nin% c("Sheltrin"),]
#forest_custplot(t1_df)
forest_custplot(noC3_df_comb_rnd2)
forest_custplot(noC3_df_comb_rnd2_filt)
forest_custplot(df_comb_rnd2)
forest_custplot(df_comb_rnd2_filt)

##ISKS_complex, ISKS_simple, ISKS_TAS; use fil_tab, PID file (genomic class)

