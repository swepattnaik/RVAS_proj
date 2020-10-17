##logistic and firth regression test
library(data.table)
library(VariantAnnotation)
library(stringr)
#library(SKAT)
library(MASS)
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

#####make genotype matrix function
make_geno_mat <- function(ftemp_file, maf_thresh, p_vec, samp_id){
  ftemp_tab_var_id <- unique(ftemp_file$VARIANT)
  samp_vec <- list()
  for(m in 1:length(ftemp_tab_var_id)){
    sam_gene_gt <- ftemp_file[ftemp_file$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,82,127:128)]
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


##Input files : genotype matrix and covariates
fil_tab_noCH <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
`%nin%` = Negate(`%in%`)

Ex_samp_id <- unique(fil_tab_noCH$SAMPLE)
##QC data
QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv", 
                      header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %in% "TRUE",]
QC2_dat_pass$isFemale <- ifelse(QC2_dat_pass$f_stat < 0.2, 1, 
                                ifelse(QC2_dat_pass$f_stat > 0.8, 0, 2))

p_Data <- read.table("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_combset_pca.scores_clustered.tsv", 
                     header = T, sep = "\t", stringsAsFactors = F)
p_Data <- p_Data[p_Data$superPopulation %in% c("ISKS", "RISC", "LIONS", "MGRB"),]
p_Data_MGRB <-  p_Data[p_Data$superPopulation %in% c("MGRB"),]
#remove duplicates
dup_samp <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/fin_samp_dup_drop.txt", sep = "", header = T,
                       stringsAsFactors = F)
dup_samp <- dup_samp$x[grepl("^[ABZ]", dup_samp$x)]
p_Data_MGRB <- p_Data_MGRB[p_Data_MGRB$sample %nin% dup_samp,]
p_Data_MGRB$rect_sam <- p_Data_MGRB$sample
p_Data_MGRB <- p_Data_MGRB[p_Data_MGRB$rect_sam %in% QC2_dat_pass$new_sampleid,]
##rename p_Data
p_Data_ISKS <-  p_Data[p_Data$superPopulation %in% c("ISKS", "RISC", "LIONS"),]
rect_sam_dat <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/ISKS_RISC_LIONS_final_freeze.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
p_Data_ISKS$rect_sam <- rect_sam_dat[match(p_Data_ISKS$sample,rect_sam_dat$JCInputID),9]
p_Data_ISKS <- p_Data_ISKS[!is.na(p_Data_ISKS$rect_sam),]
##two samples missed due to mislabelling and QC2
#rect_sam_dat$JCInputRecID[rect_sam_dat$JCInputRecID %nin% p_Data_ISKS$sample]
p_Data_ISKS <- p_Data_ISKS[p_Data_ISKS$rect_sam %in% QC2_dat_pass$new_sampleid,] ##3105 is lost(QC2 fail)
p_Data_noCH <- rbind.data.frame(p_Data_ISKS, p_Data_MGRB)
p_Data_noCH <- p_Data_noCH[as.character(p_Data_noCH$rect_sam) %in% Ex_samp_id,]

Ex_samp_id <- Ex_samp_id[match(p_Data_noCH$rect_sam, Ex_samp_id)]
##filter out QC fail cases
fil_tab_noCH <- fil_tab_noCH[fil_tab_noCH$SAMPLE %in% Ex_samp_id,]

##Add gender information
p_Data_noCH$gender <- QC2_dat_pass[match(p_Data_noCH$rect_sam, QC2_dat_pass$new_sampleid), 21]

#binary phenotype vector 
p_vec <- ifelse(!is.na(as.numeric(as.character(p_Data_noCH$rect_sam))) | grepl("^CR|^LK",as.character(p_Data_noCH$rect_sam)), 1, 0)

logistic_gene_complex <- function(gene_sym, var_file, pheno_file, p_vec, samp_id){
  
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
    Weighted_mat <- t(gene_mat) * ftemp_tab_mod$comb_score * 10
    Weighted_mat_rv_sum <- apply(Weighted_mat, 1, sum)
    Weighted_mat_df <- cbind.data.frame("p_vec" = p_vec, "weight" = Weighted_mat_rv_sum, pheno_file[,c(54,19:22)])
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

logistic_gene_complex("TP53", fil_tab_noCH,p_Data_noCH,p_vec,Ex_samp_id)

##firth logistic 
library(logistf)

logistic_firth_gene_complex <- function(gene_sym, var_file, pheno_file, p_vec, samp_id){
  
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
    Weighted_mat_df <- cbind.data.frame("p_vec" = p_vec, "weight" = Weighted_mat_rv_sum, pheno_file[,c(54,19:22)])
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
       OR <- 1
    #   z_score <- 0
    p_value <- 0.99
    wt_pval <-  cbind.data.frame(OR, p_value)
  }
  
  #return(p_val_complex)
  return(wt_pval)
}

logistic_firth_gene_complex("TP53", fil_tab_noCH,p_Data_noCH,p_vec,Ex_samp_id)
logistic_firth_gene_complex("TINF2", fil_tab_noCH,p_Data_noCH,p_vec,Ex_samp_id)
