##Note that this analysis include singletons
##changes included on Aug-06-2019:
#added cohort_MAF computation feature; 
#input is better filtered compared to Aug4 input
##Aug12: SKAT made feasible for combined ISKS and RISC samples
#the rare variants are filtered using intra cohort_MAF filter that removes variants that are found in > 3 samples
#in either case of control(added after sep05 output)
#ISKS vs EPIT+MGRB

#.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", "/home/shu/R/x86_64-redhat-linux-gnu-library/3.5/", .libPaths() ) )

##Using D data structure

library(data.table)
library(VariantAnnotation)
library(stringr)
library(SKAT)

###get input files and process
##Experiments

#fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/all_isksmgrb_combset2020_variants_filt_all.tsv",
#                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab_isks <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/all_isksmgrb_combset2020_variants_filt_all_fields.tsv",
                           sep = "\t", stringsAsFactors = F, header = T)
##fil_tab_isks has risk_factor assigned as C5
fil_tab_isks <- fil_tab_isks[!is.na(fil_tab_isks$SAMPLE),]
#fil_tab_isks_set <- fil_tab_isks[fil_tab_isks$is_case == 1,]
#dim(fil_tab_isks_set)

fil_tab_epit <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/all_epit_mgrb_combset2020_variants_filt_all_fields_risk_factor.tsv",
                           sep = "\t", stringsAsFactors = F, header = T)
fil_tab_epit <- fil_tab_epit[!is.na(fil_tab_epit$SAMPLE),]
fil_tab_epit_set <- fil_tab_epit[fil_tab_epit$is_case == 1,]
dim(fil_tab_epit_set)

fil_tab <- rbind.data.frame(fil_tab_isks, fil_tab_epit_set)

##Additive model
Ex_samp_id <- unique(fil_tab$SAMPLE)
#remove duplicates
dup_samp <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/fin_samp_dup_drop.txt", sep = "", header = T,
                       stringsAsFactors = F)
`%nin%` = Negate(`%in%`)
Ex_samp_id <- Ex_samp_id[Ex_samp_id %nin% dup_samp$x]
#Ex_var_id <- unique(fil_tab$VARIANT)

######get phenotype data to control for age and sex and PC's; 

QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv", 
                      header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %in% "TRUE",]
QC2_dat_pass$isFemale <- ifelse(QC2_dat_pass$f_stat < 0.2, 1, 
                                ifelse(QC2_dat_pass$f_stat > 0.8, 0, 2))
QC2_dat_pass <- QC2_dat_pass[,c(1,20:21)]
#QC2_dat_pass_isks <- QC2_dat_pass[!grepl("^[ABZ]", QC2_dat_pass$new_sampleid),]
#QC2_dat_pass_isks <- QC2_dat_pass_isks[,c(1,20:21)]

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

QC2_epit_isks_pass <- rbind.data.frame(QC2_epit_pass, QC2_dat_pass)

#table(ifelse(!is.na(as.numeric(as.character(QC2_epit_isks_pass$new_sampleid))) | 
#grepl("^CR|^LK", QC2_epit_isks_pass$new_sampleid), 1, 0))
##note p_Data_PC_comb does not change
#p_Data <- read.table("~/RVAS/ISKS_MGRB_gender_age_3665.tsv", header = T, sep = "\t")
p_Data <- read.delim("~/RVAS/comb_set_2020/epit_isks_mgrb_pop_PCA/MGRB_ISKS_EPIT_1000G_combset_pca.scores_clustered.tsv", 
                     header = T, sep = "\t", stringsAsFactors = F)
#p_Data$sample <- gsub("http://purl.org/net/grafli/collectedsample#", "", p_Data$sample)

#p_Data$sample <- gsub("isks", "", p_Data$sample)
#p_Data$sample <- gsub("risc", "", p_Data$sample)
p_Data_noCH <- p_Data[as.character(p_Data$sample) %in% Ex_samp_id,]
##retain QC pass samples## 66 EPIT samples failed QC
p_Data_noCH <- p_Data_noCH[p_Data_noCH$sample %in% QC2_epit_isks_pass$new_sampleid,]
##the 7 cases lost QC2_epit_isks_pass$new_sampleid[!(QC2_epit_isks_pass$new_sampleid %in% p_Data_noCH$sample)]
## are due to the DT_exome_clean step where the variants did not make it into the cleaned list perhaps due
## to VAF or comb_score filters
#[1] "02591d74-c7c0-4074-96df-5d5f171485fb" "668b5a09-a705-4074-a9e1-62010639a10a"
#[3] "68597a0f-67f3-4d7a-97fb-2afe6539e67f" "9ff2ebef-54af-4f40-8a3e-324dd4ab860b"
#[5] "ce1ce93c-3816-4e35-8e92-0edfc1a499d0" "dc3041a0-63cc-40ac-91d5-b956ab80144d"
#[7] "3083" 
##drop MGRB sample BAAUD; not in latest call
#samp_ID_match <- samp_ID_match[grep("BAAUD", samp_ID_match, invert = T)]
#p_Data_noCH <- p_Data_noCH[match(samp_ID_match, p_Data_noCH$sample),]
Ex_samp_id <- Ex_samp_id[match(p_Data_noCH$sample, Ex_samp_id)]
##filter out QC fail cases
fil_tab <- fil_tab[fil_tab$SAMPLE %in% Ex_samp_id,]
##remove from gene matrix
#col_rm <- which(colnames(D_tab) ==  "BAAUD")
#D_tab <- D_tab[,-c(col_rm)]
##Add gender information

p_Data_noCH$gender <- QC2_epit_isks_pass[match(p_Data_noCH$sample, QC2_epit_isks_pass$new_sampleid), 3]

##Use prob.NFE for future QC (remove samples with low or FALSE prob.NFE)

#binary phenotype vector 
p_vec <- ifelse(!is.na(as.numeric(as.character(p_Data_noCH$sample))) | grepl("^CR|^LK",as.character(p_Data_noCH$sample)), 1, 0)


##SKAT null function with customised covariate
SKAT_fun_null <- function(x=NULL, p_vec){
  if(is.null(x))
    obj_N <- SKAT_Null_Model(p_vec ~ 1,out_type="D")
  # else if(is.integer(x) || is.numeric()){
  else if(is.null(dim(x))){
    obj_N <- SKAT_Null_Model(p_vec ~ x,out_type="D")
  }
  else if(dim(x)[2] > 1){
    nul_for <- as.formula(paste("p_vec", paste(colnames(x), collapse = " + "), sep = " ~ "))
    obj_N <- SKAT_Null_Model(nul_for, data = p_Data_noCH, out_type="D")
  }
  return(obj_N)
}
##SKAT function for SKAT, SKATBinary, SKATO, SKAT_ERA
SKAT_run <- function(geno_mat, gene, x=NULL, p_vec, cust_weight = NULL, rho = NULL){
  
  if(dim(geno_mat)[2] > 1 ){
    null_d <- SKAT_fun_null(x, p_vec)
    pval_SKAT <- SKAT(t(geno_mat), null_d, weights = cust_weight, r.corr = rho)$p.value
    pval_SKATbin <- SKATBinary(t(geno_mat), null_d, method = "Burden", weights = cust_weight)$p.value
    ##ER: Effective resampling for best p-value calculation
    # pval_SKATbin <- SKATBinary(t(geno_mat), null_d, method.bin = "ER", weights = cust_weight, r.corr = rho)$p.value
    #  pval_SKATO <- SKAT(t(geno_mat),null_d, method='optimal.adj', weights = cust_weight)$p.value ##r.corr will be ignored
    pval_SKATO <- SKATBinary(t(geno_mat),null_d, method='SKATO', weights = cust_weight)$p.value
    pval_SKATera <- SKATBinary(t(geno_mat),null_d, method.bin="ER.A", weights = cust_weight)$p.value
    skat_pv <- cbind.data.frame("eg_ID" = gene, pval_SKAT, pval_SKATbin, pval_SKATO, pval_SKATera)
    return(skat_pv)
  }
  else if(dim(geno_mat)[2] == 1){
    null_d <- SKAT_fun_null(x, p_vec)
    pval_SKAT <- SKAT(geno_mat, null_d, weights = cust_weight, r.corr = rho)$p.value
    pval_SKATbin <- SKATBinary(t(geno_mat), null_d, method = "Burden", weights = cust_weight)$p.value
    ##ER: Effective resampling for best p-value calculation
    # pval_SKATbin <- SKATBinary(t(geno_mat), null_d, method.bin = "ER", weights = cust_weight, r.corr = rho)$p.value
    #pval_SKATO <- SKAT(geno_mat,null_d, method='optimal.adj', weights = cust_weight)$p.value ##r.corr will be ignored
    pval_SKATO <- SKATBinary(geno_mat,null_d, method='SKATO', weights = cust_weight)$p.value
    pval_SKATera <- SKATBinary(geno_mat,null_d, method.bin="ER.A", weights = cust_weight)$p.value
    skat_pv <- cbind.data.frame("eg_ID" = gene, pval_SKAT, pval_SKATbin, pval_SKATO, pval_SKATera)
    
    return(skat_pv)
  }
}

##Run SKAT on variants collapsed by gene

DT_skat_snv_str_pc123 <- list()
#genes <- unique(fil_tab$gene_symbol)

skat_pvals <- list()
skat_pvals_pc12 <- list()
skat_pvals_sex <- list()
skat_pvals_sex_pc12 <- list()
#    skat_pvals_age <- list()
#   skat_pvals_age_pc12 <- list()
#   skat_pvals_age_sex <- list()
#    skat_pvals_age_sex_pc12 <- list()

#for(i in 1:length(genes)){
# print(i)
#  print(genes[i])
##process genes in a loop

##parallel process genes
para_SKAT <- function(genes){
  ftemp_tab <- fil_tab[fil_tab$gene_symbol %in% genes,]
  ftemp_tab <- ftemp_tab[ftemp_tab$comb_score >= 5.6 & as.numeric(ftemp_tab$VAF) >= 0.35, ]
  
  print(max(ftemp_tab$comb_score))
  if(dim(ftemp_tab) == 0 ){
    next
  }
  else{
    ftemp_tab_var_id <- unique(ftemp_tab$VARIANT)
    samp_vec <- list()
    for(m in 1:length(ftemp_tab_var_id)){
      sam_gene_gt <- ftemp_tab[ftemp_tab$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,82,125:126)]
      sam_gene_gt <- unique(sam_gene_gt)
      if(dim(sam_gene_gt)[1] > 1 & length(unique(sam_gene_gt$vep_consequence)) > 1){
        # if(dim(sam_gene_gt)[1] > 1 & length(unique(sam_gene_gt$vep_consequence)) > 1 & length(unique(sam_gene_gt$comb_score)) > 1){
        vep_con <- unique(sam_gene_gt$vep_consequence)
        samp_vec_con <- list()
        for(k in 1:length(vep_con)){
          sam_gene_gt_con <- sam_gene_gt[sam_gene_gt$vep_consequence %in% vep_con[k],]
          sam_gene_gt_con <- unique(sam_gene_gt_con)
          maf_vec_cont <- sum(ifelse(is.na(as.numeric(as.character(sam_gene_gt_con$sample))) &
                                       !grepl("^CR|^LK",as.character(sam_gene_gt_con$sample)), 1, 0))/(2*length(p_vec))
          maf_vec_case <- sum(ifelse(!is.na(as.numeric(as.character(sam_gene_gt_con$sample))) & 
                                       grepl("^CR|^LK",as.character(sam_gene_gt_con$sample)), 1, 0))/(2*length(p_vec))
          ##genotype matrix  
          sam_gene_gt_con$add_mod <- as.numeric(sam_gene_gt_con$GT)
          sam10 <- ifelse(Ex_samp_id %in% sam_gene_gt_con$SAMPLE, 1, 0)
          sam10[which(sam10 != 0)] <- sam_gene_gt_con$add_mod ##additive model
          samp_vec_con[[k]] <- c(unique(sam_gene_gt_con$VARIANT),
                                 unique(sam_gene_gt_con$gene_symbol), 
                                 unique(sam_gene_gt_con$vep_consequence), 
                                 unique(as.character(sam_gene_gt_con$auto_call)),
                                 as.numeric(unique(sam_gene_gt_con$comb_score)), as.numeric(maf_vec_case),
                                 as.numeric(maf_vec_cont), sam10)
        }
        #  vep_cod <- sapply(str_extract_all(sam_gene_gt$vep_Codons, "[A-Z]+"), paste, collapse= '_')
        #  ref_alt <- unlist(lapply(strsplit(sam_gene_gt$VARIANT, split = ":"), function(x) paste(x[3], x[4], sep = "_")))
        #  sam_gene_gt <- sam_gene_gt[which(vep_cod %in% ref_alt),]
        
      }
      else{
        ##compute cohort specific MAF
        maf_vec_cont <- sum(ifelse(is.na(as.numeric(sam_gene_gt$SAMPLE)) & 
                                     !grepl("^CR|^LK", as.character(sam_gene_gt$SAMPLE)), 1, 0))/(2*length(p_vec))
        maf_vec_case <- sum(ifelse(!is.na(as.numeric(sam_gene_gt$SAMPLE)) | 
                                     grepl("^CR|^LK", as.character(sam_gene_gt$SAMPLE)), 1, 0))/(2*length(p_vec))
        ##genotype matrix  
        
        sam_gene_gt$add_mod <- as.numeric(sam_gene_gt$GT)
        sam10 <- ifelse(Ex_samp_id %in% sam_gene_gt$SAMPLE, 1, 0)
        sam10[which(sam10 != 0)] <- sam_gene_gt$add_mod ##additive model
        ##account for discrepancy in vep_consequence
        #vep_con <- names(sort(table(unlist(strsplit(sam_gene_gt$vep_consequence, split = "&"))),decreasing=TRUE)[1])
        # sam_gene_gt$vep_consequence <- vep_con
        
        
        samp_vec[[m]] <- c(ftemp_tab_var_id[m],
                           unique(sam_gene_gt$gene_symbol), 
                           unique(sam_gene_gt$vep_consequence), 
                           unique(as.character(sam_gene_gt$auto_call)),
                           as.numeric(unique(sam_gene_gt$comb_score)), as.numeric(maf_vec_case),
                           as.numeric(maf_vec_cont), sam10)
      }
    }
    samp_vec_mat_uni <- do.call("rbind.data.frame", samp_vec)
    colnames(samp_vec_mat_uni) <- c("VARIANT", "gene_symbol", "vep_consequence", "auto_call", "comb_score", "coh_MAF_case",
                                    "coh_MAF_cont", Ex_samp_id)
    if(exists("samp_vec_con")){
      samp_vec_mat_con <- do.call("rbind.data.frame", samp_vec_con)
      colnames(samp_vec_mat_con) <- c("VARIANT", "gene_symbol", "vep_consequence", "auto_call", "comb_score", "coh_MAF_case",
                                      "coh_MAF_cont", Ex_samp_id)
      samp_vec_mat <- rbind.data.frame(samp_vec_mat_uni, samp_vec_mat_con)
    }else{
      samp_vec_mat <- samp_vec_mat_uni
    }
    print(dim(samp_vec_mat))
    ##Intracohort filter  : equivalent to ~ 5/1661*2 for case ; ~ 5/3209*2 for control
    ##changed to uniform intra_cohort_MAF filter
    ## change to whole cohort MAF filter to ~ 3/(5790*2)
    ##~ 3/(5771*2) after removal of duplicate samples (15 case and 4 controls)
    #     samp_vec_mat <- samp_vec_mat[as.numeric(as.character(samp_vec_mat$coh_MAF_case)) <= 0.0015 | 
    #                                    as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) <= 0.001,]
    samp_vec_mat <- samp_vec_mat[!(as.numeric(as.character(samp_vec_mat$coh_MAF_case)) >= 0.00026 | 
                                     as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) >= 0.00026),]
    #  samp_vec_mat <- samp_vec_mat[as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) <= 0.001,]
    
    if(is.null(dim(samp_vec_mat)) | dim(samp_vec_mat)[1] == 0) {
      skat_pvals <- NULL
      skat_pvals_pc12 <- NULL
      skat_pvals_sex <- NULL
      skat_pvals_sex_pc12 <- NULL
      # skat_pvals_age[[i]] <- NULL
      # skat_pvals_age_pc12[[i]] <- NULL
      # skat_pvals_age_sex[[i]] <- NULL
      # skat_pvals_age_sex_pc12[[i]] <- NULL
      #   break
      next
    }
    ##compute maf SKAT takes care of the AFs internally by assigning higher weights to rare variants 
    
    else {
      gene_mat_comb <- as.numeric(as.character(samp_vec_mat[,5]))
      gene_mat <- as.matrix(samp_vec_mat[,-c(1:7)])
      class(gene_mat) <- "numeric"
      ##for strictly burden test set rho = 1, else set rho = 0
      ##no covariate
      skat_pvals <- SKAT_run(geno_mat = gene_mat, gene = genes, x=NULL, p_vec, cust_weight = gene_mat_comb, rho = 1)
      ##only pc1 + pc2
      skat_pvals_pc12 <- SKAT_run(geno_mat = gene_mat, gene = genes, x=p_Data_noCH[,c(19:22)], p_vec, cust_weight = gene_mat_comb, rho = 1)
      ##gender
      skat_pvals_sex <- SKAT_run(geno_mat = gene_mat, gene = genes, x=p_Data_noCH$gender, p_vec, cust_weight = gene_mat_comb, rho = 1)
      ##gender + pc1 + pc2
      skat_pvals_sex_pc12 <- SKAT_run(geno_mat = gene_mat, gene = genes, x=p_Data_noCH[,c(53,19:22)], p_vec, cust_weight = gene_mat_comb, rho = 1)
      
    }
    
  } ##end of nested else loop
  DT_skat_snv_str_pc123 <- list(skat_pvals, skat_pvals_pc12, skat_pvals_sex, skat_pvals_sex_pc12)
  return(DT_skat_snv_str_pc123)
}##end of gene loop ; i loop   

##Parallelised SKAT  
library(doParallel)
library(doMC)
registerDoMC(30)
res20 <- list()

genes <- unique(fil_tab$gene_symbol)
system.time(res20 <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
{para_SKAT(genes[i])})
# res20_skat <- lapply(res20, function(x)do.call("rbind.data.frame", x[4]))
#t2 <- do.call("rbind.data.frame", res20_skat)
res_skat_para <- list()
for(i in 1:4){
  res_skat_para[[i]] <- lapply(res20, function(x)do.call("rbind.data.frame", x[i]))
}

##SKAT ERA added along with changes to SKATO; now SKATBinary is used for SKATO and SKAT_ERA
#ERA -> Adaptive resampling for conservative p-values
saveRDS(res_skat_para, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_EPIT/round2/Exome_skat_para_result_isks_epit_mgrb_combset2020_PC1234_ver4_rect_test.rds", compress = T)
