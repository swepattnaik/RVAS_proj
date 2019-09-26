##Note that this analysis include singletons
##changes included on Aug-06-2019:
#added cohort_MAF computation feature; 
#the rare variants are filtered using cohort_MAF filter of lte 0.005
#input is better filtered compared to Aug4 input
##Aug12: SKAT made feasible for combined ISKS and RISC samples

.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )


##Using D data structure

library(data.table)
library(VariantAnnotation)
library(stringr)
library(SKAT)

###get input files and process


##Experiments

#D_tab <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/SKAT_inp_filt_var_noCH_C5eqC4_noNA_Exome_gt.rds")
#fil_tab <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_noNAs_auto_4Aug.rds")
#fil_tab <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_noNAs_auto_6Aug.rds")
#fil_tab <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_auto_7Aug.rds")
#fil_tab <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_12Aug.rds")
##DT_exome_clean3.R was used to generate 05Sept_rect file
#fil_tab <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_05Sept_rect.rds")
fil_tab <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_05Sept_rect_ASP.rds")
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)
#samp_ID <- colnames(D_tab[,-c(1:5)])
#samp_ID_match <- samp_ID
#samp_ID_match <- gsub("isks", "", samp_ID) 
##Additive model
Ex_samp_id <- unique(fil_tab$SAMPLE)
Ex_var_id <- unique(fil_tab$VARIANT)

######get phenotype data to control for age and sex and PC's; 
##note p_Data_PC_comb does not change
#p_Data <- read.table("~/RVAS/ISKS_MGRB_gender_age_3665.tsv", header = T, sep = "\t")
p_Data <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/02_ISKS_MGRB_1000G_pca.scores.tsv", header = T, sep = "\t")
p_Data$sample <- gsub("isks", "", p_Data$sample)
p_Data$sample <- gsub("risc", "", p_Data$sample)
p_Data_noCH <- p_Data[as.character(p_Data$sample) %in% Ex_samp_id,]
##drop MGRB sample BAAUD; not in latest call
#samp_ID_match <- samp_ID_match[grep("BAAUD", samp_ID_match, invert = T)]
#p_Data_noCH <- p_Data_noCH[match(samp_ID_match, p_Data_noCH$sample),]
Ex_samp_id <- Ex_samp_id[match(p_Data_noCH$sample, Ex_samp_id)]
##remove from gene matrix
#col_rm <- which(colnames(D_tab) ==  "BAAUD")
#D_tab <- D_tab[,-c(col_rm)]
##Add gender information
p_Data2 <- read.table("~/RVAS/ISKS_MGRB_gender_age_3665.tsv", header = T, sep = "\t", stringsAsFactors = F)
p_Data3 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/RISC_pheno.tsv", sep = "\t", header = T,
                      stringsAsFactors = F)
p_Data3$cohort <- c("RISC")
#p_Data3$cohort <- as.factor(p_Data3$cohort)
p_Data3 <- p_Data3[,-3]
p_Data3 <- p_Data3[,c(1,4,2,3)]
colnames(p_Data3) <- colnames(p_Data2)
p_Data3$sampleID <- gsub("risc", "", p_Data3$sampleID)
p_Data3$isFemale <- ifelse(p_Data3$isFemale %in% "Female", 1, 0)

p_Data4 <- rbind.data.frame(p_Data2, p_Data3, stringsAsFactors = F)
  
p_Data_noCH$gender <- p_Data4[match(p_Data_noCH$sample, p_Data4$sampleID), 4]

#binary phenotype vector 
p_vec <- ifelse(!is.na(as.numeric(as.character(p_Data_noCH$sample))) | grepl("^CR",as.character(p_Data_noCH$sample)), 1, 0)


#pca_pop <- read.delim("~/RVAS/pop_PCA/mgrb_isks_1000g_pca.eigenvec", sep = "\t", header = T) ##unit normal eigen vectors

##use scale eigenvectors as the enriched genes seem to be olfactory receptors; but p_Data already has pca scores which
##has the appropriate variance explained.
# pca_pop <- read.delim("~/RVAS/pop_PCA/pc_20_loading.tsv", sep = "\t", header = T)
# pca_pop_is_mg <- pca_pop[pca_pop$X %in% p_Data_noCH$sample, ]
# pca_pop_is_mg <- pca_pop_is_mg[match(p_Data_noCH$sample, pca_pop_is_mg$X),]
# p_Data_noCH <- cbind.data.frame(p_Data_noCH, pca_pop_is_mg)


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

#D_tab <- fread("~/RVAS/isks.variants.tsv", sep = "\t", header = T)
#dim(D_tab)

##list to store output
#shard_skat_pphen_pc12 <- list()
#shard_skat_snv_pc123 <- list()
DT_skat_snv_str_pc123 <- list()

##HAIL output
#fil_tab <- read.table("~/RVAS/shard_sub_tier3/CGC_genes/cgc_hail_ann_sub_pphen05_NFE_ver2.tsv", sep = "\t", header = T)
#genes <- levels(fil_tab$ENSG)
#fil_tab <- D_tab
genes <- unique(fil_tab$gene_symbol)

    skat_pvals <- list()
    skat_pvals_pc12 <- list()
    skat_pvals_sex <- list()
    skat_pvals_sex_pc12 <- list()
#    skat_pvals_age <- list()
 #   skat_pvals_age_pc12 <- list()
 #   skat_pvals_age_sex <- list()
#    skat_pvals_age_sex_pc12 <- list()
    
    for(i in 1:length(genes)){
      print(i)
      print(genes[i])
      ##process genes in a loop
      ftemp_tab <- fil_tab[fil_tab$gene_symbol %in% genes[i],]
      ftemp_tab <- ftemp_tab[ftemp_tab$comb_score >= 5 & as.numeric(ftemp_tab$VAF) >= 0.25, ]
      
      
      
      # ##filter by AF < 0.01
      # gene_mat <- gene_mat[gene_mat$af < 0.01,]
      # 
      # ##filter only SNVs
      # #gene_mat <- gene_mat[gene_mat$type == 83,] #comment to include indels
      ##remove all missense variants with score less than 5
      ##caution this may have removed indels with eigenphred scores below 5
      ##use this filter only for mis-sense mutations (this is included in the filter_HAIL script now)
      #ftemp_tab <- ftemp_tab[ftemp_tab$comb_score >= 5, ] ##comment for relaxed run
     # ftemp_tab <- ftemp_tab[ftemp_tab$comb_score > 5, ] ##for latest run; filtering out inframe indels
      
      print(max(ftemp_tab$comb_score))
      if(dim(ftemp_tab) == 0 ){
        next
      }
      else{
        ftemp_tab_var_id <- unique(ftemp_tab$VARIANT)
        samp_vec <- list()
        for(m in 1:length(ftemp_tab_var_id)){
         # sam_gene_gt <- ftemp_tab[grepl(ftemp_tab_var_id[m],
          #                               ftemp_tab$VARIANT),][,c(1:3,9,11,120,121)]
          #sam_gene_gt <- ftemp_tab[ftemp_tab$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,129:130)]
          sam_gene_gt <- ftemp_tab[ftemp_tab$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,130:131)]
          #sam_gene_gt$add_mod <- ifelse(sam_gene_gt$GT %in% "0/1", 1, ifelse(sam_gene_gt$GT %in% "0/0", 0, 2))
          sam_gene_gt <- unique(sam_gene_gt)
          ##compute cohort specific MAF
      #    maf_vec_cont <- ifelse(is.na(as.numeric(sam_gene_gt$SAMPLE)), 1, 0)
       #   maf_vec_case <- ifelse(!is.na(as.numeric(sam_gene_gt$SAMPLE)) | 
       #                            grepl("^CR", as.character(sam_gene_gt$SAMPLE)), 1, 0)
       #   maf_vec_cont <- round(sum(as.numeric(sam_gene_gt$GT)*(maf_vec_cont))/(2*1714), digits = 5)
       #   maf_vec_case <- round(sum(as.numeric(sam_gene_gt$GT)*(maf_vec_case))/(2*1095), digits = 5)
         
     #      maf_vec_cont <- round(sum(as.numeric(sam_gene_gt$GT)*(maf_vec_cont))/(2*1572), digits = 5) ##After removal of ASPREE cancer
      ##changed 1095 to 1110 on Aug12
         # maf_vec_case <- round(sum(as.numeric(sam_gene_gt$GT)*(maf_vec_case))/(2*1110), digits = 5)
         # maf_vec <- maf_vec_cont + maf_vec_case
          
           maf_vec_cont <- sum(ifelse(is.na(as.numeric(sam_gene_gt$SAMPLE)), 1, 0))/(2*1572)
           maf_vec_case <- sum(ifelse(!is.na(as.numeric(sam_gene_gt$SAMPLE)) | 
                                    grepl("^CR", as.character(sam_gene_gt$SAMPLE)), 1, 0))/(2*1110)
           #maf_vec <- (maf_vec_cont + maf_vec_case)/(2*(1572 + 1110))
          ##genotype matrix  
          sam_gene_gt$add_mod <- as.numeric(sam_gene_gt$GT)
          sam10 <- ifelse(Ex_samp_id %in% sam_gene_gt$SAMPLE, 1, 0)
          sam10[which(sam10 != 0)] <- sam_gene_gt$add_mod ##additive model
          
          # names(sam10) <- Ex_samp_id
          samp_vec[[m]] <- c(ftemp_tab_var_id[m],
                             unique(sam_gene_gt$gene_symbol), 
                             unique(sam_gene_gt$vep_consequence), 
                             unique(as.character(sam_gene_gt$auto_call)),
                             as.numeric(unique(sam_gene_gt$comb_score)), as.numeric(maf_vec_case),
                             as.numeric(maf_vec_cont), sam10)
        }
       
        samp_vec_mat <- do.call("rbind.data.frame", samp_vec)
       # colnames(samp_vec_mat) <- c("VARIANT", "gene_symbol", 
        #                            "vep_consequence", "auto_call", "comb_score", Ex_samp_id)
        
      
      ##remove singletons (DT wants to retain singletons for clinical benefit)
     # sing_row <- apply(gene_mat, 1, sing_rem_index)
     # rname <- rownames(gene_mat[sing_row == 0,])
     # gene_mat <- gene_mat[sing_row == 0,]
      # print(table(gene_mat))
     # gene_mat <- as.matrix(gene_mat)
      #print(i)
      print(dim(samp_vec_mat))
      colnames(samp_vec_mat) <- c("VARIANT", "gene_symbol", 
                                  "vep_consequence", "auto_call", "comb_score", "coh_MAF_case",
                                  "coh_MAF_cont", Ex_samp_id)
     # samp_vec_mat <- samp_vec_mat[as.numeric(as.character(samp_vec_mat$coh_MAF)) <= 0.005,]
     ##Intracohort filter  : equivalent to ~ 3/1110*2 for case ; ~ 3/1562*2
       samp_vec_mat <- samp_vec_mat[as.numeric(as.character(samp_vec_mat$coh_MAF_case)) <= 0.0015 | 
                                      as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) <= 0.001,]
    #  samp_vec_mat <- samp_vec_mat[as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) <= 0.001,]
      
      if(is.null(dim(samp_vec_mat)) | dim(samp_vec_mat)[1] == 0) {
        skat_pvals[[i]] <- NULL
        skat_pvals_pc12[[i]] <- NULL
        skat_pvals_sex[[i]] <- NULL
        skat_pvals_sex_pc12[[i]] <- NULL
       # skat_pvals_age[[i]] <- NULL
       # skat_pvals_age_pc12[[i]] <- NULL
       # skat_pvals_age_sex[[i]] <- NULL
       # skat_pvals_age_sex_pc12[[i]] <- NULL
        #   break
        next
      }
      ##compute maf SKAT takes care of the AFs internally by assigning higher weights to rare variants 
      # maf <- rowMeans(gene_genomat, na.rm = T)/2
      ##rare variant cutoff set to MAF < 0.05
      # gene_genomat <- gene_genomat[which(maf <= 0.05),]
      
      else {
        gene_mat_comb <- as.numeric(as.character(samp_vec_mat[,5]))
        gene_mat <- as.matrix(samp_vec_mat[,-c(1:7)])
        class(gene_mat) <- "numeric"
        ##for strictly burden test set rho = 1, else set rho = 0
        ##no covariate
        skat_pvals[[i]] <- SKAT_run(geno_mat = gene_mat, gene = genes[i], x=NULL, p_vec, cust_weight = gene_mat_comb, rho = 1)
        ##only pc1 + pc2
        skat_pvals_pc12[[i]] <- SKAT_run(geno_mat = gene_mat, gene = genes[i], x=p_Data_noCH[,c(19:21)], p_vec, cust_weight = gene_mat_comb, rho = 1)
        ##gender
        skat_pvals_sex[[i]] <- SKAT_run(geno_mat = gene_mat, gene = genes[i], x=p_Data_noCH$gender, p_vec, cust_weight = gene_mat_comb, rho = 1)
        ##gender + pc1 + pc2
        skat_pvals_sex_pc12[[i]] <- SKAT_run(geno_mat = gene_mat, gene = genes[i], x=p_Data_noCH[,c(39,19:21)], p_vec, cust_weight = gene_mat_comb, rho = 1)
        
      }
  
    } ##end of nested else loop
    
      }##end of gene loop ; i loop   
  
    
    DT_skat_snv_str_pc123 <- list(skat_pvals, skat_pvals_pc12, skat_pvals_sex, skat_pvals_sex_pc12)
  
    #saveRDS(DT_skat_snv_str_pc123, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_skat_wsing_load123_noCH_C5eqC4_noNA_gt.rds", compress = T)   
   # saveRDS(DT_skat_snv_str_pc123, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_skat_wsing_load123_noCH_C5eqC4_noNA_gt_4Aug.rds", compress = T)    
   # saveRDS(DT_skat_snv_str_pc123, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_skat_wsing_load123_noCH_C5eqC4_noNA_gt_6Aug.rds", compress = T)
   # saveRDS(DT_skat_snv_str_pc123, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_skat_wsing_load123_noCH_C5eqC4_nonmds_gt_7Aug.rds", compress = T)   
   # saveRDS(DT_skat_snv_str_pc123, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_skat_wsing_load123_noCH_C5eqC4_nonmds_gt_isksrisc_12Aug_rerun.rds", compress = T)
   ##SKAT ERA added along with changes to SKATO; now SKATBinary is used for SKATO and SKAT_ERA
    #ERA -> Adaptive resampling for conservative p-values
   # saveRDS(DT_skat_snv_str_pc123, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_skat_wsing_load123_noCH_C5eqC4_nonmds_gt_isksrisc_12Aug_rerun.rds", compress = T)
   # saveRDS(DT_skat_snv_str_pc123, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_skat_wsing_load123_noCH_C5eqC4_nonmds_gt_isksrisc_sept05_rect.rds", compress = T) 
    saveRDS(DT_skat_snv_str_pc123, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Exome_skat_wsing_load123_noCH_C5eqC4_nonmds_gt_isksrisc_sept05_rect_ASP_cmaf.rds", compress = T)
    