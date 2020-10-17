##cliques
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(knitr)
library(igraph)
library(ggnet)
library(intergraph)
library(network)
library(org.Hs.eg.db)
library(topGO)

`%nin%` = Negate(`%in%`)

##Cliques
##Use ASRB subtracted gene list

#ppi_res_fil_final_comb <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/ppi_res_fil_final_SKATbin_comb_str_biog_Aug31_cyto.tsv", sep = "\t", header = T,
#                                     stringsAsFactors = F)
ppi_res_fil_final_comb <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/isks_cyto_ppifisher_ASRB_comb_genes_top_new_rand_comp_score_Aug31.tsv", sep = "\t", header = T,
                                     stringsAsFactors = F)

##filter1
filt1_df <- ppi_res_fil_final_comb[ppi_res_fil_final_comb$ISKS > 1 & ppi_res_fil_final_comb$MGRB <= 14,]

##filter2: based on pvalue of ASRB genes
filt2_df <- filt1_df[filt1_df$pval_SKATbin.1 > 0.1 | is.na(filt1_df$pval_SKATbin.1),]
top_genes <- as.character(filt2_df$gene)
##filter3 : based on hypergeometric score of PPI (can be used for further fine-tuning and ranking)
##not needed for gene selection


strindb_biog_graph1 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/strindb_biog_graph_cyto.rds")

can_net1 <- igraph::as_data_frame(strindb_biog_graph1, what = "edges")
prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% top_genes 
                               & as.character(can_net1$to) %in% top_genes, ])
# prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% cep_genes 
#                                & as.character(can_net1$to) %in% cep_genes, ])

uniongraph <- igraph::graph.data.frame(prot_np_all, directed = F)

te3 <- max_cliques(uniongraph, min=4)
names(te3) <- paste("c", 1:length(te3), sep = "_")

#saveRDS(te3, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Cliques/Cliques_max101.rds", compress = T)
saveRDS(te3, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Cliques/Cliques_max64.rds", compress = T)

#####################################
##Burden test
cpx_OR_fisher <- function(ppi_res,cliq_list, case_coh_size, cont_coh_size, coh){
  ft_df <- list()
  for(i in 1:length(cliq_list)){
    ppi_res_tab <- ppi_res[ppi_res$gene %in% names(cliq_list[[i]]),]
    # ppi_res_tab[,2] <- ifelse(ppi_res_tab[,2] == 0, 1, ppi_res_tab[,2])
    inp <- c(sum(ppi_res_tab[,1]), case_coh_size - sum(ppi_res_tab[,1]) , 
             sum(ppi_res_tab[,2]), cont_coh_size - sum(ppi_res_tab[,2]))
    sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
    colnames(sim_mat) <- c("case", "cont")
    rownames(sim_mat) <- c("hits", "no_hits")
    #ft <- fisher.test(sim_mat, alternative = "greater")
    #ft <- fisher.test(sim_mat)
    ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.99)
    ft_df[[i]] <- cbind.data.frame("clique" = names(cliq_list)[i] ,"Cases" = sum(ppi_res_tab[,1]),
                                   "Controls" = sum(ppi_res_tab[,2]),
                                   "Fish_pval" = ft$p.value,"CI_lower" = ft$conf.int[1],
                                   "CI_upper" = ft$conf.int[2],
                                   "OR_Fish" = ft$estimate, "case_coh_size" = case_coh_size,
                                   "Coh" = coh, "genes" = paste(names(cliq_list[[i]]), collapse = ","))
  }
  return(ft_df)
}

##ISKS vs MGRB enrichment
isks_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_para_tab_ISKS_MGRB_C345_Aug31.rds")
isks_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Exome_para_tab_ISKS_MGRB_minusC3_Aug31.rds")


cpx_ISKS_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_mgrb_genes, te3, 1644, 3205, "ISKSvsMGRB"))
cpx_ISKS_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_mgrb_genes_noC3, te3, 1644, 3205, "ISKSvsMGRB"))
cpx_ISKS_OR_df_noC3$adj_pval <- p.adjust(cpx_ISKS_OR_df_noC3$Fish_pval, method = "bonferroni", n = length(cpx_ISKS_OR_df_noC3$Fish_pval))

############################
##Weighted Burden test(SKAT)
library(data.table)
library(VariantAnnotation)
library(stringr)
library(SKAT)

###get input files and process
##Experiments

fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_all_fields_rnd3_freeze.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)

##Additive model
Ex_samp_id <- unique(fil_tab$SAMPLE)

######get phenotype data to control for age and sex and PC's; 

QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv", 
                      header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %in% "TRUE",]
QC2_dat_pass$isFemale <- ifelse(QC2_dat_pass$f_stat < 0.2, 1, 
                                ifelse(QC2_dat_pass$f_stat > 0.8, 0, 2))
##note p_Data_PC_comb does not change
#p_Data <- read.table("~/RVAS/ISKS_MGRB_gender_age_3665.tsv", header = T, sep = "\t")
#p_Data <- read.table("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_pca_combset.scores.tsv", header = T, sep = "\t", stringsAsFactors = F)
##for future QC with prob.NFE correction use the following file and adjust the column for gender from 39 to 53
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
fil_tab <- fil_tab[fil_tab$SAMPLE %in% Ex_samp_id,]
##save table
#write.table(fil_tab, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv", 
#            sep = "\t", row.names = F, quote = F)

##Add gender information
p_Data_noCH$gender <- QC2_dat_pass[match(p_Data_noCH$rect_sam, QC2_dat_pass$new_sampleid), 21]

#binary phenotype vector 
p_vec <- ifelse(!is.na(as.numeric(as.character(p_Data_noCH$rect_sam))) | grepl("^CR|^LK",as.character(p_Data_noCH$rect_sam)), 1, 0)

##SKAT null function with customised covariate
SKAT_fun_null <- function(x=NULL, p_vec){
  if(is.null(x)){
    obj_N <- SKAT::SKAT_Null_Model(p_vec ~ 1,out_type="D")
  }
  # else if(is.integer(x) || is.numeric()){
  else if(is.null(dim(x))){
    obj_N <- SKAT::SKAT_Null_Model(p_vec ~ x,out_type="D")
  }
  else if(dim(x)[2] > 1){
    nul_for <- as.formula(paste("p_vec", paste(colnames(x), collapse = " + "), sep = " ~ "))
    obj_N <- SKAT::SKAT_Null_Model(nul_for, data = p_Data_noCH, out_type="D")
    #  obj_N <- SKAT_Null_Model(nul_for, data = x, out_type="D")
  }
  return(obj_N)
}

##SKAT function for SKAT, SKATBinary, SKATO, SKAT_ERA
SKAT_run <- function(geno_mat, gene, x=NULL, p_vec, cust_weight = NULL, rho = NULL){
  
  if(dim(geno_mat)[2] > 1 ){
    null_d <- SKAT_fun_null(x, p_vec)
    pval_SKAT <- SKAT::SKAT(t(geno_mat), null_d, weights = cust_weight, r.corr = rho)$p.value
    pval_SKATbin <- SKAT::SKATBinary(t(geno_mat), null_d, method = "Burden", weights = cust_weight)$p.value
    ##ER: Effective resampling for best p-value calculation
    pval_SKATO <- SKAT::SKATBinary(t(geno_mat),null_d, method='SKATO', weights = cust_weight)$p.value
    pval_SKATera <- SKAT::SKATBinary(t(geno_mat),null_d, method.bin="ER.A", weights = cust_weight)$p.value
    skat_pv <- cbind.data.frame("eg_ID" = gene, pval_SKAT, pval_SKATbin, pval_SKATO, pval_SKATera)
    return(skat_pv)
  }
  else if(dim(geno_mat)[2] == 1){
    null_d <- SKAT::SKAT_fun_null(x, p_vec)
    pval_SKAT <- SKAT::SKAT(geno_mat, null_d, weights = cust_weight, r.corr = rho)$p.value
    pval_SKATbin <- SKAT::SKATBinary(t(geno_mat), null_d, method = "Burden", weights = cust_weight)$p.value
    ##ER: Effective resampling for best p-value calculation
    pval_SKATO <- SKAT::SKATBinary(geno_mat,null_d, method='SKATO', weights = cust_weight)$p.value
    pval_SKATera <- SKAT::SKATBinary(geno_mat,null_d, method.bin="ER.A", weights = cust_weight)$p.value
    skat_pv <- cbind.data.frame("eg_ID" = gene, pval_SKAT, pval_SKATbin, pval_SKATO, pval_SKATera)
    
    return(skat_pv)
  }
}

DT_skat_snv_str_pc123 <- list()
skat_pvals <- list()
skat_pvals_pc12 <- list()
skat_pvals_sex <- list()
skat_pvals_sex_pc12 <- list()

##parallel process genes (remove C3)
para_SKAT <- function(genes, cpx_name){
  ftemp_tab <- fil_tab[fil_tab$gene_symbol %in% names(genes),]
  ftemp_tab <- ftemp_tab[ftemp_tab$comb_score >= 5.6 & as.numeric(ftemp_tab$VAF) >= 0.35 &
                           ftemp_tab$auto_call %nin% "C3", ]
  
  print(max(ftemp_tab$comb_score))
  if(dim(ftemp_tab) == 0 ){
    next
  }
  else{
    ftemp_tab_var_id <- unique(ftemp_tab$VARIANT)
    samp_vec <- list()
    for(m in 1:length(ftemp_tab_var_id)){
      sam_gene_gt <- ftemp_tab[ftemp_tab$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,127:128)]
      sam_gene_gt <- unique(sam_gene_gt)
      if(dim(sam_gene_gt)[1] > 1 & length(unique(sam_gene_gt$vep_consequence)) > 1){
        # if(dim(sam_gene_gt)[1] > 1 & length(unique(sam_gene_gt$vep_consequence)) > 1 & length(unique(sam_gene_gt$comb_score)) > 1){
        vep_con <- unique(sam_gene_gt$vep_consequence)
        samp_vec_con <- list()
        for(k in 1:length(vep_con)){
          sam_gene_gt_con <- sam_gene_gt[sam_gene_gt$vep_consequence %in% vep_con[k],]
          sam_gene_gt_con <- unique(sam_gene_gt_con)
          maf_vec_cont <- sum(dim(sam_gene_gt_con[grepl("^[ABZ]",sam_gene_gt_con$SAMPLE),])[1])/(2*length(p_vec))
          maf_vec_case <- sum(dim(sam_gene_gt_con[!grepl("^[ABZ]",sam_gene_gt_con$SAMPLE),])[1])/(2*length(p_vec))
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
        maf_vec_cont <- sum(ifelse(is.na(as.numeric(sam_gene_gt$SAMPLE)), 1, 0))/(2*length(p_vec))
        maf_vec_case <- sum(ifelse(!is.na(as.numeric(sam_gene_gt$SAMPLE)) | 
                                     grepl("^CR|^LK", as.character(sam_gene_gt$SAMPLE)), 1, 0))/(2*length(p_vec))
        #maf_vec <- (maf_vec_cont + maf_vec_case)/(2*(1572 + 1110))
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
    colnames(samp_vec_mat_uni) <- c("VARIANT", "cpx_name", "vep_consequence", "auto_call", "comb_score", "coh_MAF_case",
                                    "coh_MAF_cont", Ex_samp_id)
    if(exists("samp_vec_con")){
      samp_vec_mat_con <- do.call("rbind.data.frame", samp_vec_con)
      colnames(samp_vec_mat_con) <- c("VARIANT", "cpx_name", "vep_consequence", "auto_call", "comb_score", "coh_MAF_case",
                                      "coh_MAF_cont", Ex_samp_id)
      samp_vec_mat <- rbind.data.frame(samp_vec_mat_uni, samp_vec_mat_con)
    }else{
      samp_vec_mat <- samp_vec_mat_uni
    }
    print(dim(samp_vec_mat))
    ##Intracohort filter  : equivalent to ~ 5/1661*2 for case ; ~ 5/3209*2 for control
    ##changed to uniform intra_cohort_MAF filter
    ## change to entire cohort MAF filter to 3/(4849*2)##latest Aug.31
    #     samp_vec_mat <- samp_vec_mat[as.numeric(as.character(samp_vec_mat$coh_MAF_case)) <= 0.0015 | 
    #                                    as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) <= 0.001,]
    samp_vec_mat <- samp_vec_mat[!(as.numeric(as.character(samp_vec_mat$coh_MAF_case)) >= 0.00035 | 
                                     as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) >= 0.00035),]
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
      skat_pvals <- SKAT_run(geno_mat = gene_mat, gene = cpx_name, x=NULL, p_vec, cust_weight = gene_mat_comb, rho = 1)
      ##only pc1 + pc2
      skat_pvals_pc12 <- SKAT_run(geno_mat = gene_mat, gene = cpx_name, x=p_Data_noCH[,c(19:22)], p_vec, cust_weight = gene_mat_comb, rho = 1)
      ##gender
      skat_pvals_sex <- SKAT_run(geno_mat = gene_mat, gene = cpx_name, x=p_Data_noCH$gender, p_vec, cust_weight = gene_mat_comb, rho = 1)
      ##gender + pc1 + pc2
      skat_pvals_sex_pc12 <- SKAT_run(geno_mat = gene_mat, gene = cpx_name, x=p_Data_noCH[,c(54,19:22)], p_vec, cust_weight = gene_mat_comb, rho = 1)
      
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
system.time(res20 <- foreach(i=1:length(te3), .errorhandling = 'remove') %dopar% 
{para_SKAT(te3[[i]], names(te3)[i])})
res_skat_para <- list()
for(i in 1:4){
  res_skat_para[[i]] <- lapply(res20, function(x)do.call("rbind.data.frame", x[i]))
}


SKAT_res_df <- do.call("rbind.data.frame", res_skat_para[[4]])
SKAT_res_df$SKATbin_p_adj <- p.adjust(SKAT_res_df$pval_SKATbin, method = "bonferroni", n = length(SKAT_res_df$pval_SKATbin))

##combine burden(noC3) and SKAT test(C3-5)
comb_burd_SKAT <- cpx_ISKS_OR_df_noC3

comb_burd_SKAT[,12:14] <- SKAT_res_df[match(comb_burd_SKAT$clique, SKAT_res_df$eg_ID),c(1,3,6)]
comb_burd_SKAT <- comb_burd_SKAT[complete.cases(comb_burd_SKAT),]
rownames(comb_burd_SKAT) <- comb_burd_SKAT$clique
#comb_burd_SKAT[comb_burd_SKAT$SKATbin_p_adj < 0.01,-c(5:6,8:9,12)]
#write.table(comb_burd_SKAT, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Cliques/Clique_OR_SKAT.tsv",
#            sep = "\t", row.names = F, quote = F)
write.table(comb_burd_SKAT, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Cliques/Clique_OR_SKAT_968.tsv",
            sep = "\t", row.names = F, quote = F)
###Add other attributes to comb_burd_SKAT (clique similarity)(see cliques_for_paper.Rmd)

