##Make PID file
##with DT criteria: VAF >= 0.35 & Eigenphred >= 5.6
##PID file: Apr29-2020
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

summarise_var_per_samp <- function(variant_set){
  samvar_list <- list()
  sam_id <- unique(variant_set$SAMPLE)
  for(i in 1:length(sam_id)){
    v1 <- variant_set[variant_set$SAMPLE %in% sam_id[i],]
    # v1 <- v1[v1$VAF >= 0.35 & v1$comb_score >= 3,] ##comment for unfiltered variant set
    v1 <- v1[order(v1$comb_score, decreasing = T),]
    if(!is.null(v1)){
      ntot <- dim(v1)[1]
      nC5 <- sum(ifelse(v1$auto_call %in% "C5", 1, 0))
      geneC5 <- paste(v1[v1$auto_call %in% "C5",]$gene_symbol, collapse = ";")
      nC4 <- sum(ifelse(v1$auto_call %in% "C4", 1, 0))
      geneC4 <- paste(v1[v1$auto_call %in% "C4",]$gene_symbol, collapse = ";")
      nC3 <- sum(ifelse(v1$auto_call %in% "C3", 1, 0))
      geneC3 <- paste(v1[v1$auto_call %in% "C3",]$gene_symbol, collapse = ";")
      
      
    }else{    ntot = 0
    nC5 <- 0
    geneC5 <- NULL
    nC4 <- 0
    geneC4 <- NULL
    nC3 <- 0
    geneC4 <- NULL
    }
    gene_list <- list("EXT1_2" = c("EXT1", "EXT2"), 
                      "IDH1_2" = c("IDH1", "IDH2"),
                      "TP53" = c("TP53"),
                      "BRCA1_2_PALB2" = c("BRCA1", "BRCA2", "PALB2"),
                      "MMR_APC_MUTYH" = c("MMR", "APC", "MUTYH"),
                      "ATM_ATR_CHEK2_ARF_TP53" = c("ATM", "ATR", "CHEK2", "ARF", "TP53"),
                      "NF1_SDH" = c("NF1", "SDH"),
                      "Sheltrin_main" = c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "ACD"),
                      "Sheltrin_extn" = c("ACD", "POT1", "TERF1", "TERF2", "TERF2IP", "TINF2", "ATM", 
                                          "BAG3", "BLM", "BRCA1", "CALD1", "CLK3", "DCLRE1B", "FANCD2", 
                                          "FBXO4", "HSPA4", "KIAA1191", "MRE11A", "NBN", "PINX1", "PRKDC", 
                                          "RAD50", "SLX4", "STUB1", "TNKS", "TNKS2", "U2AF2", "UCHL1", 
                                          "WRN", "XRCC5", "XRCC6"))
    gene_list_match <- lapply(gene_list, function(x) ifelse(x %in% v1$gene_symbol, 1, 0))
    sum_list <- unlist(lapply(gene_list_match, function(x)sum(x)))
    sum_list <- as.data.frame(t(sum_list))
    samvar_list[[i]] <- cbind.data.frame("Total.burden" = ntot, "nC5" = nC5, "C5_genes" = geneC5,
                                         "nC4" = nC4, "C4_genes" = geneC4,
                                         "nC3" = nC3, "C3_genes" = geneC3, sum_list, "SAMPLE" = sam_id[i])
  }
  return(samvar_list)
}

##filtered set
fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/all_epit_mgrb_combset2020_variants_filt_all_fields_risk_factor.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
Ex_samp_id <- unique(fil_tab$SAMPLE)
#remove duplicates
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

p_Data <- read.delim("~/RVAS/Epi_set_2020/pop_PCA/MGRB_EPIT_1000G_combset_pca.scores_clustered.tsv", header = T, sep = "\t", stringsAsFactors = F)
p_Data <- p_Data[p_Data$sample %in% QC2_epit_mgrb_pass$new_sampleid,]
p_Data_noCH <- p_Data[as.character(p_Data$sample) %in% Ex_samp_id,]

Ex_samp_id <- Ex_samp_id[match(p_Data_noCH$sample, Ex_samp_id)]
##filter out QC fail cases
fil_tab <- fil_tab[fil_tab$SAMPLE %in% Ex_samp_id,]

dim(fil_tab)

fil_tab_epit_set <- fil_tab[fil_tab$is_case == 1,]
##use this for QC testing not in final (Discovery set topSKAT300)
#top300_SKAT <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Sep05_rect_ASP_graph_final_PPI_comb_GO_cmaf_new_score_str.tsv",
#                          sep = "\t", header = T, stringsAsFactors = F)

##use this in final
##SKATBinary
top300_SKAT <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/round2/ppi_res_fil_final_EPIT_SKATbin_wt_fisher.tsv",
                          sep = "\t", header = T, stringsAsFactors = F)
cgc_genes <- read.delim("~/RVAS/cancer_gene_census_hg37.csv", sep = ",", header = T, stringsAsFactors = F)
#remove fusion genes from cgc list
cgc_genes$Gene.Symbol[(grep("fusion", cgc_genes$Role.in.Cancer))]
sheltrin_complex <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "ACD")
Telo_extension <- c("TIMELESS", "TIPIN", "FANCM", "BRCA1", "BLM") ##replication stress is a source of telomere recombination
Sheltrin_comp_extn = c("ACD", "POT1", "TERF1", "TERF2", "TERF2IP", "TINF2", "ATM", 
                       "BAG3", "BLM", "BRCA1", "CALD1", "CLK3", "DCLRE1B", "FANCD2", 
                       "FBXO4", "HSPA4", "KIAA1191", "MRE11A", "NBN", "PINX1", "PRKDC", 
                       "RAD50", "SLX4", "STUB1", "TNKS", "TNKS2", "U2AF2", "UCHL1", 
                       "WRN", "XRCC5", "XRCC6")
DT_genes <- read.table("~/RVAS/DT_genes_sep2018.txt", sep = "", stringsAsFactors = F) ##includes POT1 interactors
gene_sub <- unique(c(top300_SKAT$gene[1:300], cgc_genes$Gene.Symbol, sheltrin_complex, Telo_extension, DT_genes$V1, Sheltrin_comp_extn))
comb_set_filt <- fil_tab_epit_set[fil_tab_epit_set$gene_symbol %in% gene_sub,]
comb_set_filt$gnomad_AF <- gsub("\\[|\\]", "", comb_set_filt$gnomad_AF)
comb_set_filt <- comb_set_filt[comb_set_filt$VAF >= 0.35 & comb_set_filt$comb_score >= 5.6,]

##Add MAF based filter too: DT May_12_2020

p_vec <- ifelse(grepl("^[ABZ]",as.character(Ex_samp_id)), 0, 1)
#####make genotype matrix function
make_geno_mat <- function(ftemp_file){
  ftemp_tab_var_id <- unique(ftemp_file$VARIANT)
  samp_vec <- list()
  for(m in 1:length(ftemp_tab_var_id)){
    sam_gene_gt <- ftemp_file[ftemp_file$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,82,125:126)]
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
  ## change to MAF filter to 3/1661*2
  #     samp_vec_mat <- samp_vec_mat[as.numeric(as.character(samp_vec_mat$coh_MAF_case)) <= 0.0015 | 
  #                                    as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) <= 0.001,]
  samp_vec_mat <- samp_vec_mat[!(as.numeric(as.character(samp_vec_mat$coh_MAF_case)) >= 0.00037 | 
                                   as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) >= 0.00037),]
  
  return(samp_vec_mat)
  
}

gene_var_maf <- function(gene_sym, var_file){
  
  #cut_off <- seq(0,10)
  #  score_scale <- c((0.001:50)^2/100,0.5,1000)
  # p_val_cutoff <- list()
  # for(i in 1:length(cut_off))  {
  ftemp_tab <- var_file[var_file$gene_symbol %in% gene_sym,] 
  ftemp_tab <- ftemp_tab[ftemp_tab$VAF >= 0.35 & ftemp_tab$comb_score >= 5.6,]
  #  if(dim(ftemp_tab)[1] > 1){
  samp_vec_mat <- make_geno_mat(ftemp_tab)
  ftemp_tab_mod <- samp_vec_mat[,c(1:7)]
  ftemp_tab_mod$comb_score <- as.numeric(as.character(ftemp_tab_mod$comb_score))
  ftemp_tab_mod$auto_call <- as.character(ftemp_tab_mod$auto_call)
  #   }
  #   else {
  #     next
  #   }
  return(ftemp_tab_mod)
}


##Parallelised SKAT  
library(doParallel)
library(doMC)
registerDoMC(30)

res20 <- list()

genes <- unique(comb_set_filt$gene_symbol)
system.time(res20 <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
{gene_var_maf(genes[i], comb_set_filt)})
gene_maf_pass_df <- do.call("rbind.data.frame", res20)
comb_set_filt1 <- comb_set_filt[comb_set_filt$VARIANT %in% gene_maf_pass_df$VARIANT,]
write.table(comb_set_filt1, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/round2/EPIT_CGC_topSKATBin300_repstress_potint_VARIANTS_filt_combset2020_rmdup.tsv",
           sep = "\t", row.names = F, quote = F)
var_per_samp_epit <- summarise_var_per_samp(comb_set_filt1)
df_sam_filt_isks_comb <- do.call("rbind.data.frame", var_per_samp_epit)

#######
##Add phenotypes

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
#QC2_epit_pass <- QC2_epit[QC2_epit$passes_qc2 %in% "TRUE",]

comb_pheno <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set/EPIT_MGRB_2020/qc0_final_all_annot.tsv",
                         sep = "\t", header = T, stringsAsFactors = F)

comb_pheno <- unique(comb_pheno)
comb_pheno <- comb_pheno[!is.na(comb_pheno$Sample_Id),]

comb_pheno[,19:20] <- QC2_epit[match(comb_pheno$Sample_Id, QC2_epit$new_sampleid), 2:3]

comb_pheno[,21:37] <- df_sam_filt_isks_comb[match(comb_pheno$Sample_Id, df_sam_filt_isks_comb$SAMPLE),c(1:17)]

write.table(comb_pheno, 
            "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/round2/EPIT_PID_combset2020_CGC_skatBin_repstress_potint_predNFE_rmdup.tsv",
            sep = "\t", row.names = F, quote = F)
