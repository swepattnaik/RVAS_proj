##Make PID file
##with DT criteria: VAF >= 0.35 & Eigenphred >= 5.6
##PID file: Apr29-2020; amended by denia in Aug2020
##include mitotic checkpoint genes
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

`%nin%` = Negate(`%in%`)

library(dplyr)
library(stringr)
library(readxl)
mirabello_genes <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/mirabello_genes.txt",
                              sep = "", header = F, stringsAsFactors = F)
##remove trailing and leading white spaces
mira_genes <- mirabello_genes %>% 
  mutate(V1 = str_trim(mirabello_genes$V1, side = "both"))

##mitotic check point genes
mito_chk_point = read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/mitotic_checkpoint/mitocheck_point.txt",
sep = "", header = F, skip = 1, stringsAsFactors = F)

##centrosome maturation complex
centrosome_mat <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/Centrosome maturation.xlsx",
                           sheet = 1)
centrosome_mat <- as.data.frame(centrosome_mat)
centrosome_mat <- centrosome_mat[centrosome_mat$MoleculeType %nin% "Chemical Compounds",]

##cep_haus C3C4C5
cep_haus_all <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CEP_HAUS_C345.txt",
                           sep = "", header = F, stringsAsFactors = F)

summarise_var_per_samp <- function(variant_set){
  samvar_list <- list()
  sam_id <- unique(variant_set$SAMPLE)
  for(i in 1:length(sam_id)){
    v1 <- variant_set[variant_set$SAMPLE %in% sam_id[i],]
    # v1 <- v1[v1$VAF >= 0.35 & v1$comb_score >= 3,] ##comment for unfiltered variant set
    v1 <- v1[order(v1$comb_score, decreasing = T),]
    if(!is.null(v1)){
      ntot <- dim(v1)[1]
      nC5 <- sum(ifelse(grepl("C5", v1$auto_call), 1, 0))
      geneC5 <- paste(v1[grepl("^C5$", v1$auto_call),]$gene_symbol, collapse = ";")
      if(sum(ifelse(grepl("C5_ar", v1$auto_call), 1, 0)) > 0){
        geneC5_ar <- paste0(v1[grepl("C5_ar", v1$auto_call),]$gene_symbol, "_C5_ar", collapse = ";") 
        geneC5 <- paste(c(geneC5,geneC5_ar), collapse = ";")
      }
      nC4 <- sum(ifelse(grepl("C4", v1$auto_call), 1, 0))
      geneC4 <- paste(v1[grepl("^C4$", v1$auto_call, perl = T),]$gene_symbol, collapse = ";")
      if(sum(ifelse(grepl("C4_ar", v1$auto_call), 1, 0)) > 0){
        geneC4_ar <- paste0(v1[grepl("C4_ar", v1$auto_call),]$gene_symbol, "_C4_ar", collapse = ";")
        geneC4 <- paste(c(geneC4,geneC4_ar), collapse = ";")
      }
      nC3 <- sum(ifelse(grepl("C3", v1$auto_call), 1, 0))
      geneC3 <- paste(v1[grepl("C3", v1$auto_call),]$gene_symbol, collapse = ";")
      
      
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
                                          "WRN", "XRCC5", "XRCC6", "SMARCAL1", "STAG3"),
                      "mitotic_chkpt" = mito_chk_point$V1,
                      "mirabello_genes" = mira_genes$V1,
                      "centrosome_genes" = centrosome_mat$X__1,
                      "cep_haus_genes" = cep_haus_all$V1)
    gene_list_match <- lapply(gene_list, function(x) ifelse(x %in% v1$gene_symbol, 1, 0))
    sum_list <- unlist(lapply(gene_list_match, function(x)sum(x)))
    sum_list <- as.data.frame(t(sum_list))
    gene_list_C4C5 <- lapply(gene_list, function(x) paste(v1[v1$gene_symbol %in% x & 
                                                         v1$auto_call %in% c("C4","C5","C4_ar","C5_ar"),9], collapse = ";"))
    gene_list_all <- lapply(gene_list, function(x) paste(v1[v1$gene_symbol %in% x & 
                                                               v1$auto_call %in% c("C4","C5","C4_ar","C5_ar", "C3"),9], collapse = ";"))
    C3C4C5_list <- unlist(lapply(gene_list_all, function(x)ifelse(x %in% "", "None", x)))
    C3C4C5_list <- as.data.frame(t(C3C4C5_list))
    
    C4C5_list <- unlist(lapply(gene_list_C4C5, function(x)ifelse(x %in% "", "None", x)))
    C4C5_list <- as.data.frame(t(C4C5_list))
    colnames(C4C5_list) <- paste(colnames(C4C5_list),"C4C5", sep = "_")
    samvar_list[[i]] <- cbind.data.frame("Total.burden" = ntot, "nC5" = nC5, "C5_genes" = geneC5,
                                         "nC4" = nC4, "C4_genes" = geneC4,
                                         "nC3" = nC3, "C3_genes" = geneC3, sum_list, C4C5_list,
                                         "cep_haus_C345" = C3C4C5_list[,13], "SAMPLE" = sam_id[i])
  }
  return(samvar_list)
}

##filtered set
print("Read combined ISKS MGRB input with AD and AR variants")
fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/comb_clin_isks_mgrb_2020_C3C4C5_NFE0002_AD_AR_all_fields_rnd3_Aug31.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)
fil_tab$is_case <- ifelse(grepl("^ISKS", fil_tab$set), 1, 0)
fil_tab_isks_set <- fil_tab[fil_tab$is_case == 1,]
fil_tab_isks_set$auto_call <- ifelse(fil_tab_isks_set$is_AR == 1 & fil_tab_isks_set$GT == 2, paste0(fil_tab_isks_set$auto_call, "_ar"), fil_tab_isks_set$auto_call)
##use this for QC testing not in final (Discovery set topSKAT300)
#top300_SKAT <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Sep05_rect_ASP_graph_final_PPI_comb_GO_cmaf_new_score_str.tsv",
#                          sep = "\t", header = T, stringsAsFactors = F)

##use this in final
##SKATBinary
top_SKAT_cgc <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/DT_skatbin_cgc.txt",
                           sep = "", header = F, stringsAsFactors = F)
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

##Add mito_chk_point
print("Add mito_chk_point genes")
gene_sub <- unique(c(top_SKAT_cgc$V1, cgc_genes$Gene.Symbol, sheltrin_complex, Telo_extension, Sheltrin_comp_extn, 
                     mito_chk_point$V1, mira_genes$V1, centrosome_mat$X__1, cep_haus_all$V1))
comb_set_filt <- fil_tab_isks_set[fil_tab_isks_set$gene_symbol %in% gene_sub,]
comb_set_filt_AR <- comb_set_filt[comb_set_filt$is_AR == 1,]
comb_set_filt_AD <- comb_set_filt[comb_set_filt$is_AR == 0,]
##Add MAF based filter too: DT May_12_2020:
##Run this for AD variants and NOT AR as for SKAT
fil_tab_AD <- fil_tab[fil_tab$is_AR == 0,]
Ex_samp_id <- unique(fil_tab_AD$SAMPLE)

p_vec <- ifelse(!is.na(as.numeric(as.character(Ex_samp_id))) | grepl("^CR|^LK",Ex_samp_id), 1, 0)
#####make genotype matrix function
make_geno_mat <- function(ftemp_file){
  ftemp_tab_var_id <- unique(ftemp_file$VARIANT)
  samp_vec <- list()
  for(m in 1:length(ftemp_tab_var_id)){
    sam_gene_gt <- ftemp_file[ftemp_file$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,82,127:128)]
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
  samp_vec_mat <- samp_vec_mat[!(as.numeric(as.character(samp_vec_mat$coh_MAF_case)) >= 0.00035 | 
                                   as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) >= 0.00035),]
  
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


genes <- unique(comb_set_filt_AD$gene_symbol)
print("save PID genes")
#write.table(genes, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CGC_topSKATBin300_repstress_potint_mito_chkpt_PID_unique_genes_cluego_Aug31.txt",
#            sep = "\t", row.names = F, quote = F)
write.table(genes, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CGC_topSKATBin300_repstress_potint_mito_chkpt_centrosome_PID_unique_genes_cluego_Sep19.txt",
            sep = "\t", row.names = F, quote = F)

print("compute AD variant intra cohort frequencies and filter based on SKAT threshold")
##Parallelised   
library(doParallel)
library(doMC)
registerDoMC(30)

res20 <- list()

system.time(res20 <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
{gene_var_maf(genes[i], comb_set_filt_AD)})
gene_maf_pass_df <- do.call("rbind.data.frame", res20)
comb_set_filt1_AD <- comb_set_filt_AD[comb_set_filt_AD$VARIANT %in% gene_maf_pass_df$VARIANT,]
comb_set_filt1 <- rbind.data.frame(comb_set_filt1_AD, comb_set_filt_AR)

print("save ISKS AD variants")
#write.table(comb_set_filt1, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CGC_topSKATBin300_ISKS_repstress_potint_mito_chkpt_VARIANTS_filt_combset2020_clueGOplus_Aug31.tsv",
#            sep = "\t", row.names = F, quote = F)
write.table(comb_set_filt1, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CGC_topSKATBin300_ISKS_repstress_potint_mito_chkpt_centrosome_VARIANTS_filt_combset2020_clueGOplus_Sep19.tsv",
            sep = "\t", row.names = F, quote = F)
var_per_samp_isks <- summarise_var_per_samp(comb_set_filt1)
df_sam_filt_isks_comb <- do.call("rbind.data.frame", var_per_samp_isks)


##Save all ISKS MGRB variants that are in PID genes (superset of comb_set_filt1)
print("save ISKS MGRB AD variants mapped to PID genes")
fil_tab_PID <- fil_tab[fil_tab$gene_symbol %in% gene_sub,]
fil_tab_PID$auto_call <- ifelse(fil_tab_PID$is_AR == 1 & fil_tab_PID$GT == 2, paste0(fil_tab_PID$auto_call, "_ar"), fil_tab_PID$auto_call)
#write.table(fil_tab_PID, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CGC_topSKATBin_ISKS_MRGB_repstress_potint_mito_chkpt_VARIANTS_filt_combset2020_clueGOplus_Aug31.tsv",
#            sep = "\t", row.names = F, quote = F)
write.table(fil_tab_PID, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CGC_topSKATBin_ISKS_MRGB_repstress_potint_mito_chkpt_centrosome_VARIANTS_filt_combset2020_clueGOplus_Sep19.tsv",
            sep = "\t", row.names = F, quote = F)

#######
##Add phenotypes
##Import telseq, CHIP and clinical information file to VM 
##Also add LION's clinical data
#detach("package:gdata", unload=TRUE)
print("Add clinical correlates")
library(readxl)
##Add telomere length
telseq <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/AllCohorts_telseq_telomeres.xlsx",
                     sheet = 1)
telseq <- as.data.frame(telseq)
##telseq lions
telseq_LK <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/Lions_telseq_67samples.tsv",
                        sep = "\t", header = F, stringsAsFactors = F)
telseq_LK1 <- cbind.data.frame("Cohort" = "lions", telseq_LK) 
colnames(telseq_LK1)[2:3] <- c("Sample", "LENGTH_ESTIMATE")

telseq_all <- rbind.data.frame(telseq, telseq_LK1)
##Chip
isks_chip <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/ISKS_chip1_chip2.chip2_all_categories_nb2_depth_and_vaf_and_position_tallies.tsv",
                        sep = "\t", header = T, stringsAsFactors = F)
risc_chip <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/riscSarcoma_chip1_chip2.chip2_all_categories_nb2_depth_and_vaf_and_position_tallies.tsv",
                        sep = "\t", header = T, stringsAsFactors = F)
risc_chip$sample_id <- paste0("CR", risc_chip$sample_id)
lions_chip <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/lionsSarcoma_chip1_chip2.chip2_all_categories_nb2_depth_and_vaf_and_position_tallies.tsv",
                         sep = "\t", header = T, stringsAsFactors = F)

isks_risc_lions_chip <- rbind.data.frame(isks_chip, risc_chip, lions_chip)
#telseq$Sample <- as.character(telseq$Sample)
##Mandy phenotype
#comb_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/Copy of ISKS_RisC_PID file 170120.xlsx", sheet = 2)
comb_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/PID_Master_file_290420AgeBloodTaken_Aug5.xlsx", sheet = 1, col_types = c("list"))
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

##Additional phenotypes added on Sept.8
add_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/extraPIDSdataset_100920.xlsx", sheet = 1, col_types = c("list"))
add_pheno <- as.data.frame(add_pheno)
add_pheno1 <- sapply(add_pheno, unlist)
colnames(add_pheno1) <- colnames(add_pheno)
add_pheno <- add_pheno1
add_pheno <- as.data.frame(add_pheno, stringsAsFactors = F)
add_pheno <- unique(add_pheno)
add_pheno <- add_pheno[!is.na(add_pheno$pid),]
add_pheno$`age at dateExtracted` <- as.numeric(add_pheno$`age at dateExtracted`)
add_pheno$AgeatSarcoma <- as.numeric(add_pheno$AgeatSarcoma)
add_pheno$SubjectAgeCancer <- as.numeric(add_pheno$SubjectAgeCancer)

colnames(add_pheno) <- colnames(comb_pheno)
##########
##Collate all phenotypes
comb_ALL_phen <- rbind.data.frame(comb_pheno, add_pheno)

comb_ALL_phen <- comb_ALL_phen[comb_ALL_phen$pmn %in% Ex_samp_id,]

##Apr29-2020
##Aug5-2020 new column was added to distinguish between age at blood extracted and age of diagnosis

##Add QC2 fail check
##Doesn't need rectification since Emma was aware of M012 mislabelling and taken care of it
QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/joint_calls_2019jan_2019nov.final_qc_output.tsv", header = T, sep = "\t", stringsAsFactors = F)

QC2_dat_fail <- QC2_dat[QC2_dat$passes_qc2 %in% "FALSE",]$new_sampleid
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %nin% "FALSE",]

`%nin%` = Negate(`%in%`)
comb_ALL_phen$QC2 <- ifelse(as.character(comb_ALL_phen$pmn) %in% QC2_dat_fail, "Fail", 
                         ifelse(as.character(comb_ALL_phen$pmn) %nin% QC2_dat$new_sampleid, "Unknown", "Pass"))


##these step has to be done after running summarise_var_per_samp function 
##Added Sheltrin_extn column; Added Centrosome maturation; Added CEP/HAUS attributes

comb_ALL_phen[,c(68:101)] <- df_sam_filt_isks_comb[match(as.character(comb_ALL_phen$pmn), df_sam_filt_isks_comb$SAMPLE), c(1:34)]

comb_ALL_phen$telomere_length <- telseq_all[match(comb_ALL_phen$pmn, telseq_all$Sample), 3]

#comb_ALL_phen[,c(62:67)] <- isks_risc_lions_chip[match(comb_ALL_phen$pmn, isks_risc_lions_chip$sample_id), c(1:6)]

comb_ALL_phen[,c(103:108)] <- isks_risc_lions_chip[match(comb_ALL_phen$pmn, isks_risc_lions_chip$sample_id), c(1:6)]


####Prob_NFE
# p_Data_comb <- read.table("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_combset_pca.scores_clustered.tsv", 
#                           header = T, sep = "\t", stringsAsFactors = F)
# p_Data_comb <- p_Data_comb[p_Data_comb$superPopulation %in% c("ISKS", "RISC", "LIONS"),]
# p_Data_comb <- p_Data_comb[!duplicated(p_Data_comb$sample),]

p_Data <- read.table("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_combset_pca.scores_clustered.tsv", 
                     header = T, sep = "\t", stringsAsFactors = F)
##rename p_Data since Shyam had not accounted for the error in M012 ISKS in the joint call
p_Data_ISKS <-  p_Data[p_Data$superPopulation %in% c("ISKS", "RISC", "LIONS"),]
rect_sam_dat <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/ISKS_RISC_LIONS_final_freeze.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
p_Data_ISKS$rect_sam <- rect_sam_dat[match(p_Data_ISKS$sample,rect_sam_dat$JCInputID),9]
p_Data_ISKS <- p_Data_ISKS[!is.na(p_Data_ISKS$rect_sam),]
##two samples missed due to mislabelling and QC2
#rect_sam_dat$JCInputRecID[rect_sam_dat$JCInputRecID %nin% p_Data_ISKS$sample]
p_Data_ISKS <- p_Data_ISKS[p_Data_ISKS$rect_sam %in% QC2_dat_pass$new_sampleid,] ##3105 is lost(QC2 fail)
p_Data_comb <- p_Data_ISKS[as.character(p_Data_ISKS$rect_sam) %in% Ex_samp_id,]


##3137 and 3468 are present in replication set and are removed from PC clustering after relatedness check.
table(comb_ALL_phen$pmn %in% p_Data_comb$rect_sam)
comb_ALL_phen$pred.NFE <- p_Data_comb[match(comb_ALL_phen$pmn, p_Data_comb$rect_sam),51]
comb_ALL_phen$prob.NFE <- p_Data_comb[match(comb_ALL_phen$pmn, p_Data_comb$rect_sam),52]

comb_ALL_phen_filt <- comb_ALL_phen[comb_ALL_phen$QC2 %in% "Pass" | comb_ALL_phen$QC2 %in% "Unknown",]
# ##mark duplicate samples
# deg2_isks <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/dup_pairs_and_deg2_isks.tsv", sep = "\t",
#                         header = T, stringsAsFactors = F)
# 
# deg2_isks <- deg2_isks[deg2_isks$degree == 0,]
# 
# comb_ALL_phen_filt$duplicate <- deg2_isks[match(comb_ALL_phen_filt$pmn, deg2_isks$name2), 7]

comb_ALL_phen_filt[is.na(comb_ALL_phen_filt$Total.burden),c(68:69,71,73,75:87)] <- 0

print("Save PID file")
#write.table(comb_ALL_phen_filt[,-c(75:81)], 
#            "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/PID_combset2020_CGC_skatBin_repstress_potint_mito_chkpt_predNFE_clueGO_Sep102020_AD_Aug31.tsv",
#            sep = "\t", row.names = F, quote = F)
#write.table(comb_ALL_phen_filt[,-c(75:81,86:92)], 
#            "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/PID_combset2020_CGC_skatBin_repstress_potint_mito_chkpt_predNFE_clueGO_Sep102020_AD_Aug31_addC4C5.tsv",
#            sep = "\t", row.names = F, quote = F)
write.table(comb_ALL_phen_filt[,-c(75:81,88:94)], 
            "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/PID_combset2020_CGC_skatBin_repstress_potint_mito_chkpt_centrosome_predNFE_clueGO_Sep192020_AD_addC4C5.tsv",
            sep = "\t", row.names = F, quote = F)
