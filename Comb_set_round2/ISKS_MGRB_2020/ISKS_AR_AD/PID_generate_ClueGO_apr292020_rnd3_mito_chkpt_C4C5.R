##comb_set_filt1 for mito_variants with C4 and C5 variants only

.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

`%nin%` = Negate(`%in%`)

mito_chk_point = read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/mitotic_checkpoint/mitocheck_point.txt",
                            sep = "", header = F, skip = 1, stringsAsFactors = F)

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
                                          "WRN", "XRCC5", "XRCC6"),
                      "mitotic_chkpt" = mito_chk_point$V1)
    gene_list_match <- lapply(gene_list, function(x) ifelse(x %in% v1$gene_symbol, 1, 0))
    sum_list <- unlist(lapply(gene_list_match, function(x)sum(x)))
    sum_list <- as.data.frame(t(sum_list))
    samvar_list[[i]] <- cbind.data.frame("Total.burden" = ntot, "nC5" = nC5, "C5_genes" = geneC5,
                                         "nC4" = nC4, "C4_genes" = geneC4,
                                         "nC3" = nC3, "C3_genes" = geneC3, sum_list, "SAMPLE" = sam_id[i])
  }
  return(samvar_list)
}

t1 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/PID_combset2020_CGC_skatBin_repstress_potint_mito_chkpt_predNFE_clueGO_Apr292020_AD_Aug19.tsv", sep = "\t", header = T, stringsAsFactors = F)
t1 <- comb_ALL_phen_filt[,-c(75:81)]
t2 <- t1[t1$mitotic_chkpt > 0,]
t2 <- t2[!is.na(t2$pmn),]
comb_set_filt1 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CGC_topSKATBin300_repstress_potint_mito_chkpt_VARIANTS_filt_combset2020_clueGOplus_Aug19.tsv",
                             sep = "\t", header = T, stringsAsFactors = F)
comb_set_filt1_var <- comb_set_filt1[comb_set_filt1$SAMPLE %in% t2$pmn,]
comb_set_filt1_var <- comb_set_filt1_var[comb_set_filt1_var$auto_call %in% c("C4", "C5", "C5_ar", "C4_ar"),]

var_per_samp_isks_mito <- summarise_var_per_samp(comb_set_filt1_var)
df_sam_filt_isks_comb_mito <- do.call("rbind.data.frame", var_per_samp_isks_mito)
df_sam_filt_isks_comb_mito <- df_sam_filt_isks_comb_mito[df_sam_filt_isks_comb_mito$mitotic_chkpt > 0,]

#######
##Add phenotypes
##Import telseq, CHIP and clinical information file to VM 
##Also add LION's clinical data
#detach("package:gdata", unload=TRUE)
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


##Apr29-2020
##Aug5-2020 new column was added to distinguish between age at blood extracted and age of diagnosis

##Add QC2 fail check
QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/joint_calls_2019jan_2019nov.final_qc_output.tsv", header = T, sep = "\t", stringsAsFactors = F)

QC2_dat_fail <- QC2_dat[QC2_dat$passes_qc2 %in% "FALSE",]$new_sampleid

`%nin%` = Negate(`%in%`)
comb_pheno$QC2 <- ifelse(as.character(comb_pheno$pmn) %in% QC2_dat_fail, "Fail", 
                         ifelse(as.character(comb_pheno$pmn) %nin% QC2_dat$new_sampleid, "Unknown", "Pass"))

##########
##Collate all phenotypes

#comb_ALL_phen <- rbind.data.frame(comb_pheno, lion_pheno1)
comb_ALL_phen <- comb_pheno
##these step has to be done after running summarise_var_per_samp function 
##Added Sheltrin_extn column

comb_ALL_phen[,c(68:84)] <- df_sam_filt_isks_comb_mito[match(as.character(comb_ALL_phen$pmn), df_sam_filt_isks_comb_mito$SAMPLE), c(1:17)]

comb_ALL_phen$telomere_length <- telseq_all[match(comb_ALL_phen$pmn, telseq_all$Sample), 3]

#comb_ALL_phen[,c(62:67)] <- isks_risc_lions_chip[match(comb_ALL_phen$pmn, isks_risc_lions_chip$sample_id), c(1:6)]

comb_ALL_phen[,c(86:91)] <- isks_risc_lions_chip[match(comb_ALL_phen$pmn, isks_risc_lions_chip$sample_id), c(1:6)]


####Prob_NFE
p_Data_comb <- read.table("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_combset_pca.scores_clustered.tsv", 
                          header = T, sep = "\t", stringsAsFactors = F)
p_Data_comb <- p_Data_comb[p_Data_comb$superPopulation %in% c("ISKS", "RISC", "LIONS"),]
p_Data_comb <- p_Data_comb[!duplicated(p_Data_comb$sample),]


##3137 and 3468 are present in replication set and are removed from PC clustering after relatedness check.
table(comb_ALL_phen$pmn %in% p_Data_comb$sample)
comb_ALL_phen$pred.NFE <- p_Data_comb[match(comb_ALL_phen$pmn, p_Data_comb$sample),51]
comb_ALL_phen$prob.NFE <- p_Data_comb[match(comb_ALL_phen$pmn, p_Data_comb$sample),52]

comb_ALL_phen_filt <- comb_ALL_phen[comb_ALL_phen$QC2 %in% "Pass" | comb_ALL_phen$QC2 %in% "Unknown",]
##mark duplicate samples
deg2_isks <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/dup_pairs_and_deg2_isks.tsv", sep = "\t",
                        header = T, stringsAsFactors = F)

deg2_isks <- deg2_isks[deg2_isks$degree == 0,]

comb_ALL_phen_filt$duplicate <- deg2_isks[match(comb_ALL_phen_filt$pmn, deg2_isks$name2), 7]
comb_ALL_phen_filt <- comb_ALL_phen_filt[comb_ALL_phen_filt$mitotic_chkpt > 0,]
comb_ALL_phen_filt <- comb_ALL_phen_filt[!is.na(comb_ALL_phen_filt$mitotic_chkpt),]
write.table(comb_ALL_phen_filt[,-c(75:81)], 
            "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/PID_combset2020_C4C5_mito_chkpt_predNFE_clueGO_Apr292020_AD_Aug19.tsv",
            sep = "\t", row.names = F, quote = F)
