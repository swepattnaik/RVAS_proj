##Make PID file
##with DT criteria: VAF >= 0.35 & Eigenphred >= 5.6
##PID file: Jan17-2020
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
fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/all_isksmgrb_combset2020_variants_filt_all_fields_rmdup.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)
fil_tab_isks_set <- fil_tab[fil_tab$is_case == 1,]
##use this for QC testing not in final (Discovery set topSKAT300)
#top300_SKAT <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Sep05_rect_ASP_graph_final_PPI_comb_GO_cmaf_new_score_str.tsv",
#                          sep = "\t", header = T, stringsAsFactors = F)

##use this in final
##SKATO
#top300_SKAT <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/ppi_res_fil_final.tsv",
#                          sep = "\t", header = T, stringsAsFactors = F)
##SKATBinary
top300_SKAT <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/ppi_res_fil_final_SKATbin.tsv",
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
comb_set_filt <- fil_tab_isks_set[fil_tab_isks_set$gene_symbol %in% gene_sub,]
comb_set_filt$gnomad_AF <- gsub("\\[|\\]", "", comb_set_filt$gnomad_AF)
comb_set_filt <- comb_set_filt[comb_set_filt$VAF >= 0.35 & comb_set_filt$comb_score >= 5.6,]

write.table(comb_set_filt, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/CGC_topSKATBin300_repstress_potint_VARIANTS_filt_combset2020_rmdup.tsv",
            sep = "\t", row.names = F, quote = F)
var_per_samp_isks <- summarise_var_per_samp(comb_set_filt)
df_sam_filt_isks_comb <- do.call("rbind.data.frame", var_per_samp_isks)

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
comb_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/Copy of ISKS_RisC_PID file 170120.xlsx", sheet = 2)
comb_pheno <- as.data.frame(comb_pheno, stringsAsFactors = F)
comb_pheno <- unique(comb_pheno)
comb_pheno <- comb_pheno[,-c(44:55)]
comb_pheno$pid <- as.character(comb_pheno$pid)
comb_pheno$`age at dateExtracted` <- as.numeric(comb_pheno$`age at dateExtracted`)
comb_pheno$AgeatSarcoma <- as.numeric(comb_pheno$AgeatSarcoma)
comb_pheno$SubjectAgeCancer <- as.numeric(comb_pheno$SubjectAgeCancer)


##rectify pmns based on pids: Mandy's input Jan17-2020

comb_pheno$pmn <- ifelse(grepl("^901.|^902.",comb_pheno$pid), paste0("CR", comb_pheno$pmn), comb_pheno$pmn)

##Add QC2 fail check
QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/joint_calls_2019jan_2019nov.final_qc_output.tsv", header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_fail <- QC2_dat[QC2_dat$passes_qc2 %in% "FALSE",]$new_sampleid

`%nin%` = Negate(`%in%`)
comb_pheno$QC2 <- ifelse(as.character(comb_pheno$pmn) %in% QC2_dat_fail, "Fail", 
                         ifelse(as.character(comb_pheno$pmn) %nin% QC2_dat$new_sampleid, "Unknown", "Pass"))

########David's LION phenotype
lion_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/Lions_ZCC_germlines_M084.xlsx",
                         sheet = 1)
lion_pheno <- as.data.frame(lion_pheno)

lion_pheno1 <- lion_pheno

##For Jan17
lion_pheno1[,c(12:55)] <- colnames(comb_pheno)
colnames(lion_pheno1)[12:55] <- colnames(comb_pheno)
lion_pheno1[,c(12:55)] <- NA

lion_pheno1$pid <- lion_pheno1$matched_normal_id
lion_pheno1$pmn <- lion_pheno1$matched_normal_id
lion_pheno1$gender1 <- ifelse(lion_pheno1$sex %in% "Female", "F", "M")
lion_pheno1$CaseControl <- "Proband"
lion_pheno1$`age at dateExtracted` <- lion_pheno1$age_at_sample
lion_pheno1$SubjectCancers <- lion_pheno1$final_diagnosis
lion_pheno1$SarcomaTopography <- lion_pheno1$cancer_category
lion_pheno1$sarcomatype <- lion_pheno1$cancer_type

lion_pheno1$QC2 <- ifelse(as.character(lion_pheno1$pmn) %in% QC2_dat_fail, "Fail", "Pass")
lion_pheno1 <- lion_pheno1[,-c(1:11)]
lion_pheno1 <- lion_pheno1[lion_pheno1$pmn %in% QC2_dat$new_sampleid,]
##########
##Collate all phenotypes

comb_ALL_phen <- rbind.data.frame(comb_pheno, lion_pheno1)
##these step has to be done after running summarise_var_per_samp function 
##Added Sheltrin_extn column
comb_ALL_phen[,c(45:60)] <- df_sam_filt_isks_comb[match(as.character(comb_ALL_phen$pmn), df_sam_filt_isks_comb$SAMPLE), c(1:16)]

comb_ALL_phen$telomere_length <- telseq_all[match(comb_ALL_phen$pmn, telseq_all$Sample), 3]

comb_ALL_phen[,c(62:67)] <- isks_risc_lions_chip[match(comb_ALL_phen$pmn, isks_risc_lions_chip$sample_id), c(1:6)]



##Variant file merged with phenotypes(not necessary, but generate later if needed)
#comb_set_filt[,c(133:198)] <- comb_ALL_phen[match(comb_ALL_phen$SAMPLE, comb_ALL_phen$pmn),c(1:66)]
#comb_set_filt_merged <- comb_set_filt[,-c(178:191)]
#write.table(var_set_merged, 
#            "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/Variants_CGC_topskat_repstress_genes_combet2020_clin_all_unfilt.tsv",
#            sep = "\t", row.names = F, quote = F)

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

#comb_ALL_phen_filt$duplicate <- ifelse(comb_ALL_phen_filt$pmn %in% deg2_isks$name1, deg2_isks$name2,
#                                       ifelse(comb_ALL_phen_filt$pmn %in% deg2_isks$name2, deg2_isks$name1, "none"))

#write PID file

# write.table(comb_ALL_phen_filt[,-c(52:58)], 
#             "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/PID_combset2020_CGC_skatBin_repstress_potint_predNFE_rmdup.tsv",
#             sep = "\t", row.names = F, quote = F)
write.table(comb_ALL_phen_filt[,-c(52:58)], 
            "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/PID_combset2020_CGC_skatBin_repstress_potint_predNFE_rmdup_Apr292020.tsv",
            sep = "\t", row.names = F, quote = F)
