##Augment Replication set with different proportions of ISKS data (n = 100, 200, 400) 
##make table for ISKS and Repset for each condition 
##compute ORs to assess the relative enrichment in separate complexes

.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

`%nin%` = Negate(`%in%`)
##ISKS+MGRB skat input
Ex_tab_isks_mgrb <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv",
                               sep = "\t", header = T, stringsAsFactors = F)
##ISKS
Ex_tab_isks_disc <- Ex_tab_isks_mgrb[Ex_tab_isks_mgrb$set %in% "ISKS_AR_AD",]
##cohort MAF filter (~ 3/(2*1663): joint called number)
Ex_tab_isks_disc <- Ex_tab_isks_disc[Ex_tab_isks_disc$cohort_MAF <= 0.001,]

isks_filter <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/isks_int_var.bed",
                          sep = "\t", header = F, stringsAsFactors = F)
Ex_tab_isks <- Ex_tab_isks_disc[Ex_tab_isks_disc$VARIANT %in% isks_filter$V4,]

Ex_samp_id <- unique(Ex_tab_isks$SAMPLE)
set.seed(3782921)
Ex_samp_id_100 <- Ex_samp_id[sample(1:length(Ex_samp_id), size = 100, replace = F)]
set.seed(37829234)
Ex_samp_id_200 <- Ex_samp_id[sample(1:length(Ex_samp_id), size = 200, replace = F)]
set.seed(43782921)
Ex_samp_id_400 <- Ex_samp_id[sample(1:length(Ex_samp_id), size = 400, replace = F)]

Ex_tab_isks_100 <- Ex_tab_isks[Ex_tab_isks$SAMPLE %in% Ex_samp_id_100,]
Ex_tab_isks_200 <- Ex_tab_isks[Ex_tab_isks$SAMPLE %in% Ex_samp_id_200,]
Ex_tab_isks_400 <- Ex_tab_isks[Ex_tab_isks$SAMPLE %in% Ex_samp_id_400,]

##MGRB
Ex_tab_mgrb_disc <- Ex_tab_isks_mgrb[Ex_tab_isks_mgrb$set %in% "MGRB_AR_AD",]
##cohort MAF filter (~ 3/(2*3209): joint called number)
Ex_tab_mgrb_disc <- Ex_tab_mgrb_disc[Ex_tab_mgrb_disc$cohort_MAF <= 0.0005,]

mgrb_filter <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/MGRB/mgrb_int_var.bed",
                          sep = "\t", header = F, stringsAsFactors = F)
Ex_tab_mgrb <- Ex_tab_mgrb_disc[Ex_tab_mgrb_disc$VARIANT %in% mgrb_filter$V4,]

##TCGA(253)
Ex_tab_sarc <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/TCGA_SARC_AR_AD/all_tcga_sarc_combset2020_variants_AR_AD_all_fields_clingene_rnd3.tsv", sep = "\t", header = T, stringsAsFactors = F)
sarc_filter <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/TCGA_SARC_AR_AD/Sarc_tcga_int_var.bed",
                          sep = "\t", header = F, stringsAsFactors = F)
Ex_tab_sarc <- Ex_tab_sarc[Ex_tab_sarc$VARIANT %in% sarc_filter$V4,]

##cohort MAF filter (4/(2*253))
Ex_tab_sarc <- Ex_tab_sarc[Ex_tab_sarc$cohort_MAF <= 0.008,]

##BERG(96)
Ex_tab_BERG <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/SARC_BERG_AR_AD/all_sarc_berg_combset2020_variants_AR_AD_all_fields_clingene_rnd3.tsv", sep = "\t", header = T, stringsAsFactors = F)
berg_filter <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/SARC_BERG_AR_AD/sarc_berg_int_var.bed",
                          sep = "\t", header = F, stringsAsFactors = F)
Ex_tab_BERG <- Ex_tab_BERG[Ex_tab_BERG$VARIANT %in% berg_filter$V4,]
##cohort MAF filter (4/(2*96))
Ex_tab_BERG <- Ex_tab_BERG[Ex_tab_BERG$cohort_MAF <= 0.021,]


##NOSARC(214)
Ex_tab_nosarc <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/NoSarc_nov2018_AR_AD/all_nosarc_nov18_combset2020_variants_AR_AD_all_fields_clingene_rnd3.tsv", sep = "\t", header = T, stringsAsFactors = F)
nosarc_filter <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/NoSarc_nov2018_AR_AD/NoSarc_Nov2018_int_var.bed",
                            sep = "\t", header = F, stringsAsFactors = F)
Ex_tab_nosarc <- Ex_tab_nosarc[Ex_tab_nosarc$VARIANT %in% nosarc_filter$V4,]
##cohort MAF filter (4/(2*214))
Ex_tab_nosarc <- Ex_tab_nosarc[Ex_tab_nosarc$cohort_MAF <= 0.010,]

Ex_tab <- rbind.data.frame(Ex_tab_sarc,Ex_tab_BERG, Ex_tab_nosarc)
Ex_tab$gnomad_AF <- gsub("\\[", "", Ex_tab$gnomad_AF)
Ex_tab$gnomad_AF <- gsub("\\]", "", Ex_tab$gnomad_AF)
class(Ex_tab$gnomad_AF) <- "numeric"
Ex_tab$gnomad_AF_NFE <- gsub("\\[", "", Ex_tab$gnomad_AF_NFE)
Ex_tab$gnomad_AF_NFE <- gsub("\\]", "", Ex_tab$gnomad_AF_NFE)
class(Ex_tab$gnomad_AF_NFE) <- "numeric"
Ex_tab$gnomad_AF <- ifelse(is.na(Ex_tab$gnomad_AF), 0, Ex_tab$gnomad_AF)
Ex_tab$swegen_AF <- ifelse(is.na(Ex_tab$swegen_AF), 0, Ex_tab$swegen_AF)
Ex_tab$gnomad_AF_NFE <- ifelse(is.na(Ex_tab$gnomad_AF_NFE), 0, Ex_tab$gnomad_AF_NFE)
##GQ filter

Ex_tab <- Ex_tab[Ex_tab$GQ >= 80,]

##AD
##C4/C5; gnomad_AF_NFE <=  0.0002
#Ex_tab_C4C5 <- Ex_tab[Ex_tab$auto_call %in% c("C4", "C5") & Ex_tab$gnomad_AF < 0.01,]
Ex_tab_C4C5 <- Ex_tab[Ex_tab$auto_call %in% c("C4", "C5") & Ex_tab$gnomad_AF_NFE <= 0.0002,]
#Ex_tab_C4C5$cohort_MAF <- ifelse(Ex_tab_C4C5$set %in% "TCGA_SARC_AR_AD", (2*253*Ex_tab_C4C5$cohort_MAF)/(2*563),
#                                 ifelse(Ex_tab_C4C5$set %in% "SARC_BERG_AR_AD", (2*96*Ex_tab_C4C5$cohort_MAF)/(2*563),
#                                        ifelse(Ex_tab_C4C5$set %in% "NOSARC_nov18_AR_AD", (2*214*Ex_tab_C4C5$cohort_MAF)/(2*563), Ex_tab_C4C5$cohort_MAF)))
##Add combined intra_cohort MAF filter cohort_MAF <= 0.01 for combined SARC & 0.0008 for  MGRB ##Added later
Ex_tab_C4C5_sarc <- Ex_tab_C4C5[Ex_tab_C4C5$swegen_AF <= 0.001,]
#Ex_tab_C4C5_sarc <- Ex_tab_C4C5[Ex_tab_C4C5$cohort_MAF <= 0.005 & Ex_tab_C4C5$swegen_AF <= 0.001,]
#Ex_tab_C4C5_sarc <- Ex_tab_C4C5[Ex_tab_C4C5$cohort_MAF <= 0.005 & Ex_tab_C4C5$swegen_AF <= 0.001 & Ex_tab_C4C5$set %nin% "MGRB_AR_AD",]
#Ex_tab_C4C5_mrgb <- Ex_tab_C4C5[Ex_tab_C4C5$cohort_MAF <= 0.0008 & Ex_tab_C4C5$swegen_AF <= 0.001 & Ex_tab_C4C5$set %in% "MGRB_AR_AD",]
#Ex_tab_C4C5 <- rbind.data.frame(Ex_tab_C4C5_sarc, Ex_tab_C4C5_mrgb)
##C3; gnomad_AF == 0
#Ex_tab_C3 <- Ex_tab[Ex_tab$auto_call %in% "C3" & Ex_tab$gnomad_AF <= 0.01 & Ex_tab$gnomad_AF_NFE == 0,]
Ex_tab_C3 <- Ex_tab[Ex_tab$auto_call %in% "C3" & Ex_tab$gnomad_AF_NFE  == 0 & Ex_tab$swegen_AF == 0,]
##cohort_MAF <= 0.001
# Ex_tab_C3$cohort_MAF <- ifelse(Ex_tab_C3$set %in% "TCGA_SARC_AR_AD", (2*253*Ex_tab_C3$cohort_MAF)/(2*563),
#                                ifelse(Ex_tab_C3$set %in% "SARC_BERG_AR_AD", (2*96*Ex_tab_C3$cohort_MAF)/(2*563),
#                                       ifelse(Ex_tab_C3$set %in% "NOSARC_nov18_AR_AD", (2*214*Ex_tab_C3$cohort_MAF)/(2*563), Ex_tab_C3$cohort_MAF)))
# Ex_tab_C3 <- Ex_tab_C3[Ex_tab_C3$cohort_MAF <= 0.005 & Ex_tab_C3$swegen_AF <= 0.001,]
#Ex_tab_C3_sarc <- Ex_tab_C3[Ex_tab_C3$cohort_MAF <= 0.005 & Ex_tab_C3$swegen_AF <= 0.001 & Ex_tab_C3$set %nin% "MGRB_AR_AD",]
#Ex_tab_C3_mrgb <- Ex_tab_C3[Ex_tab_C3$cohort_MAF <= 0.0008 & Ex_tab_C3$swegen_AF <= 0.001 & Ex_tab_C3$set %in% "MGRB_AR_AD",]
#Ex_tab_C3 <- rbind.data.frame(Ex_tab_C3_sarc, Ex_tab_C3_mrgb)

Ex_tab_C3 <- Ex_tab_C3[Ex_tab_C3$comb_score >= 5.6,]

##Compound Het Filter(FS, 5bp) for Ex_tab_C4C5
all_genes <- unique(Ex_tab_C4C5$gene_symbol)
sam_sel_AD <- list()
for(i in 1:length(all_genes)){
  # gene_df <- Ex_tab[Ex_tab$gene_symbol %in% all_genes[i],]
  gene_df <- Ex_tab_C4C5[Ex_tab_C4C5$gene_symbol %in% all_genes[i],]
  if(dim(gene_df)[1] > 1){
    sam_freq <- as.data.frame(table(gene_df$SAMPLE))
    sam_sel <- as.character(sam_freq[sam_freq$Freq > 1,1])
    if(length(sam_sel) > 0){
      sam_sel_AD[[i]] <- gene_df[gene_df$SAMPLE %in% sam_sel,]
    }
    else{sam_sel_AD[[i]] <- NULL}
  }
  else { next }
  
}

AD_cpd_het_C45 <- do.call("rbind.data.frame", sam_sel_AD)
AD_cpd_het_C45 <- unique(AD_cpd_het_C45)

var_coods_all <- lapply(strsplit(AD_cpd_het_C45$VARIANT, ":"),
                        function(x)cbind.data.frame("chr" = x[1], "pos" = x[2]))

var_coods_df_all <-  do.call("rbind.data.frame", var_coods_all) 
var_coods_df_all$pos <- as.numeric(as.character(var_coods_df_all$pos))
AD_cpd_het_pos_C45 <- cbind.data.frame(AD_cpd_het_C45, var_coods_df_all)
#rm(var_coods_all,var_coods_df_all)

all_genes <- unique(AD_cpd_het_pos_C45$gene_symbol)
#extract hits for removal; all variants within 5 bps of the fs variant
sam_sel_ns_rm <- list()
sam_sel_ns_rm_ls <- list()
for(i in 1:length(all_genes)){
  gene_df <- AD_cpd_het_pos_C45[AD_cpd_het_pos_C45$gene_symbol %in% all_genes[i],]
  sam_freq <- as.data.frame(table(gene_df$SAMPLE))
  sam_sel <- as.character(sam_freq[sam_freq$Freq >1,1])
  if(length(sam_sel) > 0){
    for(k in 1:length(sam_sel)){
      sam_sel_df <- gene_df[gene_df$SAMPLE %in% sam_sel[k],]
      #  pos_fs <- sam_sel_df[sam_sel_df$vep_consequence %in% "frameshift_variant",]$pos
      pos_fs <- sam_sel_df[grepl("frameshift_variant", sam_sel_df$vep_consequence),]$pos
      #   pos_nfs <- sam_sel_df[sam_sel_df$vep_consequence %nin% "frameshift_variant",]$pos
      pos_nfs <- sam_sel_df[sam_sel_df$pos != pos_fs,]$pos
      sel_var <- ifelse(abs(pos_fs - pos_nfs) <= 5, 1, 0) 
      pos_nfs <- pos_nfs[which(sel_var == 1)]
      sam_sel_ns_rm[[k]] <-  sam_sel_df[sam_sel_df$pos %in% pos_nfs,]
    }
    sam_sel_ns_rm_ls[[i]] <- do.call("rbind.data.frame",sam_sel_ns_rm)
  }
  else{sam_sel_ns_rm_ls[[i]] <- NULL}
}
sam_rm_het_df <- do.call("rbind.data.frame",sam_sel_ns_rm_ls)
sam_rm_het_df <- unique(sam_rm_het_df)

##Ex_tab compound het filter: removes all variants within 5 bps of the FS variant
AD_cpd_het_pos_C45 <- AD_cpd_het_pos_C45[AD_cpd_het_pos_C45$filt_tab %nin% sam_rm_het_df$filt_tab,]

Ex_tab_C4C5 <- Ex_tab_C4C5[Ex_tab_C4C5$filt_tab %nin% sam_rm_het_df$filt_tab,]

##For SKAT(use the same set of variants used above, compound het filtered, gnomad_AF_NFE filtered)
Ex_tab_C3C4C5_AD_filt <- unique(rbind.data.frame(Ex_tab_C4C5, Ex_tab_C3))
rm_X_GT2_AD <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$GT == 2 & grepl("^X", Ex_tab_C3C4C5_AD_filt$VARIANT), ]
Ex_tab_C3C4C5_AD_filt <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$filt_tab %nin% rm_X_GT2_AD$filt_tab,]
Ex_tab_C3C4C5_AD_filt <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$comb_score >= 5.6,]
Ex_tab_C3C4C5_AD_filt_100 <- rbind.data.frame(Ex_tab_C3C4C5_AD_filt, Ex_tab_isks_100, Ex_tab_mgrb)
Ex_tab_C3C4C5_AD_filt_200 <- rbind.data.frame(Ex_tab_C3C4C5_AD_filt, Ex_tab_isks_200, Ex_tab_mgrb)
Ex_tab_C3C4C5_AD_filt_400 <- rbind.data.frame(Ex_tab_C3C4C5_AD_filt, Ex_tab_isks_400, Ex_tab_mgrb)
write.table(Ex_tab_C3C4C5_AD_filt_100, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_tcga_nosarc_berg_isks100_mgrb.tsv", 
            row.names = F, quote = F, sep = "\t")
write.table(Ex_tab_C3C4C5_AD_filt_200, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_tcga_nosarc_berg_isks200_mgrb.tsv", 
            row.names = F, quote = F, sep = "\t")
write.table(Ex_tab_C3C4C5_AD_filt_400, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_tcga_nosarc_berg_isks400_mgrb.tsv", 
            row.names = F, quote = F, sep = "\t")

##save isks without the augmented samples(only for PID genes)
##test in discovery set
Ex_tab_isks_minus_100 <- Ex_tab_isks_disc[Ex_tab_isks_disc$SAMPLE %nin% Ex_samp_id_100,]
Ex_tab_isks_minus_100_mgrb <- rbind.data.frame(Ex_tab_isks_minus_100, Ex_tab_mgrb_disc)
Ex_tab_isks_minus_200 <- Ex_tab_isks_disc[Ex_tab_isks_disc$SAMPLE %nin% Ex_samp_id_200,]
Ex_tab_isks_minus_200_mgrb <- rbind.data.frame(Ex_tab_isks_minus_200, Ex_tab_mgrb_disc)
Ex_tab_isks_minus_400 <- Ex_tab_isks_disc[Ex_tab_isks_disc$SAMPLE %nin% Ex_samp_id_400,]
Ex_tab_isks_minus_400_mgrb <- rbind.data.frame(Ex_tab_isks_minus_400, Ex_tab_mgrb_disc)
Ex_tab_isks_mgrb_all <- rbind.data.frame(Ex_tab_isks_disc, Ex_tab_mgrb_disc)
##AD PID variants for DT
#PIDgenes
top_SKAT_cgc <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/DT_skatbin_cgc.txt",
                           sep = "", header = F, stringsAsFactors = F)
cgc_genes <- read.delim("~/RVAS/cancer_gene_census_hg37.csv", sep = ",", header = T, stringsAsFactors = F)

sheltrin_complex <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "ACD")
Telo_extension <- c("TIMELESS", "TIPIN", "FANCM", "BRCA1", "BLM") ##replication stress is a source of telomere recombination
Sheltrin_comp_extn = c("ACD", "POT1", "TERF1", "TERF2", "TERF2IP", "TINF2", "ATM",
                       "BAG3", "BLM", "BRCA1", "CALD1", "CLK3", "DCLRE1B", "FANCD2",
                       "FBXO4", "HSPA4", "KIAA1191", "MRE11A", "NBN", "PINX1", "PRKDC",
                       "RAD50", "SLX4", "STUB1", "TNKS", "TNKS2", "U2AF2", "UCHL1",
                       "WRN", "XRCC5", "XRCC6", "SMARCAL1", "STAG3")
gene_sub <- unique(c(top_SKAT_cgc$V1, cgc_genes$Gene.Symbol, sheltrin_complex, Telo_extension, Sheltrin_comp_extn))

Ex_tab_C3C4C5_AD_filt_100_PID <- Ex_tab_C3C4C5_AD_filt_100[Ex_tab_C3C4C5_AD_filt_100$gene_symbol %in% gene_sub,]
Ex_tab_C3C4C5_AD_filt_200_PID <- Ex_tab_C3C4C5_AD_filt_200[Ex_tab_C3C4C5_AD_filt_200$gene_symbol %in% gene_sub,]
Ex_tab_C3C4C5_AD_filt_400_PID <- Ex_tab_C3C4C5_AD_filt_400[Ex_tab_C3C4C5_AD_filt_400$gene_symbol %in% gene_sub,]

Ex_tab_isks_mgrb_all_PID <- Ex_tab_isks_mgrb_all[Ex_tab_isks_mgrb_all$gene_symbol %in% gene_sub, ]
Ex_tab_isks_minus_100_mgrb_PID <- Ex_tab_isks_minus_100_mgrb[Ex_tab_isks_minus_100_mgrb$gene_symbol %in% gene_sub,]
Ex_tab_isks_minus_200_mgrb_PID <- Ex_tab_isks_minus_200_mgrb[Ex_tab_isks_minus_200_mgrb$gene_symbol %in% gene_sub,]
Ex_tab_isks_minus_400_mgrb_PID <- Ex_tab_isks_minus_400_mgrb[Ex_tab_isks_minus_400_mgrb$gene_symbol %in% gene_sub,]

write.table(Ex_tab_C3C4C5_AD_filt_100_PID, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_tcga_nosarc_berg_isks100_mgrb_PID.tsv", 
            row.names = F, quote = F, sep = "\t")
write.table(Ex_tab_C3C4C5_AD_filt_200_PID, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_tcga_nosarc_berg_isks200_mgrb_PID.tsv", 
            row.names = F, quote = F, sep = "\t")
write.table(Ex_tab_C3C4C5_AD_filt_400_PID, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_tcga_nosarc_berg_isks400_mgrb_PID.tsv", 
            row.names = F, quote = F, sep = "\t")

write.table(Ex_tab_isks_minus_100_mgrb_PID, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_isks_mgrb_minus100.tsv", 
            row.names = F, quote = F, sep = "\t")
write.table(Ex_tab_isks_minus_200_mgrb_PID, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_isks_mgrb_minus200.tsv", 
            row.names = F, quote = F, sep = "\t")
write.table(Ex_tab_isks_minus_400_mgrb_PID, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_isks_mgrb_minus400.tsv", 
            row.names = F, quote = F, sep = "\t")
write.table(Ex_tab_isks_mgrb_all_PID, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/all_isks_mgrb_maf_filt.tsv", 
            row.names = F, quote = F, sep = "\t")
