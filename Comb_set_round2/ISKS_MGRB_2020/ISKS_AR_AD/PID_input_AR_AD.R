###combine relaxed gnomad_AF ISKS and MGRB
##AD will be used for SKAT
##Filtering flowchart for both cohorts
##1. gnomad_AF < 0.05, 
##2. AD : For C4/C5; gnomad_AF < 0.01 & For C3; gnomad_AF == 0
##3. AD : intra_cohort MAF filter cohort_MAF < 0.001 
##4. AD : swegen < 0.001 and gnomad_NFE_AF < 0.0002 | gnomad_NFE_AF < 0.0005
##New filtering criteria: June 12 2020
##HAIL output (gnomad_AF_NFE <= 0.01)
## VAF >= 0.35 & DP >= 10 & GQ > 80

##AD
#C4/5

#C3

##AR
#C4/5

##Comp het filter within 5 bps of the best FS.

.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

`%nin%` = Negate(`%in%`)

#PIDgenes
top_SKAT_cgc <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/DT_skatbin_cgc.txt",
                           sep = "", header = F, stringsAsFactors = F)
cgc_genes <- read.delim("~/RVAS/cancer_gene_census_hg37.csv", sep = ",", header = T, stringsAsFactors = F)
#remove fusion genes from cgc list
#cgc_genes$Gene.Symbol[(grep("fusion", cgc_genes$Role.in.Cancer))]
sheltrin_complex <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "ACD")
Telo_extension <- c("TIMELESS", "TIPIN", "FANCM", "BRCA1", "BLM") ##replication stress is a source of telomere recombination
Sheltrin_comp_extn = c("ACD", "POT1", "TERF1", "TERF2", "TERF2IP", "TINF2", "ATM", 
                       "BAG3", "BLM", "BRCA1", "CALD1", "CLK3", "DCLRE1B", "FANCD2", 
                       "FBXO4", "HSPA4", "KIAA1191", "MRE11A", "NBN", "PINX1", "PRKDC", 
                       "RAD50", "SLX4", "STUB1", "TNKS", "TNKS2", "U2AF2", "UCHL1", 
                       "WRN", "XRCC5", "XRCC6")
gene_sub <- unique(c(top_SKAT_cgc$V1, cgc_genes$Gene.Symbol, sheltrin_complex, Telo_extension, Sheltrin_comp_extn))

Ex_tab_mgrb <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/MGRB/all_mgrb_combset2020_variants_AR_AD_all_fields_clingene_rnd3.tsv", sep = "\t", header = T, stringsAsFactors = F)
Ex_tab_isks <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/all_isks_combset2020_variants_AR_AD_all_fields_clingene_rnd3.tsv", sep = "\t", header = T, stringsAsFactors = F)
Ex_tab_mgrb$set <- gsub("ISKS", "MGRB", Ex_tab_mgrb$set)

Ex_tab <- rbind.data.frame(Ex_tab_mgrb, Ex_tab_isks)
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
dim(Ex_tab[Ex_tab$gnomad_AF_NFE <= 0.0002,])
##AD
##C4/C5; gnomad_AF_NFE <=  0.0002
#Ex_tab_C4C5 <- Ex_tab[Ex_tab$auto_call %in% c("C4", "C5") & Ex_tab$gnomad_AF < 0.01,]
Ex_tab_C4C5 <- Ex_tab[Ex_tab$auto_call %in% c("C4", "C5") & Ex_tab$gnomad_AF_NFE <= 0.0002,]
##Add intra_cohort MAF filter cohort_MAF < 0.001 ##Added later
Ex_tab_C4C5 <- Ex_tab_C4C5[Ex_tab_C4C5$cohort_MAF <= 0.001 & Ex_tab_C4C5$swegen_AF <= 0.001,]
##C3; gnomad_AF == 0
#dim(Ex_tab[Ex_tab$auto_call %in% c("C3") & Ex_tab$gnomad_AF_NFE <= 0.0002,])
#Ex_tab_C3 <- Ex_tab[Ex_tab$auto_call %in% "C3" & Ex_tab$gnomad_AF <= 0.01 & Ex_tab$gnomad_AF_NFE == 0,]
Ex_tab_C3 <- Ex_tab[Ex_tab$auto_call %in% "C3" & Ex_tab$gnomad_AF_NFE  == 0 & Ex_tab$swegen_AF == 0,]
##cohort_MAF <= 0.001
Ex_tab_C3 <- Ex_tab_C3[Ex_tab_C3$cohort_MAF <= 0.001 & Ex_tab_C3$comb_score >= 5.6,]


##AR
#dim(Ex_tab[Ex_tab$gnomad_AF_NFE > 0.0002 & Ex_tab$gnomad_AF_NFE <= 0.01,])
Ex_tab_C4C5_AR <- Ex_tab[Ex_tab$auto_call %in% c("C4", "C5") & Ex_tab$gnomad_AF_NFE > 0.0002 & Ex_tab$gnomad_AF_NFE <= 0.01,]

Ex_tab_C4C5_AR_PID <- Ex_tab_C4C5_AR[Ex_tab_C4C5_AR$gene_symbol %in% gene_sub,]

##Add new filters
##GT =2; no X chromosome; gnomad < 0.05(already applied in HAIL); swegen <= 0.05, cohort_MAF <= 0.001 not needed
Ex_tab_C4C5_AR_PID_noX_GT2 <- Ex_tab_C4C5_AR_PID[Ex_tab_C4C5_AR_PID$GT == 2,]
Ex_tab_C4C5_AR_PID_noX_GT2 <- Ex_tab_C4C5_AR_PID_noX_GT2[!grepl("^X", Ex_tab_C4C5_AR_PID_noX_GT2$VARIANT),]
##latest changes; for GT == 1
Ex_tab_C4C5_AR_PID_GT1 <- Ex_tab_C4C5_AR_PID[Ex_tab_C4C5_AR_PID$GT == 1 & Ex_tab_C4C5_AR_PID$cohort_MAF <= 0.01 & 
                                               Ex_tab_C4C5_AR_PID$swegen_AF <= 0.01,]

###combined C4/5 from AD and AR; identify compound heterozygotes and retain the high confidence FS
Ex_tab_C4C5 <- Ex_tab_C4C5[Ex_tab_C4C5$comb_score >= 5.6,]
Ex_tab_C4C5_AR_PID_GT1_AD <- rbind.data.frame(Ex_tab_C4C5,Ex_tab_C4C5_AR_PID_GT1)

##Add same_gene_same_sample filter for complex indels with other variants within <= 5bp of it
#all_genes <- unique(Ex_tab$gene_symbol)
all_genes <- unique(Ex_tab_C4C5_AR_PID_GT1_AD$gene_symbol)
sam_sel_AR <- list()
for(i in 1:length(all_genes)){
  # gene_df <- Ex_tab[Ex_tab$gene_symbol %in% all_genes[i],]
  gene_df <- Ex_tab_C4C5_AR_PID_GT1_AD[Ex_tab_C4C5_AR_PID_GT1_AD$gene_symbol %in% all_genes[i],]
  if(dim(gene_df)[1] > 1){
    sam_freq <- as.data.frame(table(gene_df$SAMPLE))
    sam_sel <- as.character(sam_freq[sam_freq$Freq > 1,1])
    if(length(sam_sel) > 0){
      sam_sel_AR[[i]] <- gene_df[gene_df$SAMPLE %in% sam_sel,]
    }
    else{sam_sel_AR[[i]] <- NULL}
  }
  else { next }
  
}

All_cpd_het <- do.call("rbind.data.frame", sam_sel_AR)
All_cpd_het <- unique(All_cpd_het)

var_coods_all <- lapply(strsplit(All_cpd_het$VARIANT, ":"),
                        function(x)cbind.data.frame("chr" = x[1], "pos" = x[2]))

var_coods_df_all <-  do.call("rbind.data.frame", var_coods_all) 
var_coods_df_all$pos <- as.numeric(as.character(var_coods_df_all$pos))
All_cpd_het_pos <- cbind.data.frame(All_cpd_het, var_coods_df_all)

all_genes <- unique(All_cpd_het_pos$gene_symbol)
##########extract hits for removal; all variants within 5 bps of the fs variant
sam_sel_ns_rm <- list()
sam_sel_ns_rm_ls <- list()
for(i in 1:length(all_genes)){
  gene_df <- All_cpd_het_pos[All_cpd_het_pos$gene_symbol %in% all_genes[i],]
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
All_cpd_het_pos <- All_cpd_het_pos[All_cpd_het_pos$filt_tab %nin% sam_rm_het_df$filt_tab,]
##remove variants that are no longer compound het after filtering
library(dplyr)
s1 <- All_cpd_het_pos %>% group_by(SAMPLE, gene_symbol) %>% summarise(n = n()) 
s2 <- as.data.frame(s1)
s3 <- s2[s2$n > 1,]
All_cpd_het_pos_fin <- All_cpd_het_pos[All_cpd_het_pos$SAMPLE %in% s3$SAMPLE,]
All_cpd_het_pos_fin <- All_cpd_het_pos_fin[,-c(131:132)]


###AD: remove spurious compound hets
Ex_tab_C4C5 <- Ex_tab_C4C5[Ex_tab_C4C5$filt_tab %nin% sam_rm_het_df$filt_tab, ]

##remove compound heterozygotes including C3 variants too; this gets rid of the PTEN C3 variant
#Ex_tab_C3C4C5_AD <- rbind.data.frame(Ex_tab_C4C5, Ex_tab_C3)

##The compound hets from C3s removed with a cutoff of 3 missense variants/gene/sample in the DT_exome_clingene...R script

##Combine C3,C4 and C5 from AD and C4, C5 from AR
##latest file
Ex_tab_C3C4C5_AD_filt <- unique(rbind.data.frame(Ex_tab_C4C5, Ex_tab_C3))
rm_X_GT2_AD <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$GT == 2 & grepl("^X", Ex_tab_C3C4C5_AD_filt$VARIANT), ]
Ex_tab_C3C4C5_AD_filt <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$filt_tab %nin% rm_X_GT2_AD$filt_tab,]
Ex_tab_C3C4C5_AD_filt <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$comb_score >= 5.6,]

#AR should only contain compound het and homozygous variants, remove all AR GT == 1 that are not compound hets (DT Aug.26 meeting)
#Ex_tab_comb_C3C4C5 <- unique(rbind.data.frame(Ex_tab_C3C4C5_AD_filt, All_cpd_het_pos_fin, Ex_tab_C4C5_AR_PID_GT1, Ex_tab_C4C5_AR_PID_noX_GT2))
Ex_tab_comb_C3C4C5 <- unique(rbind.data.frame(Ex_tab_C3C4C5_AD_filt, All_cpd_het_pos_fin, Ex_tab_C4C5_AR_PID_noX_GT2))
Ex_tab_comb_C3C4C5$is_AR <- ifelse(Ex_tab_comb_C3C4C5$gnomad_AF_NFE > 0.0002, 1, 0)
##remove chromosome X variants from GT = 2
rm_X_GT2 <- Ex_tab_comb_C3C4C5[Ex_tab_comb_C3C4C5$GT == 2 & grepl("^X", Ex_tab_comb_C3C4C5$VARIANT), ]
Ex_tab_comb_C3C4C5 <- Ex_tab_comb_C3C4C5[Ex_tab_comb_C3C4C5$filt_tab %nin% rm_X_GT2$filt_tab, ]
#Ex_tab_comb_C3C4C5$is_AR <- ifelse(Ex_tab_comb_C3C4C5$gnomad_AF_NFE > 0.0002 | Ex_tab_comb_C3C4C5$GT == 2, 1, 0)
write.table(Ex_tab_comb_C3C4C5, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/comb_clin_isks_mgrb_2020_C3C4C5_NFE0002_AD_AR_all_fields_rnd3.tsv", 
            row.names = F, quote = F, sep = "\t")
Ex_tab_comb_C3C4C5_PID <- Ex_tab_comb_C3C4C5[Ex_tab_comb_C3C4C5$gene_symbol %in% gene_sub, ]
write.table(Ex_tab_comb_C3C4C5_PID, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/comb_clin_isks_mgrb_2020_C3C4C5_NFE0002_AD_AR_all_fields_rnd3_PID_genes.tsv", 
            row.names = F, quote = F, sep = "\t")
All_cpd_het_pos_fin$is_AR <- ifelse(All_cpd_het_pos_fin$gnomad_AF_NFE > 0.0002, 1, 0)
write.table(All_cpd_het_pos_fin, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/cmp_het_clin_isks_mgrb_2020_C4C5_NFE0002_AD_AR_all_fields_rnd3.tsv", 
            row.names = F, quote = F, sep = "\t")

# ##For SKAT(use the same set of variants used above, compound het filtered, gnomad_AF_NFE filtered)
# Ex_tab_C3C4C5_AD_filt <- unique(rbind.data.frame(Ex_tab_C4C5, Ex_tab_C3))
# rm_X_GT2_AD <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$GT == 2 & grepl("^X", Ex_tab_C3C4C5_AD_filt$VARIANT), ]
# Ex_tab_C3C4C5_AD_filt <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$filt_tab %nin% rm_X_GT2_AD$filt_tab,]
# Ex_tab_C3C4C5_AD_filt <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$comb_score >= 5.6,]
# write.table(Ex_tab_C3C4C5_AD_filt, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_all_fields_rnd3.tsv", 
#             row.names = F, quote = F, sep = "\t")
# Ex_tab_C3C4C5_AD_filt_PID <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$gene_symbol %in% gene_sub,]
# write.table(Ex_tab_C3C4C5_AD_filt_PID, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/all_isksmgrb_skatinp_PIDgenes_combset2020_clin_C3C4C5_NFE0002_AD_all_fields_rnd3.tsv", 
#             row.names = F, quote = F, sep = "\t")
