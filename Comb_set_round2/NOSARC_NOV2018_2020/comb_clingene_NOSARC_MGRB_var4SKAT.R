###combine relaxed gnomad_AF ISKS and MGRB for SKAT input
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

Ex_tab_mgrb <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/MGRB/all_mgrb_combset2020_variants_AR_AD_all_fields_clingene_rnd3.tsv", sep = "\t", header = T, stringsAsFactors = F)
Ex_tab_sarc <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/TCGA_SARC_AR_AD/all_tcga_sarc_combset2020_variants_AR_AD_all_fields_clingene_rnd3.tsv", sep = "\t", header = T, stringsAsFactors = F)
Ex_tab_mgrb$set <- gsub("ISKS", "MGRB", Ex_tab_mgrb$set)

Ex_tab <- rbind.data.frame(Ex_tab_mgrb, Ex_tab_sarc)
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
##Add intra_cohort MAF filter cohort_MAF < 0.001 ##Added later
Ex_tab_C4C5 <- Ex_tab_C4C5[Ex_tab_C4C5$cohort_MAF <= 0.001 & Ex_tab_C4C5$swegen_AF <= 0.001,]
##C3; gnomad_AF == 0
#Ex_tab_C3 <- Ex_tab[Ex_tab$auto_call %in% "C3" & Ex_tab$gnomad_AF <= 0.01 & Ex_tab$gnomad_AF_NFE == 0,]
Ex_tab_C3 <- Ex_tab[Ex_tab$auto_call %in% "C3" & Ex_tab$gnomad_AF_NFE  == 0 & Ex_tab$swegen_AF == 0,]
##cohort_MAF <= 0.001
Ex_tab_C3 <- Ex_tab_C3[Ex_tab_C3$cohort_MAF <= 0.001 & Ex_tab_C3$comb_score >= 5.6,]

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
write.table(Ex_tab_C3C4C5_AD_filt, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/TCGA_SARC_AR_AD/all_tcgasarc_mgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_all_fields_rnd3_Aug31.tsv", 
            row.names = F, quote = F, sep = "\t")
