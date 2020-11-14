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
Ex_tab_isks <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/all_isks_combset2020_variants_AR_AD_all_fields_clingene_rnd3_freeze.tsv", sep = "\t", header = T, stringsAsFactors = F)
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
Ex_tab <- Ex_tab[!is.na(Ex_tab$SAMPLE),]

##remove duplicates
print("removing duplicates and QC2 fail samples")
fil_tab_noCH <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
Ex_tab <- Ex_tab[Ex_tab$SAMPLE %in% fil_tab_noCH$SAMPLE, ]
Ex_tab <- Ex_tab[!is.na(Ex_tab$SAMPLE),]
rm(fil_tab_noCH)
##AD
##C4/C5; gnomad_AF_NFE <=  0.0002
#Ex_tab_C4C5 <- Ex_tab[Ex_tab$auto_call %in% c("C4", "C5") & Ex_tab$gnomad_AF < 0.01,]
dim(Ex_tab[Ex_tab$gnomad_AF_NFE <= 0.0002,])[1]
table(Ex_tab[Ex_tab$gnomad_AF_NFE <= 0.0002,]$set)
Ex_tab_C4C5 <- Ex_tab[Ex_tab$auto_call %in% c("C4", "C5") & Ex_tab$gnomad_AF_NFE <= 0.0002,]
dim(Ex_tab_C4C5)[1]
table(Ex_tab_C4C5$set)
##Add intra_cohort MAF filter cohort_MAF < 0.001 ##Added later
Ex_tab_C4C5 <- Ex_tab_C4C5[Ex_tab_C4C5$cohort_MAF <= 0.001 & Ex_tab_C4C5$swegen_AF <= 0.001,]
dim(Ex_tab_C4C5)[1]
table(Ex_tab_C4C5$set)
##C3; gnomad_AF == 0
#Ex_tab_C3 <- Ex_tab[Ex_tab$auto_call %in% "C3" & Ex_tab$gnomad_AF <= 0.01 & Ex_tab$gnomad_AF_NFE == 0,]
dim(Ex_tab[Ex_tab$auto_call %in% c("C3") & Ex_tab$gnomad_AF_NFE <= 0.0002,])[1]
table(Ex_tab[Ex_tab$auto_call %in% c("C3") & Ex_tab$gnomad_AF_NFE <= 0.0002,]$set)
Ex_tab_C3 <- Ex_tab[Ex_tab$auto_call %in% "C3" & Ex_tab$gnomad_AF_NFE  == 0 & Ex_tab$swegen_AF == 0,]
##cohort_MAF <= 0.001
Ex_tab_C3 <- Ex_tab_C3[Ex_tab_C3$cohort_MAF <= 0.001 & Ex_tab_C3$comb_score >= 5.6,]
dim(Ex_tab_C3)[1]
table(Ex_tab_C3$set)
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

dim(Ex_tab_C4C5)[1]
table(Ex_tab_C4C5$set)

##For SKAT(use the same set of variants used above, compound het filtered, gnomad_AF_NFE filtered)
Ex_tab_C3C4C5_AD_filt <- unique(rbind.data.frame(Ex_tab_C4C5, Ex_tab_C3))
dim(Ex_tab_C3C4C5_AD_filt)[1]
table(Ex_tab_C3C4C5_AD_filt$set)
rm_X_GT2_AD <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$GT == 2 & grepl("^X", Ex_tab_C3C4C5_AD_filt$VARIANT), ]
Ex_tab_C3C4C5_AD_filt <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$filt_tab %nin% rm_X_GT2_AD$filt_tab,]
Ex_tab_C3C4C5_AD_filt <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$comb_score >= 5.6,]
dim(Ex_tab_C3C4C5_AD_filt)[1]
table(Ex_tab_C3C4C5_AD_filt$set) ##Matched SKAT input ~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv

#write.table(Ex_tab_C3C4C5_AD_filt, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_all_fields_rnd3_freeze.tsv", 
#            row.names = F, quote = F, sep = "\t")
#Ex_tab_C3C4C5_AD_filt_PID <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$gene_symbol %in% gene_sub,]
#write.table(Ex_tab_C3C4C5_AD_filt_PID, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/all_isksmgrb_skatinp_PIDgenes_combset2020_clin_C3C4C5_NFE0002_AD_all_fields_rnd3.tsv", 
#            row.names = F, quote = F, sep = "\t")

###################For PIDfile
#PIDgenes
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

##AD; GT = 2; For DT, July 6
Ex_tab_C4C5_noX_GT2_AD <- Ex_tab_C4C5[Ex_tab_C4C5$GT == 2 & !grepl("^X", Ex_tab_C4C5$VARIANT),]
Ex_tab_C4C5_noX_GT2_AD_PID <- Ex_tab_C4C5_noX_GT2_AD[Ex_tab_C4C5_noX_GT2_AD$gene_symbol %in% gene_sub,]
Ex_tab_C3_noX_GT2_AD <- Ex_tab_C3[Ex_tab_C3$GT == 2 & !grepl("^X", Ex_tab_C3$VARIANT),]
Ex_tab_C3_noX_GT2_AD_PID <- Ex_tab_C3_noX_GT2_AD[Ex_tab_C3_noX_GT2_AD$gene_symbol %in% gene_sub,]
# write.table(rbind.data.frame(Ex_tab_C4C5_noX_GT2_AD_PID, Ex_tab_C3_noX_GT2_AD_PID), "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/AD_clin_C3C4C5_GT2_noX_PID_genes.tsv",
#             row.names = F, quote = F, sep = "\t")
# write.table(rbind.data.frame(Ex_tab_C4C5_noX_GT2_AD, Ex_tab_C3_noX_GT2_AD), "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/AD_clin_C3C4C5_GT2_noX_all_genes.tsv",
#             row.names = F, quote = F, sep = "\t")

##AR
print("Processing AR")
Ex_tab_C4C5_AR <- Ex_tab[Ex_tab$auto_call %in% c("C4", "C5") & Ex_tab$gnomad_AF_NFE > 0.0002 & Ex_tab$gnomad_AF_NFE <= 0.01,]
dim(Ex_tab[Ex_tab$gnomad_AF_NFE > 0.0002 & Ex_tab$gnomad_AF_NFE <= 0.01,])[1]
table(Ex_tab[Ex_tab$gnomad_AF_NFE > 0.0002 & Ex_tab$gnomad_AF_NFE <= 0.01,]$set)

dim(Ex_tab_C4C5_AR)[1]
table(Ex_tab_C4C5_AR$set)

Ex_tab_C4C5_AR_PID <- Ex_tab_C4C5_AR[Ex_tab_C4C5_AR$gene_symbol %in% gene_sub,]
dim(Ex_tab_C4C5_AR_PID)[1]
table(Ex_tab_C4C5_AR_PID$set)

##Add new filters for AR for PID list
#AR should only contain compound het and homozygous variants, remove all AR GT == 1 that are not compound hets (DT Aug.26 meeting)
##GT =2; no X chromosome; gnomad < 0.05(already applied in HAIL); swegen <= 0.05, cohort_MAF <= 0.001 not needed
Ex_tab_C4C5_AR_PID_noX_GT2 <- Ex_tab_C4C5_AR_PID[Ex_tab_C4C5_AR_PID$GT == 2,]
Ex_tab_C4C5_AR_PID_noX_GT2 <- Ex_tab_C4C5_AR_PID_noX_GT2[!grepl("^X", Ex_tab_C4C5_AR_PID_noX_GT2$VARIANT),]
dim(Ex_tab_C4C5_AR_PID_noX_GT2)[1]
table(Ex_tab_C4C5_AR_PID_noX_GT2$set)
##latest changes (for compound het generation)
Ex_tab_C4C5_AR_PID_GT1 <- Ex_tab_C4C5_AR_PID[Ex_tab_C4C5_AR_PID$GT == 1 & Ex_tab_C4C5_AR_PID$cohort_MAF <= 0.01 & 
                                               Ex_tab_C4C5_AR_PID$swegen_AF <= 0.01,]
dim(Ex_tab_C4C5_AR_PID_GT1)[1]
table(Ex_tab_C4C5_AR_PID_GT1$set)
##combined C4/5 from AD and AR; identify compound heterozygotes and retain the high confidence FS
Ex_tab_C4C5_AR_PID_GT1_AD <- rbind.data.frame(Ex_tab_C4C5,Ex_tab_C4C5_AR_PID_GT1)
dim(Ex_tab_C4C5_AR_PID_GT1_AD)[1]
table(Ex_tab_C4C5_AR_PID_GT1_AD$set)####continue from here

##Add same_gene_same_sample filter for complex indels with other variants within <= 5bp of it
#all_genes <- unique(Ex_tab$gene_symbol)
print("Filtering GATK FS artefacts")
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

print("Extracting Compound hets")
var_coods_all <- lapply(strsplit(All_cpd_het$VARIANT, ":"),
                        function(x)cbind.data.frame("chr" = x[1], "pos" = x[2]))

var_coods_df_all <-  do.call("rbind.data.frame", var_coods_all) 
var_coods_df_all$pos <- as.numeric(as.character(var_coods_df_all$pos))
All_cpd_het_pos <- cbind.data.frame(All_cpd_het, var_coods_df_all)

all_genes <- unique(All_cpd_het_pos$gene_symbol)
#extract hits for removal; all variants within 5 bps of the fs variant
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
##remove variants that are not compound het after filtering
s1 <- All_cpd_het_pos %>% group_by(SAMPLE, gene_symbol) %>% summarise(n = n()) 
s2 <- as.data.frame(s1)
s3 <- s2[s2$n > 1,]
All_cpd_het_pos_fin <- All_cpd_het_pos[All_cpd_het_pos$SAMPLE %in% s3$SAMPLE,]
All_cpd_het_pos_fin <- All_cpd_het_pos_fin[,-c(131:132)]
dim(All_cpd_het_pos)[1]
table(All_cpd_het_pos$set)
#Ex_tab <- Ex_tab[Ex_tab$filt_tab %nin% sam_rm_het_df$filt_tab,]
#Ex_tab_C4C5_AD <- Ex_tab_C4C5[Ex_tab_C4C5$filt_tab %nin% rm_fs$filt_tab, ]
##AD
Ex_tab_C4C5 <- Ex_tab_C4C5[Ex_tab_C4C5$filt_tab %nin% sam_rm_het_df$filt_tab, ]

##remove compound heterozygotes including C3 variants too; this gets rid of the PTEN C3 variant
#Ex_tab_C3C4C5_AD <- rbind.data.frame(Ex_tab_C4C5, Ex_tab_C3)


##Combine C3,C4 and C5
##latest file
print("Saving PID file input")
Ex_tab_comb_C3C4C5 <- unique(rbind.data.frame(Ex_tab_C4C5, Ex_tab_C3, All_cpd_het_pos_fin, Ex_tab_C4C5_AR_PID_noX_GT2))
##remove chromosome X variants from GT = 2
rm_X_GT2 <- Ex_tab_comb_C3C4C5[Ex_tab_comb_C3C4C5$GT == 2 & grepl("^X", Ex_tab_comb_C3C4C5$VARIANT), ]
Ex_tab_comb_C3C4C5 <- Ex_tab_comb_C3C4C5[Ex_tab_comb_C3C4C5$filt_tab %nin% rm_X_GT2$filt_tab, ]
Ex_tab_comb_C3C4C5$is_AR <- ifelse(Ex_tab_comb_C3C4C5$gnomad_AF_NFE > 0.0002 | Ex_tab_comb_C3C4C5$GT == 2, 1, 0)
All_cpd_het_pos_fin$is_AR <- ifelse(All_cpd_het_pos_fin$gnomad_AF_NFE > 0.0002, 1, 0)
#write.table(All_cpd_het_pos_fin, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/cmp_het_clin_isks_mgrb_2020_C4C5_NFE0002_AD_AR_all_fields_rnd3_Aug31.tsv", 
#            row.names = F, quote = F, sep = "\t")

##for PID file(numbers for flowchart)
Ex_tab_comb_C3C4C5_PID <- Ex_tab_comb_C3C4C5[Ex_tab_comb_C3C4C5$gene_symbol %in% gene_sub,]
dim(Ex_tab_comb_C3C4C5_PID)[1]
table(Ex_tab_comb_C3C4C5_PID$set)
C45_GT2_var <- rbind.data.frame(All_cpd_het_pos_fin[,-131], Ex_tab_C4C5_AR_PID_noX_GT2)
dim(C45_GT2_var)[1]
table(C45_GT2_var$set)
##SKAT input
dim(fil_tab_noCH[fil_tab_noCH$gene_symbol %in% gene_sub,])
table(fil_tab_noCH[fil_tab_noCH$gene_symbol %in% gene_sub,]$set)


##########new approach

Ex_tab_01 <- Ex_tab[Ex_tab$gnomad_AF_NFE <= 0.01,]
##Add PID gene filter
Ex_tab_01_PID <- Ex_tab_01[Ex_tab_01$gene_symbol %in%  gene_sub,]
##C45 filter
Ex_tab_01_PID_C45 <- Ex_tab_01_PID[Ex_tab_01_PID$auto_call %in% c("C4", "C5"),]
##Get homozygous and not X
Ex_tab_01_PID_C45_noX_GT2 <- Ex_tab_01_PID_C45[Ex_tab_01_PID_C45$GT == 2 & 
                                                 !grepl("^X", Ex_tab_01_PID_C45$VARIANT),]
##Get GT == 1 with swegen <= 0.01 and cohort_MAF <= 0.01
Ex_tab_01_PID_C45_GT1 <- Ex_tab_01_PID_C45[Ex_tab_01_PID_C45$GT == 1 & 
                                             Ex_tab_01_PID_C45$swegen_AF <= 0.01 & 
                                             Ex_tab_01_PID_C45$cohort_MAF <= 0.01,]
##Extract compound hets
#function to remove GATK artefacts
all_genes <- unique(Ex_tab_01_PID_C45_GT1$gene_symbol)

##assign coordinates
var_coods_all <- lapply(strsplit(Ex_tab_01_PID_C45_GT1$VARIANT, ":"),
                        function(x)cbind.data.frame("chr" = x[1], "pos" = x[2]))

var_coods_df_all <-  do.call("rbind.data.frame", var_coods_all) 
var_coods_df_all$pos <- as.numeric(as.character(var_coods_df_all$pos))
Ex_tab_01_PID_C45_GT1_pos <- cbind.data.frame(Ex_tab_01_PID_C45_GT1, var_coods_df_all)

rm_artefacts <- function(gene_list, variant_file){
#extract hits for removal; all variants within 5 bps of the fs variant
sam_sel_ns_rm <- list()
sam_sel_ns_rm_ls <- list()
for(i in 1:length(gene_list)){
  gene_df <- variant_file[variant_file$gene_symbol %in% gene_list[i],]
  sam_freq <- as.data.frame(table(gene_df$SAMPLE))
  sam_sel <- as.character(sam_freq[sam_freq$Freq >1,1])
  if(length(sam_sel) > 0){
    print(gene_list[i])
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
return(sam_rm_het_df)
}

art_rm_df <- rm_artefacts(all_genes, Ex_tab_01_PID_C45_GT1_pos)
Ex_tab_01_PID_C45_GT1_pos <- Ex_tab_01_PID_C45_GT1_pos[Ex_tab_01_PID_C45_GT1_pos$filt_tab %nin% art_rm_df$filt_tab,]

##Compound het extraction
sam_sel_AR <- list()
for(i in 1:length(all_genes)){
  # gene_df <- Ex_tab[Ex_tab$gene_symbol %in% all_genes[i],]
  gene_df <- Ex_tab_01_PID_C45_GT1_pos[Ex_tab_01_PID_C45_GT1_pos$gene_symbol %in% all_genes[i],]
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



####From AD side; add only C3 variants with the right filters

Ex_tab_C3 <- Ex_tab[Ex_tab$auto_call %in% "C3" & Ex_tab$gnomad_AF_NFE  == 0 & Ex_tab$swegen_AF == 0,]
##cohort_MAF <= 0.001
Ex_tab_C3 <- Ex_tab_C3[Ex_tab_C3$cohort_MAF <= 0.001 & Ex_tab_C3$comb_score >= 5.6,]
##PID filter
Ex_tab_C3_PID <- Ex_tab_C3[Ex_tab_C3$gene_symbol %in% gene_sub,]

##C45 :AD
Ex_tab_C45 <- fil_tab_noCH[fil_tab_noCH$gene_symbol %in% gene_sub & fil_tab_noCH$auto_call %in% c("C4", "C5"),]

