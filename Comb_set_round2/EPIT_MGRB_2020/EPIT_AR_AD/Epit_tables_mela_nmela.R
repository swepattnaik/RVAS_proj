##Separate Melanoma table; Melanoma may harbour Shelterin complex related variants and hence need to be 
## dealt separately; improve relative gain of Shelterin complex in ISKS
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)

fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/all_epitmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_all_fields_rnd3.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)

##Additive model
Ex_samp_id <- unique(fil_tab$SAMPLE)
#remove duplicates
dup_samp <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/fin_samp_dup_drop.txt", sep = "", header = T,
                       stringsAsFactors = F)

Ex_samp_id <- Ex_samp_id[Ex_samp_id %nin% dup_samp$x]

#Ex_var_id <- unique(fil_tab$VARIANT)

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
##note p_Data_PC_comb does not change
#p_Data <- read.table("~/RVAS/ISKS_MGRB_gender_age_3665.tsv", header = T, sep = "\t")
p_Data <- read.delim("~/RVAS/Epi_set_2020/pop_PCA/MGRB_EPIT_1000G_combset_pca.scores_clustered.tsv", header = T, sep = "\t", stringsAsFactors = F)
p_Data <- p_Data[p_Data$sample %in% QC2_epit_mgrb_pass$new_sampleid,]
#p_Data$sample <- gsub("isks", "", p_Data$sample)
#p_Data$sample <- gsub("risc", "", p_Data$sample)
p_Data_noCH <- p_Data[as.character(p_Data$sample) %in% Ex_samp_id,]

##Add gender information

p_Data_noCH$gender <- QC2_epit_mgrb_pass[match(p_Data_noCH$sample, QC2_epit_mgrb_pass$new_sampleid), 3]

##Add epithelial cancer phenotype
epit_pheno <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set/EPIT_MGRB_2020/qc0_final_all_annot_tel_length_20200623.tsv",
                         sep = "\t", header = T, stringsAsFactors = F)

#p_Data_noCH$study <- ifelse(p_Data_noCH$sample %in% epit_pheno$Sample_Id, epit_pheno$Study, "MGRB")
p_Data_noCH$study <- epit_pheno[match(p_Data_noCH$sample, epit_pheno$Sample_Id),19]
p_Data_noCH$study <- ifelse(is.na(p_Data_noCH$study), "MGRB", p_Data_noCH$study)

##remove 78 Familial breast cases
p_Data_noCH <- p_Data_noCH[p_Data_noCH$study %nin% "Familial breast",]
##Exlude 7 oesophageal cancer cases (withdrawn consent)
exclude_oeso <- epit_pheno[epit_pheno$EXCLUDE %in% "yes",1]
p_Data_noCH <- p_Data_noCH[p_Data_noCH$sample %nin% exclude_oeso,]
Ex_samp_id <- Ex_samp_id[match(p_Data_noCH$sample, Ex_samp_id)]
##filter out QC fail, 78 Familial breast cases and 7 oesophageal cases
fil_tab <- fil_tab[fil_tab$SAMPLE %in% Ex_samp_id,]

write.table(fil_tab,"~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/Tabs/EPIT_skat_inp_tab.tsv",
            sep = "\t", row.names = F, quote = F)

##melanoma table(does not include mucosal melanoma)
#p_Data_mela <- p_Data_noCH[p_Data_noCH$study %in% c("Melanoma","MGRB"),]
##melanoma table(include mucosal melanoma)Nov.2.2020
p_Data_mela <- p_Data_noCH[p_Data_noCH$study %in% c("Melanoma","Mucosal melanoma", "MGRB"),]
mela_samp_id <- Ex_samp_id[match(p_Data_mela$sample, Ex_samp_id)]
##filter out QC fail and 78 Familial breast cases
fil_tab_mela <- fil_tab[fil_tab$SAMPLE %in% mela_samp_id,]
write.table(fil_tab_mela,"~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/Tabs/EPIT_melanoma_skat_inp_tab.tsv",
            sep = "\t", row.names = F, quote = F)

##non-melanoma table(excludes melanoma and mucosal melanoma) Nov.2.2020
p_Data_nmela <- p_Data_noCH[p_Data_noCH$study %nin% c("Melanoma", "Mucosal melanoma"),]
nmela_samp_id <- Ex_samp_id[match(p_Data_nmela$sample, Ex_samp_id)]
##filter out QC fail and 78 Familial breast cases
fil_tab_nmela <- fil_tab[fil_tab$SAMPLE %in% nmela_samp_id,]
write.table(fil_tab_nmela,"~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/Tabs/EPIT_non_melanoma_skat_inp_tab.tsv",
            sep = "\t", row.names = F, quote = F)
epit_pheno_mod_age <- epit_pheno[epit_pheno$Sample_Id %in% p_Data_nmela$sample, ]
epit_pheno_mod_age$gender <- p_Data_nmela[match(epit_pheno_mod_age$Sample_Id, p_Data_nmela$sample),53]
epit_pheno_mod_age_gender <- epit_pheno_mod_age[,c(1,19,21,24,27)]

write.table(epit_pheno_mod_age_gender,"~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/Tabs/EPIT_age_gender_telo_non_melanoma.tsv",
            sep = "\t", row.names = F, quote = F)
