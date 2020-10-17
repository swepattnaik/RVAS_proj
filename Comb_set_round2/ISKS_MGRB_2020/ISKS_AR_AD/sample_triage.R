##samples with no genomic data
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)

##PID file
t1 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/PID_combset2020_CGC_skatBin_repstress_potint_mito_chkpt_predNFE_clueGO_Apr292020_AD_Aug19.tsv", 
                 sep = "\t", header = T, stringsAsFactors = F)

na_geno <- t1[is.na(t1$nC5),]$pmn
##duplicated
dup_sam <- t1[!is.na(t1$duplicate),]$pmn
dup_sam_dup <- t1[!is.na(t1$duplicate),]$duplicate
##Add batch data to t1
t1$batch <- QC2_dat[match(t1$pmn, QC2_dat$new_sampleid),2]

t1[!is.na(t1$duplicate),c(1,3,85:87),]
t1[!is.na(t1$duplicate) & t1$duplicate %in% dup_ne,c(1,3,85:87),]
t1[!is.na(t1$duplicate) & t1$duplicate %in% dup_ne,c(1,3,87:88),] ##uses QC2_dat, dup_ne objects
##unknown QC2
unk_sam <- t1[t1$QC2 %in% "Unknown",]$pmn
##no PID variants

na_geno_pid <- na_geno[na_geno %nin% c(dup_sam, unk_sam)]

##joint called files not present in PID files

jc_isks <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/ISKS.samples.txt",
                      sep = "", stringsAsFactors = F, header = F)
jc_no_pid <- jc_isks$V1[jc_isks$V1 %nin% t1$pmn]
pid_no_jc <- t1$pmn[t1$pmn %nin% jc_isks$V1]
table(jc_isks$V1 %in% t1$pmn)

##joint called QC file
QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/joint_calls_2019jan_2019nov.final_qc_output.tsv", header = T, sep = "\t", stringsAsFactors = F)

QC2_dat[QC2_dat$new_sampleid %in% dup_sam,]$joint_call
QC2_dat[QC2_dat$new_sampleid %in% jc_no_pid,]$joint_call ##all joint called ISKS were from nov2019 batch
QC2_dat[QC2_dat$new_sampleid %in% pid_no_jc,]$joint_call ##all joint called ISKS were from nov2019 batch


#not joint called; the samples with unknown QC2 were not joint called; only 1643/1661 were called
##1661 not 1663 as 3083,CR57(not sure) were removed
##All except 3083 has variants; 3083 was not joint called as the gvcf was not found(Shyam's email confirms it)
##Note 3083 is not European (A8_Ethnic group1 is Caribean Black)
table(t1$pmn %in% jc_isks$V1) 
table(unk_sam %in% jc_isks$V1)


##PID files need to be renamed using Mandy's file

ren_sam <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/samp_name_rect.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
table(t1[!is.na(t1$duplicate),]$duplicate %in% ren_sam$Incorrect_UDF.External.ID)

##filtered ISKS set
fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/comb_clin_isks_mgrb_2020_C3C4C5_NFE0002_AD_AR_all_fields_rnd3.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)
fil_tab$is_case <- ifelse(grepl("^ISKS", fil_tab$set), 1, 0)
fil_tab_isks_set <- fil_tab[fil_tab$is_case == 1,]
fil_tab_isks_set$auto_call <- ifelse(fil_tab_isks_set$is_AR == 1 & fil_tab_isks_set$GT == 2, paste0(fil_tab_isks_set$auto_call, "_ar"), fil_tab_isks_set$auto_call)

isks_sam <- unique(fil_tab_isks_set$SAMPLE)
table(t1$pmn %in% isks_sam)
table(na_geno_pid %in% isks_sam)
table(unk_sam %in% isks_sam)


##All except 3083 has variants; 3083 was not joint called as the gvcf was not found(Shyam's email confirms it)
##Note 3083 is not European (A8_Ethnic group1 is Caribean Black)

##PID file rename
t1$rect_name <- ren_sam[match(t1$pmn,ren_sam$Incorrect_UDF.External.ID), 2]
t1$rect_pid <- t1[match(t1$rect_name, t1$pmn),1]

all_ren_sam <- unique(c(ren_sam$Incorrect_UDF.External.ID, ren_sam$rect_ID))

rect_t1 <- t1[!is.na(t1$rect_name),c(1,3,87:89)]
t1[t1$pmn %in% all_ren_sam & t1$Sheltrin_main == 1,c(68:79, 87:89)]
t1[t1$pmn %in% all_ren_sam & t1$Sheltrin_extn == 1,c(68:79, 87:89)]

##dup samples not explained
dup_ne <- dup_sam_dup[dup_sam_dup %nin% ren_sam[ren_sam$Incorrect_UDF.External.ID %in% dup_sam_dup,]$Incorrect_UDF.External.ID]
dup_ne_dup <- dup_sam[dup_sam_dup %nin% ren_sam[ren_sam$Incorrect_UDF.External.ID %in% dup_sam_dup,]$Incorrect_UDF.External.ID]
t1[!is.na(t1$duplicate),c(1,3,87:89),]


rect_jan2019 <- ren_sam[ren_sam$Incorrect_UDF.External.ID %in% dup_sam_dup,]
ren_sam[ren_sam$Incorrect_UDF.External.ID %in% dup_sam,]

table(QC2_dat[QC2_dat$new_sampleid %in% ren_sam$Incorrect_UDF.External.ID,2])





##check denia's file

comb_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/PID_Master_file_290420.xlsx", sheet = 1, col_types = c("list"))
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

comb_pheno_1 <- comb_pheno
comb_pheno_2 <- comb_pheno


##PID variant set


comb_set_filt1 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CGC_topSKATBin300_repstress_potint_mito_chkpt_VARIANTS_filt_combset2020_clueGOplus_Aug19.tsv",
                             sep = "\t", header = T, stringsAsFactors = F)
pid_filt_sam <- unique(comb_set_filt1$SAMPLE)
table(pid_filt_sam %in% isks_sam)
table(pid_filt_sam %in% jc_no_pid)
table(t1$pmn %nin% pid_filt_sam)
table(t1$pmn %nin% pid_filt_sam)


######collate all data
##joint called QC file
QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/joint_calls_2019jan_2019nov.final_qc_output.tsv", header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_fail <- QC2_dat[QC2_dat$passes_qc2 %in% "FALSE",]$new_sampleid
##all joint called data (n = 1661)
all_jc_sam <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/ISKS.samples.txt",
                         sep = "", stringsAsFactors = F, header = F) ## n = 1661
all_jc_sam <- all_jc_sam$V1

##duplicates
deg2_isks <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/dup_pairs_and_deg2_isks.tsv", sep = "\t",
                        header = T, stringsAsFactors = F)
deg2_isks <- deg2_isks[deg2_isks$degree == 0,]
deg2_isks <- deg2_isks[,7:8]

##ISKS samples
fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/comb_clin_isks_mgrb_2020_C3C4C5_NFE0002_AD_AR_all_fields_rnd3.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)
fil_tab$is_case <- ifelse(grepl("^ISKS", fil_tab$set), 1, 0)
fil_tab_isks_set <- fil_tab[fil_tab$is_case == 1,]
all_isks_sam <- unique(fil_tab_isks_set$SAMPLE) ## n = 1646 (minus 15 duplicates)

##samples in enriched set (CGC+SKATbin+O-PINT) 
comb_set_filt1 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CGC_topSKATBin300_repstress_potint_mito_chkpt_VARIANTS_filt_combset2020_clueGOplus_Aug19.tsv",
                             sep = "\t", header = T, stringsAsFactors = F)
#d1 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CGC_topSKATBin300_repstress_potint_mito_chkpt_VARIANTS_filt_combset2020_clueGOplus_Aug19.tsv",
#              sep = "\t", header = T, stringsAsFactors = F)
pid_filt_sam <- unique(comb_set_filt1$SAMPLE) ## n = 1575 

##samples that didn't make it to PID set
isks_no_pid <- all_isks_sam[all_isks_sam %nin% pid_filt_sam] ## n = 71

#####

##All PID file samples
pid_fil <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/PID_combset2020_CGC_skatBin_repstress_potint_mito_chkpt_predNFE_clueGO_Apr292020_AD_Aug19.tsv", 
           sep = "\t", header = T, stringsAsFactors = F)
pid_all_sam <- pid_fil$pmn ## n = 1663

##samples with genomic variants 
geno_sam <- pid_fil[!is.na(pid_fil$nC5),]$pmn ## n = 1559

##Explore why 104 samples didn't have genomic coordinates
##no genomic variants samples
pid_na_geno <- pid_fil[is.na(pid_fil$nC5),]$pmn ## n = 104
pid_na_geno_no_var <- pid_na_geno[pid_na_geno %in% isks_no_pid] ## n = 69
##Add batch data to pid_fil
pid_fil$batch <- QC2_dat[match(pid_fil$pmn, QC2_dat$new_sampleid),2]
pid_fil$batch_dup <- QC2_dat[match(pid_fil$duplicate, QC2_dat$new_sampleid),2]
##duplicated samples; 15 are duplicates (these were removed, hence no genomic variants)
dup_sam <- pid_fil[!is.na(pid_fil$duplicate),]$pmn ## n = 15,
dup_sam_dup <- pid_fil[!is.na(pid_fil$duplicate),]$duplicate ## n = 15
pid_fil_dup <- pid_fil[!is.na(pid_fil$duplicate),c(79,87:89)]
pid_fil_dup$ren_pmn <- ren_sam[match(pid_fil_dup$sample_id, ren_sam$Incorrect_UDF.External.ID),2]
pid_fil_dup$ren_dup <- ren_sam[match(pid_fil_dup$duplicate, ren_sam$Incorrect_UDF.External.ID),2]
table(dup_sam_dup %in% isks_no_pid)
dup_sam_dup_no_geno <- dup_sam_dup[dup_sam_dup %in% isks_no_pid] ## n = 3; duplicated sample with no genomic variants 
##20/104 samples still needs to be explained
#check overlap between joint called file and pid file
table(pid_all_sam %in% all_jc_sam)
pid_com_jc <- pid_all_sam[pid_all_sam %in% all_jc_sam] ## n = 1643
pid_com_jc_rmdup <- pid_com_jc[pid_com_jc %nin% dup_sam] ## n = 1628
#pid files not joint called; n = 20
pid_no_jc <- pid_all_sam[pid_all_sam %nin% all_jc_sam] ##n = 20 (19 CRs + 3083 case) 
##note 3083 gvcf was missing
##done
##Summary of PID NAs:
##66 No_variants + 12 duplicates + 3 duplicates_no_variants + 20_not_joint_called

##Add samples that are exclusive to ISKS joint call to PID file
excl_isks_sam <- all_isks_sam[all_isks_sam %nin% pid_all_sam] #n = 18 exclusively ISKS.
table(excl_isks_sam %in% dup_sam) ##none are duplicated
table(excl_isks_sam %in% isks_no_pid) ## 16 have variants for PID file while 2 have none.

excl_isks_sam_16 <- excl_isks_sam[excl_isks_sam %nin% isks_no_pid]

##duplicate analysis

dup_sam <- pid_fil[!is.na(pid_fil$duplicate),]$pmn ## n = 15,
dup_sam_dup <- pid_fil[!is.na(pid_fil$duplicate),]$duplicate ## n = 15
pid_fil_dup <- pid_fil[!is.na(pid_fil$duplicate),c(79,87:89)]
pid_fil_dup$ren_pmn <- ren_sam[match(pid_fil_dup$sample_id, ren_sam$Incorrect_UDF.External.ID),2]
pid_fil_dup$ren_dup <- ren_sam[match(pid_fil_dup$duplicate, ren_sam$Incorrect_UDF.External.ID),2]

######
##Joint call input files
##M012 files
rect_ren_sam <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/ISKS_samples_JC4_freeze.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)

rect_ren_sam_M012 <- rect_ren_sam[rect_ren_sam$Manifest %in% "R_160922_MANBAL_ISKSDNA_M012",]$JCInputID
rect_ren_sam[rect_ren_sam$Same. == "FALSE" & rect_ren_sam$check == "TRUE" & rect_ren_sam$Manifest %in% "R_160922_MANBAL_ISKSDNA_M012",5:10]

table(ren_sam$Incorrect_UDF.External.ID %in% rect_ren_sam_M012)
table(ren_sam$rect_ID %in% rect_ren_sam_M012)
ren_sam$rect_ID[ren_sam$rect_ID %nin% rect_ren_sam_M012]
table(ren_sam$rect_ID[ren_sam$rect_ID %nin% rect_ren_sam_M012] %in% QC2_dat_fail)
QC2_dat_fail[QC2_dat_fail %in% ren_sam$rect_ID[ren_sam$rect_ID %nin% rect_ren_sam_M012]]
table(ren_sam$rect_ID %in% QC2_dat_fail)



#rect_ren_sam[duplicated(rect_ren_sam$SampleName)|duplicated(rect_ren_sam$SampleName, fromLast=TRUE),5:10]
##remove duplicated cases from joint called from M012 due to mislabelling
rect_ren_sam_rmdup1 <- rect_ren_sam[!duplicated(rect_ren_sam$SampleName, fromLast=TRUE),]
##remove duplicated cases due to sequencing as a part of two studies;
##retain the case that is from a more recent manifest if they are from the same study
##else retain the cases that pass both QC1 and QC2
library(readxl)
isks_risc_sam <- read_xlsx("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/duplicate_analysis.xlsx",
                           sheet = 3)
isks_risc_sam <- as.data.frame(isks_risc_sam)
isks_risc_sam <- isks_risc_sam[,-5]
#table(rect_ren_sam_rmdup1$JCInputID %in% isks_risc_sam$`RisC Blood ID`)
##fails QC1; remove "CR184", "CR693"
isks_risc_sam_rm <- c("CR184", "CR693")
rect_ren_sam_rmdup2 <- rect_ren_sam_rmdup1[rect_ren_sam_rmdup1$JCInputID %nin% c("CR184", "CR693"),]

##duplicates not explained by M012: c("2085","1371", "1762", "2349") 
##retain: c("2085","1371")
##remove: c("1762", "2349") ##have been sequenced twice in ISKS cohort, see ISKS Dup sheet (Google sheet)
rect_ren_sam_rmdup3 <- rect_ren_sam_rmdup2[rect_ren_sam_rmdup2$JCInputID %nin% c("1762", "2349"),]
##Add LIONS
Lions_sam <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/ISKS_samples_LIONS_freeze.tsv",
                        sep = "\t", header = T, stringsAsFactors = F)
Lions_sam_rmdup <- Lions_sam[!duplicated(Lions_sam$SampleName),]
##remove LIONs cases from duplicated samples
##remove c("LKCGP-P001604-266610-02-04-07-G1", "LKCGP-P000505-266108-02-04-07-G1")
Lions_sam_rmdup1 <-Lions_sam_rmdup[Lions_sam_rmdup$ExtID %nin% c("LKCGP-P001604-266610-02-04-07-G1", "LKCGP-P000505-266108-02-04-07-G1"),] 

##combine ISKS+RISC+LIONS
rect_ren_sam_all <- rbind.data.frame(rect_ren_sam_rmdup3,Lions_sam_rmdup1)
##remove CR57 from rectified ID column; total 1645 cases
rect_ren_sam_all <- rect_ren_sam_all[rect_ren_sam_all$JCInputRecID %nin% "CR57",]
##share this list with Denia for generating PID
write.table(rect_ren_sam_all, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/ISKS_RISC_LIONS_final_freeze.tsv",
            sep = "\t", row.names = F, quote = F)