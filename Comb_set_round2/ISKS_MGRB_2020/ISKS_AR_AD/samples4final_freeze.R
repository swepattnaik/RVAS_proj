######
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)
##adopted from sample_triage.R 
##Joint call input files
##M012 files
rect_ren_sam <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/ISKS_samples_JC4_freeze.tsv",
                           sep = "\t", stringsAsFactors = F, header = T)

# rect_ren_sam_M012 <- rect_ren_sam[rect_ren_sam$Manifest %in% "R_160922_MANBAL_ISKSDNA_M012",]$JCInputID
# rect_ren_sam[rect_ren_sam$Same. == "FALSE" & rect_ren_sam$check == "TRUE" & rect_ren_sam$Manifest %in% "R_160922_MANBAL_ISKSDNA_M012",5:10]
# 
# table(ren_sam$Incorrect_UDF.External.ID %in% rect_ren_sam_M012)
# table(ren_sam$rect_ID %in% rect_ren_sam_M012)
# ren_sam$rect_ID[ren_sam$rect_ID %nin% rect_ren_sam_M012]
# table(ren_sam$rect_ID[ren_sam$rect_ID %nin% rect_ren_sam_M012] %in% QC2_dat_fail)
# QC2_dat_fail[QC2_dat_fail %in% ren_sam$rect_ID[ren_sam$rect_ID %nin% rect_ren_sam_M012]]
# table(ren_sam$rect_ID %in% QC2_dat_fail)



#rect_ren_sam[duplicated(rect_ren_sam$SampleName)|duplicated(rect_ren_sam$SampleName, fromLast=TRUE),5:10]
##remove duplicated cases from joint called from M012 due to mislabelling
rect_ren_sam_rmdup1 <- rect_ren_sam[!duplicated(rect_ren_sam$SampleName, fromLast=TRUE),] 
#rect_ren_sam_rmdup1 <- rect_ren_sam[!duplicated(rect_ren_sam$SampleName),]##this removes samples
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
#isks_risc_sam_rm <- c("CR184", "CR693")
rect_ren_sam_rmdup2 <- rect_ren_sam_rmdup1[rect_ren_sam_rmdup1$JCInputID %nin% c("CR184", "CR693"),]

##duplicates not explained by M012: c("2085","1371", "1762", "2349") 
##retain: c("2349","1371") ##2349 in present in M014 and should be retained and removed from M012(ID: 2085)
##remove: c("1762", "2085") ##have been sequenced twice in ISKS cohort, see ISKS Dup sheet (Google sheet)
rect_ren_sam_rmdup3 <- rect_ren_sam_rmdup2[rect_ren_sam_rmdup2$JCInputID %nin% c("1762", "2085"),]
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
##share this list with Denia for generating PID minus 3105 which fails QC2 ad hence must be left out. 3105 is addressed
##in the SKAT script
write.table(rect_ren_sam_all, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/ISKS_RISC_LIONS_final_freeze.tsv",
            sep = "\t", row.names = F, quote = F)
