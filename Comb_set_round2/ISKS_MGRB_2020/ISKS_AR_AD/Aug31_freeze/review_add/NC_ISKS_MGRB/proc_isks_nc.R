##NC variants in Sarcoma risk genes
`%nin%` = Negate(`%in%`)
setwd("~/RVAS/ISKS_MGRB_NC/")
print("Reading annotated files")
files = grep("AD_filt_Comb_vep_tab_isks", dir("~/RVAS/ISKS_MGRB_NC/"), value = T)
comb_list = list()
for(i in 1:length(files)){
  f1 = readRDS(files[i])
  comb_list[[i]] = f1
}
library(dplyr)

fix_data_frame = function(x) {
  x[is.na(x[,58]),58] = 0
  x[is.na(x[,64]),64] = 0
  # What ever else you want to do
  return(x)
}

filt_data_frame = function(df) {
  df <- df[df[,9] >= 10,]
  df$VAF <- df[,10]/df[,9]
  ##somatic filter
  df = df[df$VAF >= 0.35,]
  df = df[df[,64] <= 0.0002,]
  return(df)
}

library(tidyr)
library(dplyr)
coh_MAF_filt = function(df){
tt1 <- df %>% group_by(ID) %>% summarise(n = n()) 

tt2 <- as.data.frame(tt1)

tt3 <- tt2[tt2$n < 4,]

df1 = df[df[,3] %in% tt3[,1],]
return(df1)

}

print("Filtering")
comb_list1 <- lapply(comb_list, function(df) {df[c("gnomAD_AF", "gnomAD_NFE_AF")] <- lapply(df[c("gnomAD_AF", "gnomAD_NFE_AF")], as.numeric); df})
comb_list1 = lapply(comb_list1, function(x)fix_data_frame(x))
comb_list1 = lapply(comb_list1, function(x)filt_data_frame(x))
comb_list2 = lapply(comb_list1, function(x)coh_MAF_filt(x))
saveRDS(comb_list2, "~/RVAS/ISKS_MGRB_NC/results/filt_ISKS_dflist.rds", compress = T)
#isks_nc_var <- readRDS("~/RVAS/ISKS_MGRB_NC/vep_annot_comb.rds") 
isks_nc_var_filt3 <- do.call("rbind.data.frame", comb_list2)

# ##gnomAD_NFE_AF filter
# isks_nc_var$gnomAD_NFE_AF <- as.numeric(isks_nc_var$gnomAD_NFE_AF)
# isks_nc_var$gnomAD_NFE_AF <- ifelse(is.na(isks_nc_var$gnomAD_NFE_AF), 0, isks_nc_var$gnomAD_NFE_AF)
# isks_nc_var_filt1 <- isks_nc_var[isks_nc_var$gnomAD_NFE_AF <= 0.0002,]
# # colnames(isks_nc_var_filt1)[1:10] <- c("CHROM", "POS", "ID", "REF", "ALT", "GT", "GQ", "PL", "DP", "AD")
# # isks_nc_var_filt1$AD <- as.numeric(unlist(lapply(strsplit(isks_nc_var_filt1$AD, split = ","), function(x)x[2])))
# # isks_nc_var_filt2 <- isks_nc_var_filt1[isks_nc_var_filt1$AD > 0, ]
# isks_nc_var_filt2 <- isks_nc_var_filt2[isks_nc_var_filt2$DP >= 10,]
# isks_nc_var_filt2$VAF <- isks_nc_var_filt2$AD/isks_nc_var_filt2$DP
# ##somatic filter
# isks_nc_var_filt2 = isks_nc_var_filt2[isks_nc_var_filt2$VAF >= 0.35,]
# 
# 
# 
# library(tidyr)
# library(dplyr)
# 
# t1 <- isks_nc_var_filt2 %>% group_by(sel_ID) %>% summarise(n = n()) 
# 
# t2 <- as.data.frame(t1)
# 
# t3 <- t2[t2$n < 4,]
# 
# ##remove variants with a frequency of greater than 3 in ISKS (IntraCohort MAF 3/(1*1644))
# isks_nc_var_filt3 = isks_nc_var_filt2[isks_nc_var_filt2$sel_ID %in% t3$sel_ID,]

print("Removal of Duplicated and QC failed samples")
##remove duplicates and rename samples
rect_sam_dat <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/ISKS_RISC_LIONS_final_freeze.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
isks_nc_var_filt3$rect_sam <- rect_sam_dat[match(isks_nc_var_filt3$sid,rect_sam_dat$JCInputID),9]
isks_nc_var_filt4 <- isks_nc_var_filt3[!is.na(isks_nc_var_filt3$rect_sam),]
QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv", 
                      header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %in% "TRUE",]

isks_nc_var_filt4 = isks_nc_var_filt4[isks_nc_var_filt4$rect_sam %in% QC2_dat_pass$new_sampleid,]
##GQ filter
isks_nc_var_filt4 = isks_nc_var_filt4[isks_nc_var_filt4$GQ >= 80,] 
write.table(isks_nc_var_filt4, "~/RVAS/ISKS_MGRB_NC/results/isks_nc_var_filt_rm_dup.tsv", sep = "\t", row.names = F,
            quote = F)
##Add NC scores
print("Scoring Variants")
NC_scores = read.delim("~/RVAS/ISKS_MGRB_NC/Noncoding_ncScore.txt", sep = "\t", header = F, stringsAsFactors = F)

NC_scores$sel1 = gsub("^chr", "",paste(NC_scores$V1, NC_scores$V2, sep = ":"))
#NC_scores$sel2 = gsub("^chr", "",paste(NC_scores$V1, NC_scores$V3, sep = ":"))
isks_nc_var_filt4$NC_scores = NC_scores[match(isks_nc_var_filt4$Location, NC_scores$sel1), 4]
write.table(isks_nc_var_filt4, "~/RVAS/ISKS_MGRB_NC/results/isks_nc_var_filt_rm_dup_NC_scored.tsv", sep = "\t", row.names = F,
            quote = F)
