##NC variants in Sarcoma risk genes
`%nin%` = Negate(`%in%`)

`%nin%` = Negate(`%in%`)
setwd("~/RVAS/ISKS_MGRB_NC/")
print("Reading annotated files")
files = grep("mgrb_complexes_pos_neg", dir("~/RVAS/ISKS_MGRB_NC/"), value = T)
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
saveRDS(comb_list2, "~/RVAS/ISKS_MGRB_NC/results/filt_MGRB_dflist.rds", compress = T)
mgrb_nc_var_filt3 <- do.call("rbind.data.frame", comb_list2)

# print("Filtering")
# ##gnomAD_NFE_AF filter
# mgrb_nc_var$gnomAD_NFE_AF <- as.numeric(mgrb_nc_var$gnomAD_NFE_AF)
# mgrb_nc_var$gnomAD_NFE_AF <- ifelse(is.na(mgrb_nc_var$gnomAD_NFE_AF), 0, mgrb_nc_var$gnomAD_NFE_AF)
# mgrb_nc_var_filt1 <- mgrb_nc_var[mgrb_nc_var$gnomAD_NFE_AF <= 0.0002,]
# colnames(mgrb_nc_var_filt1)[1:10] <- c("CHROM", "POS", "ID", "REF", "ALT", "GT", "GQ", "PL", "DP", "AD")
# mgrb_nc_var_filt1$AD <- as.numeric(unlist(lapply(strsplit(mgrb_nc_var_filt1$AD, split = ","), function(x)x[2])))
# mgrb_nc_var_filt2 <- mgrb_nc_var_filt1[mgrb_nc_var_filt1$AD > 0, ]
# mgrb_nc_var_filt2 <- mgrb_nc_var_filt2[mgrb_nc_var_filt2$DP >= 10,]
# mgrb_nc_var_filt2$VAF <- mgrb_nc_var_filt2$AD/mgrb_nc_var_filt2$DP
# ##somatic filter
# mgrb_nc_var_filt2 = mgrb_nc_var_filt2[mgrb_nc_var_filt2$VAF >= 0.35,]
# ##GQ filter
# mgrb_nc_var_filt2 = mgrb_nc_var_filt2[mgrb_nc_var_filt2$GQ >= 80,]
# 
# 
# library(tidyr)
# library(dplyr)
# 
# t1 <- mgrb_nc_var_filt2 %>% group_by(sel_ID) %>% summarise(n = n()) 
# 
# t2 <- as.data.frame(t1)
# 
# t3 <- t2[t2$n < 4,]
# 
# ##remove variants with a frequency of greater than 3 in ISKS (IntraCohort MAF 3/(1*1644))
# mgrb_nc_var_filt3 = mgrb_nc_var_filt2[mgrb_nc_var_filt2$sel_ID %in% t3$sel_ID,]

print("Removal of Duplicated and QC failed samples")
##remove duplicates and rename samples

QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv", 
                      header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %in% "TRUE",]
dup_samp <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/fin_samp_dup_drop.txt", sep = "", header = T,
                       stringsAsFactors = F)
dup_samp <- dup_samp$x[grepl("^[ABZ]", dup_samp$x)]

mgrb_nc_var_filt4 = mgrb_nc_var_filt3[mgrb_nc_var_filt3$sid %in% QC2_dat_pass$new_sampleid,]
mgrb_nc_var_filt4 = mgrb_nc_var_filt4[mgrb_nc_var_filt4$sid %nin% dup_samp,]

##GQ filter
#mgrb_nc_var_filt4 = mgrb_nc_var_filt4[mgrb_nc_var_filt4$GQ >= 80,] 
write.table(mgrb_nc_var_filt4, "~/RVAS/ISKS_MGRB_NC/results/mgrb_nc_var_filt_rm_dup.tsv", sep = "\t", row.names = F,
            quote = F)
##Add NC scores
print("Scoring Variants")
NC_scores = read.delim("~/RVAS/ISKS_MGRB_NC/Noncoding_ncScore.txt", sep = "\t", header = F, stringsAsFactors = F)

NC_scores$sel1 = gsub("^chr", "",paste(NC_scores$V1, NC_scores$V2, sep = ":"))
#NC_scores$sel2 = gsub("^chr", "",paste(NC_scores$V1, NC_scores$V3, sep = ":"))
mgrb_nc_var_filt4$NC_scores = NC_scores[match(mgrb_nc_var_filt4$Location, NC_scores$sel1), 4]
write.table(mgrb_nc_var_filt4, "~/RVAS/ISKS_MGRB_NC/results/mgrb_nc_var_filt_rm_dup_NC_scored.tsv", sep = "\t", row.names = F,
            quote = F)
print("Done!")