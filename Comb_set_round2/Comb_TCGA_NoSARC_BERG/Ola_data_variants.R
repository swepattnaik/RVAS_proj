##Ola data
#Script:/home/shu/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/Comb_TCGA_NoSARC_BERG/comb_clingene_TCGA_NoSARC_BERG_MGRB_var4SKAT.R:
#line 156

`%nin%` = Negate(`%in%`)
comb_skat_inp <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/all_tcga_nosarc_berg_mgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_all_fields_rnd3_Aug31.tsv", 
           sep = "\t", header = T, stringsAsFactors = F)
##Extract NoSarc and TCGA for Ola
comb_skat_inp_ola <- comb_skat_inp[comb_skat_inp$set %in% c("NOSARC_nov18_AR_AD", "SARC_BERG_AR_AD"),]

comb_skat_inp_ola <- comb_skat_inp_ola[comb_skat_inp_ola$EigenPhred >= 5.6,]

##exclude cohort_MAF
comb_skat_inp_ola <- comb_skat_inp_ola[,c(1:9,22,20:21,23:24,30:31,53,57,76,83,112,127,130)]

write.table(comb_skat_inp_ola, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/Germline_variants_NoSARC_BERG.tsv",
            row.names = F, quote = F, sep = "\t")
