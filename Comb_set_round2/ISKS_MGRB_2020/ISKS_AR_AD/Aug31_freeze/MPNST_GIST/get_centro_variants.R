##CEP63, CEP72, HAUS4, HAUS5, PCM1 variants

SKAT_inp <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv", 
                       sep = "\t", header = T, stringsAsFactors = F)
SKAT_inp_centro <- SKAT_inp[SKAT_inp$gene_symbol %in% c("CEP63", "CEP72", "HAUS4", "HAUS5", "PCM1"),]
SKAT_inp_centro_isks <- SKAT_inp_centro[SKAT_inp_centro$set %in% "ISKS_AR_AD",c(1:6,9,11,126:128,130)]
write.table(SKAT_inp_centro_isks, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/ISKS_centro_C345.tsv",
            sep = "\t", row.names = F, quote = F)
