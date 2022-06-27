

library(dplyr)
library(stringr)
library(readxl)

`%nin%` = Negate(`%in%`)
####PID generate script
fil_tab_PID = read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CGC_topSKATBin_ISKS_MRGB_repstress_potint_mito_chkpt_centrosome_VARIANTS_filt_combset2020_clueGOplus_Sep19.tsv",
sep = "\t", header = T, stringsAsFactors = F)


fil_tab_PID_C45 = fil_tab_PID[fil_tab_PID$auto_call %in% c("C4", "C5", "C4_ar", "C4_ar") & fil_tab_PID$is_case == 1,]
fil_tab_PID_C3 = fil_tab_PID[fil_tab_PID$auto_call %nin% c("C4", "C5", "C4_ar", "C4_ar") & fil_tab_PID$is_case == 1,]
 
write.table(fil_tab_PID_C45, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CGC_topSKATBin_ISKS_repstress_potint_mito_chkpt_centrosome_C4_C5_VARIANTS_filt_combset2020_clueGOplus_Feb17.tsv",
            sep = "\t", row.names = F, quote = F)
write.table(fil_tab_PID_C3, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CGC_topSKATBin_ISKS_repstress_potint_mito_chkpt_centrosome_C3_VARIANTS_filt_combset2020_clueGOplus_Feb17.tsv",
            sep = "\t", row.names = F, quote = F)

##########
##latest PID file

##For LOH analysis
LOH_genes = c("POT1", "TERF1", "TINF2", "TERF2IP", "SMARCAL1", 
              "TIMELESS", "STAG3", "CEP63", "CEP72", "HAUS4", "HAUS5", "TP53", 
              "NF1", "LZTR1", "SDHA", "SDHB", "SDHC", "SDHD", "BRCA1", "BRCA2", 
              "PALB2", "CTC1")

fil_tab_PID_LOH_genes =  fil_tab_PID[fil_tab_PID$gene_symbol %in% LOH_genes,]

##subset GIST/MPNST samples that have tumour data
path = "/home/shu/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/review_2"
GIST_MPNST = read.delim(paste(path, "/", "LOH_GIST_MPNST_samp.txt", sep = ""), header = F, stringsAsFactors = F, sep = "")

####PID generate script with PMN
fil_tab_PID_2020 = read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/PID_combset2020_CGC_skatBin_repstress_potint_mito_chkpt_centrosome_predNFE_clueGO_Sep192020_AD_addC4C5.tsv",
                              sep = "\t", header = T, stringsAsFactors = F)
PMN_ID = fil_tab_PID_2020[match(GIST_MPNST$V1,fil_tab_PID_2020$pid),3]
PMN_ID = PMN_ID[!is.na(PMN_ID)]
fil_tab_PID_LOH_genes_gist_mpnst = fil_tab_PID_LOH_genes[fil_tab_PID_LOH_genes$SAMPLE %in% PMN_ID,]
fil_tab_PID_LOH_genes_gist_mpnst$PID = fil_tab_PID_2020[match(fil_tab_PID_LOH_genes_gist_mpnst$SAMPLE,fil_tab_PID_2020$pmn),1]

write.table(fil_tab_PID_LOH_genes_gist_mpnst[,c(1:9,11, 127:128,133)], paste(path, "/", "for_LOH_GIST_MPNST_var.tsv", sep = ""), sep = "\t", row.names = F, quote = F)
