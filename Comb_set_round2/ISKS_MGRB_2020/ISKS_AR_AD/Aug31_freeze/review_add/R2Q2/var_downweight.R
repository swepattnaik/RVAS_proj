##R2Q2
##why were the variants in the last 5% of the protein sequence scaled down
##demonstrate analytically the choice for scoring ESS as C3

`%nin%` = Negate(`%in%`)

##start gridsearch
Ex_tab_sarc <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/R2Q8/Ex_tab_sarc_genes_rnd3_freeze.tsv",
                          sep = "\t", header = T, stringsAsFactors = F)
#variants in last 5% of the protein sequence
table(Ex_tab_sarc[Ex_tab_sarc$last_five_perc == 1 & Ex_tab_sarc$vep_consequence %in% "stop_gained",]$clinvar_Clinical_Significance)
table(Ex_tab_sarc[Ex_tab_sarc$last_five_perc == 1 & Ex_tab_sarc$vep_consequence %in% "stop_gained",]$set)
table(Ex_tab_sarc[Ex_tab_sarc$vep_consequence %in% "stop_gained",]$clinvar_Clinical_Significance)
table(Ex_tab_sarc[Ex_tab_sarc$vep_consequence %in% "stop_gained",]$set)
