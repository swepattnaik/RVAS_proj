rnd3_var <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round3/all_isksmgrb_combset2020_variants_filt_all_fields_rnd3.tsv",
                       sep = "\t", header = T, stringsAsFactors = F)
`%nin%` = Negate(`%in%`)
table(rnd3_var$GT)
table(rnd3_var[rnd3_var$GT == 2, ]$auto_call)
table(rnd3_var[rnd3_var$GT == 2 & rnd3_var$auto_call %in% c("C5", "C4"), ]$gene_symbol)
boxplot(rnd3_var[rnd3_var$GT == 2 & rnd3_var$auto_call %in% c("C5", "C4"), ]$VAF)
dim(rnd3_var[rnd3_var$GT == 2 & rnd3_var$auto_call %in% c("C5", "C4") & rnd3_var$VAF == 1, c(1:3,9,11,82,125:126)])
head(rnd3_var[rnd3_var$GT == 2 & rnd3_var$auto_call %in% c("C5", "C4") & rnd3_var$VAF < 0.9,c(1:3,9,11,82,125:126)])
length(unique(rnd3_var[rnd3_var$GT == 2 & rnd3_var$auto_call %in% c("C5", "C4"), ]$gene_symbol))
rnd3_var[rnd3_var$GT == 2 & rnd3_var$auto_call %in% "C5" & rnd3_var$gene_symbol %in% "SERPINA7", c(1:3,9,11,82,125:126)]
rnd3_var[rnd3_var$GT == 2 & rnd3_var$auto_call %in% "C4" & rnd3_var$gene_symbol %in% "TPRG1L", c(1:3,9,11,82,125:126)]
head(rnd3_var[rnd3_var$GT == 2 & rnd3_var$auto_call %in% "C4" & rnd3_var$vep_consequence %nin% "frameshift_variant", c(1:3,9,11,82,125:126)])

dim(rnd3_var[rnd3_var$SAMPLE %in% "2371",])
dim(rnd3_var[rnd3_var$SAMPLE %in% "2371" & rnd3_var$vep_CLIN_SIG %in% "pathogenic",])
dim(rnd3_var[rnd3_var$gene_symbol %in% "ABCB6",])

#####
table(rnd3_var[rnd3_var$GT == 2, ]$auto_call)
#check with gene_sub
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
DT_genes <- read.table("~/RVAS/DT_genes_sep2018.txt", sep = "", stringsAsFactors = F) ##includes POT1 interactors
DT_ClueGO <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/ClueGO_list.txt", sep = "", stringsAsFactors = F)
DT_ClueGO_plus <- c("BUB1", "CEP57", "CDC20", "TRIP13")
DT_new_add <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/add_PID_May31.txt", sep = "", stringsAsFactors = F)
DT_new_add_June <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/add_PID_June1.txt", sep = "", stringsAsFactors = F)
#gene_sub <- unique(c(top300_SKAT$gene[1:300], cgc_genes$Gene.Symbol, sheltrin_complex, Telo_extension, DT_genes$V1, DT_ClueGO$V1, Sheltrin_comp_extn))
gene_sub <- unique(c(cgc_genes$Gene.Symbol, sheltrin_complex, Telo_extension, DT_genes$V1, DT_ClueGO$V1, DT_ClueGO_plus, DT_new_add$V1, DT_new_add_June$V1, Sheltrin_comp_extn))

table(rnd3_var[rnd3_var$GT == 2 & rnd3_var$gene_symbol %in% gene_sub & rnd3_var$auto_call %in% c("C4", "C5"), ])

rnd3_var[rnd3_var$GT == 2 & rnd3_var$gene_symbol %in% gene_sub & rnd3_var$auto_call %in% c("C4", "C5"), c(1:3,9,11,82,125:126)]

dim(rnd3_var[rnd3_var$GT == 2 & rnd3_var$auto_call %in% c("C4", "C5"), c(1:6,8:9,11,15,17:18,20,23:25,30,53,57,59:68,72,125:132)])
var_C4C5_rnd3_hom_all <- rnd3_var[rnd3_var$GT == 2 & rnd3_var$auto_call %in% c("C4", "C5"), c(1:6,8:9,11,15,17:18,20,23:25,30,53,57,59:68,72,125:132)]
write.table(var_C4C5_rnd3_hom_all, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round3/var_C4C5_rnd3_hom_all.tsv",
            sep = "\t", quote = F, row.names = F)

dim(rnd3_var[rnd3_var$GT == 2 & rnd3_var$gene_symbol %in% gene_sub & rnd3_var$auto_call %in% c("C4", "C5"), c(1:6,8:9,11,15,17:18,20,23:25,30,53,57,59:68,72,125:132)])
var_C4C5_rnd3_hom_PID_genes <- rnd3_var[rnd3_var$GT == 2 & rnd3_var$gene_symbol %in% gene_sub & rnd3_var$auto_call %in% c("C4", "C5"), c(1:6,8:9,11,15,17:18,20,23:25,30,53,57,59:68,72,125:132)]
write.table(var_C4C5_rnd3_hom_PID_genes, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round3/var_C4C5_rnd3_hom_PID_genes.tsv",
            sep = "\t", quote = F, row.names = F)

##PTEN
dim(rnd3_var[rnd3_var$gene_symbol %in% "PTEN",])
table(rnd3_var[rnd3_var$gene_symbol %in% "PTEN",]$auto_call)
rnd3_var[rnd3_var$gene_symbol %in% "PTEN",c(1:3,9,11,82,125:126)]
# rnd2_var <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/all_isksmgrb_combset2020_variants_filt_all_fields_rmdup.tsv",
#                        sep = "\t", header = T, stringsAsFactors = F)
# t1 <- rnd2_var[rnd2_var$SAMPLE %in% "2371",]
# dim(rnd2_var[rnd2_var$SAMPLE %in% "2371",])
# dim(rnd2_var[rnd2_var$SAMPLE %in% "2371" & rnd2_var$vep_CLIN_SIG %in% "pathogenic",])


##PID file

PID_jun1 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/PID_combset2020_CGC_skatBin_repstress_potint_predNFE_rmdup_maf_clueGO_Apr292020_June1.tsv",
sep = "\t", header = T, stringsAsFactors = F)

##PID gene list does not contain ABCB6, the C5 variant containing gene
##SEAVE captured HNF1A which is part of gene_sub list