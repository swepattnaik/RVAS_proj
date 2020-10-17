##preprocess variant input file
##1. remove duplicate samples
##2. remove related individuals

fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/all_isksmgrb_combset2020_variants_filt_all_fields.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)

#remove duplicates
dup_samp <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/fin_samp_dup_drop.txt", sep = "", header = T,
                       stringsAsFactors = F)
`%nin%` = Negate(`%in%`)
fil_tab_rmdup <- fil_tab[fil_tab$SAMPLE %nin% dup_samp$x, ]
write.table(fil_tab_rmdup, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/all_isksmgrb_combset2020_variants_filt_all_fields_rmdup.tsv",
            sep = "\t", row.names = F, quote = F)
##remove degree2 relatives (do this if needed)
deg2_isks <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/dup_pairs_and_deg2_isks.tsv", sep = "\t",
                        header = T, stringsAsFactors = F)
deg2_mgrb <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/dup_pairs_and_deg2_mgrb.tsv", sep = "\t",
                        header = T, stringsAsFactors = F)
deg2_isks <- deg2_isks[deg2_isks$degree > 0,]
deg2_mgrb <- deg2_mgrb[deg2_mgrb$degree > 0,]
Ex_samp_id <- unique(fil_tab_rmdup$SAMPLE)
Ex_samp_id <- Ex_samp_id[Ex_samp_id %nin% deg2_isks$name1]
Ex_samp_id <- Ex_samp_id[Ex_samp_id %nin% deg2_mgrb$name1]
fil_tab_rmdeg2 <- fil_tab_rmdup[fil_tab_rmdup$SAMPLE %in% Ex_samp_id,]
write.table(fil_tab_rmdeg2, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/all_isksmgrb_combset2020_variants_filt_all_fields_rmdup_rmdeg2.tsv",
            sep = "\t", row.names = F, quote = F)

######remove degree 2 relatives: workflow
##check ~/RVAS/comb_set_2020/PC_relate/code/pcrelate_sample_drop_check.R for details
##scratch work below inspired from pcrelate_sample_drop_check.R
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

library(SNPRelate)
library(gdsfmt)
library(GENESIS)
library(GWASTools)
pcrelate_results <- readRDS("~/RVAS/comb_set_2020/PC_relate/ldpruned/MGRBISKS.merged.combset2020.SNPtier12.match.vqsr.minrep.ldpruned.pcrelate.rds")
######MGRB
kinship_MGRB = pcrelate_results$kinship[grepl("^[ABZ]", pcrelate_results$kinship$name1) & grepl("^[ABZ]", pcrelate_results$kinship$name2),]
kinship_MGRB$degree = "5+"
for (degree_i in 0:4)
{
  kinship_MGRB$degree[kinship_MGRB$kin >= 2^(-(degree_i+3/2)) & kinship_MGRB$kin < 2^(-(degree_i+1/2))] = as.character(degree_i)
}
kinship_MGRB$degree = ordered(as.character(kinship_MGRB$degree), levels = c(as.character(0:4), "5+"))

######ISKS
kinship_ISKS = pcrelate_results$kinship[!(grepl("^[ABZ]", pcrelate_results$kinship$name1)) & !(grepl("^[ABZ]", pcrelate_results$kinship$name2)),]
##ISKS
kinship_ISKS$degree = "5+"
for (degree_i in 0:4)
{
  kinship_ISKS$degree[kinship_ISKS$kin >= 2^(-(degree_i+3/2)) & kinship_ISKS$kin < 2^(-(degree_i+1/2))] = as.character(degree_i)
}
kinship_ISKS$degree = ordered(as.character(kinship_ISKS$degree), levels = c(as.character(0:4), "5+"))

##combine MGRB and ISKS

dim(kinship_MGRB[kinship_MGRB$kin >= 0.2 & kinship_MGRB$degree <= 2, ])
dim(kinship_ISKS[kinship_ISKS$kin >= 0.2 & kinship_ISKS$degree <= 2, ])

dim(kinship_MGRB[kinship_MGRB$degree == 1, ])
dim(kinship_ISKS[kinship_ISKS$degree == 1, ])
##removal of degree1 relatives leads to a loss of C4 variant in POT1; The power loss is not much.

dim(fil_tab[fil_tab$gene_symbol %in% c("TP53","NF1", "TINF2", "POT1", "TERF1", "TERF2IP", "BRCA1", "BRCA2", "SDHB", "SDHA", "SDHD", "EXT1", "EXT2") &
              fil_tab$SAMPLE %in% kinship_ISKS[kinship_ISKS$degree == 1, ]$name1,c(1:2,9, 125:126)])

dim(fil_tab[fil_tab$gene_symbol %in% c("TP53","NF1", "TINF2", "POT1", "TERF1", "TERF2IP", "BRCA1", "BRCA2", "SDHB", "SDHA", "SDHD", "EXT1", "EXT2") &
              fil_tab$SAMPLE %in% kinship_MGRB[kinship_MGRB$degree == 1, ]$name2,c(1:2,9, 125:126)])
##check this out