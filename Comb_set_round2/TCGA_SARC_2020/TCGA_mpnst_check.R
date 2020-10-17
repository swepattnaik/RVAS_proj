
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

`%nin%` = Negate(`%in%`)
Ex_tab_sarc <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/TCGA_SARC_AR_AD/all_tcga_sarc_combset2020_variants_AR_AD_all_fields_clingene_rnd3.tsv", sep = "\t", header = T, stringsAsFactors = F)
sarc_pheno <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/TCGA_SARC_AR_AD/TCGA_SARC_age_sex_info.tsv",
                         sep = "\t", header = T, stringsAsFactors = F)
library(readxl)

tcga_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/TCGA_pheno_sarc_S0092867417312035_mmc1.xlsx", sheet = 2)
tcga_pheno <- as.data.frame(tcga_pheno)
hdr <- as.character(tcga_pheno[1,])
colnames(tcga_pheno) <- hdr
tcga_pheno <- tcga_pheno[-1,]
tcga_pheno_mpnst_ID <- tcga_pheno[tcga_pheno$`short histo` %in% "MPNST",]$`TCGA barcode`
tcga_pheno_mpnst_ID_mod <- unlist(lapply(strsplit(tcga_pheno_mpnst_ID,split = "-"), function(x)paste(x[1],x[2], x[3], sep = "-")))

sarc_pheno_mpnst_ID <- sarc_pheno[match(sarc_pheno$Patient.ID, tcga_pheno_mpnst_ID_mod),2]
sarc_pheno_mpnst_ID <- sarc_pheno_mpnst_ID[!is.na(sarc_pheno_mpnst_ID)]

Extab_tcga_mpnst_genes <- unique(Ex_tab_sarc[Ex_tab_sarc$SAMPLE %in% sarc_pheno_mpnst_ID,]$gene_symbol)

Ex_tab_sarc[Ex_tab_sarc$SAMPLE %in% sarc_pheno_mpnst_ID & Ex_tab_sarc$gene_symbol %in% c("CEP85L","NUP155", "NUPL2", "HAUS4", "HAUS5",
                                                                                         "CEP63", "CEP72", "PCM1", "MZT1", "SSNA1", "MACC1", "NF1"),c(1:3,9,11,127:128)]

##unfiltered sarc_tcga

Ex_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/TCGA_SARC_AR_AD/all_SARC_TCGA_gnomad05.tsv", sep = "\t",
                     header = T, stringsAsFactors = F)
Extab_unfilt_tcga_mpnst_genes <- unique(Ex_tab[Ex_tab$SAMPLE %in% sarc_pheno_mpnst_ID,]$gene_symbol)

pcm1_var <- Ex_tab[Ex_tab$SAMPLE %in% sarc_pheno_mpnst_ID & Ex_tab$gene_symbol %in% "PCM1",c(1:3,9,11,30,31)]

Ex_tab[Ex_tab$SAMPLE %in% sarc_pheno_mpnst_ID & Ex_tab$gene_symbol %in% c("CEP85L","NUP155", "NUPL2", "HAUS4", "HAUS5",
                                                                          "CEP63", "CEP72", "PCM1", "MZT1", "SSNA1", "MACC1", "NF1"),c(1:3,9,11,30,31)]
##also check after running all regular filters (VAF, depth etc)