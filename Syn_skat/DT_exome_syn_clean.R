##Synonymous variant effect processing
##Uses MAF for scoring
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
Ex_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/comb_isks_risc_latest05Sept_rect1.tsv", sep = "\t",
                     header = T, stringsAsFactors = F)
Ex_tab <- Ex_tab[!is.na(Ex_tab$SAMPLE),]
Ex_tab$SAMPLE <- gsub(".variant.*$", "", Ex_tab$SAMPLE) ##check this in the next run
#remove variants ending with Asterisk; these donot have vep annotation
Ex_tab <- Ex_tab[!is.na(Ex_tab$gene_symbol),]

Ex_tab <- Ex_tab[!is.na(Ex_tab$DP),] ##There are records with NA in DP attribute ; file: sept05_rect


Ex_tab <- Ex_tab[!(Ex_tab$VAF %in% "NaN"),]
class(Ex_tab$VAF) <- "numeric" ##remove headers from shards in the next run
Ex_tab <- Ex_tab[!is.na(Ex_tab$VAF),] ##some VAFs are NAs in sept05_rect file


DT_gene_tab <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/DT_gene_tab_all.rds")
load("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/ASPREE_can.RData")

#DT_gene_tab_ISKS_MGRB_nCH <-DT_gene_tab[DT_gene_tab$cohort %in% c("ISKS", "MGRB no CH"),]

MGRB_noCH <- unique(DT_gene_tab[(DT_gene_tab$cohort) %in% "MGRB no CH",]$SAMPLE)

MGRB_all <- unique(DT_gene_tab[(DT_gene_tab$cohort) %in% c("MGRB no CH", "MGRB"),]$SAMPLE)

Ex_tab$is_CH <- ifelse(Ex_tab$SAMPLE %in% MGRB_noCH, 1, 0)


#Ex_tab1 <- Ex_tab[Ex_tab$VAF >= 0.25 & Ex_tab$VAF < 1,] ##OR use the command below to not lose highly significant homozygous variants
Ex_tab1 <- Ex_tab[as.numeric(Ex_tab$VAF) >= 0.25 & as.numeric(Ex_tab$DP) >= 10,] ##Depth filter
Ex_tab1$gnomad_AF_NFE <- gsub("\\[", "", Ex_tab1$gnomad_AF_NFE)
Ex_tab1$gnomad_AF_NFE <- gsub("\\]", "", Ex_tab1$gnomad_AF_NFE)
class(Ex_tab1$gnomad_AF_NFE) <- "numeric"
Ex_tab1$gnomad_AF_NFE <- ifelse(is.na(Ex_tab1$gnomad_AF_NFE), 0, Ex_tab1$gnomad_AF_NFE)

#Ex_tab_filt1 <- Ex_tab1[Ex_tab1$gnomad_AF_NFE <= 0.0001, ] ##check this; loses some variants from Aug12 on Sept05 run(don't use)

##remove SAMPLE CR57
Ex_tab_filt1 <- Ex_tab1[!(Ex_tab1$SAMPLE %in% "CR57"),]

#toMatch <- c("synonymous_variant")

Ex_tab_filt2 <- Ex_tab_filt1[grepl("synonymous_variant", Ex_tab_filt1$vep_consequence),]

Ex_tab_filt2$is_MGRB <- ifelse(Ex_tab_filt2$SAMPLE %in% MGRB_all, 1, 0)

Ex_tab_filt2$is_ASPC <- ifelse(Ex_tab_filt2$SAMPLE %in% ASPREE_can$samp, 1, 0)

##filter out MGRB samples that have CH

Ex_tab_noCH_filt3 <- Ex_tab_filt2[!(Ex_tab_filt2$is_CH == 0 & Ex_tab_filt2$is_MGRB == 1),]
Ex_tab_noCH_filt3 <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$is_ASPC == 0,]

##select synonymous variants based on delta in preferred vs mutated.. 
##select variants if abs(preferred vs mutated) >= 2
codon_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Syn/codon_usage_table_genscript.xlsx.txt",
                        sep = "\t", header = T, stringsAsFactors = F)
vep_codon <- do.call("rbind.data.frame", strsplit(Ex_tab_noCH_filt3$vep_Codons, "/")) 
colnames(vep_codon) <- c("ref_cod", "alt_cod")
vep_codon <- as.data.frame(apply(vep_codon,2,toupper), stringsAsFactors = F)
vep_codon$ref_frac <- codon_tab[match(vep_codon$ref_cod, codon_tab$Triplet),3]
vep_codon$alt_frac <- codon_tab[match(vep_codon$alt_cod, codon_tab$Triplet),3]
vep_codon$delta <- abs(vep_codon$ref_frac - vep_codon$alt_frac)/vep_codon$ref_frac
Ex_tab_noCH_filt4 <- cbind.data.frame(Ex_tab_noCH_filt3, vep_codon)
Ex_tab_noCH_filt4 <- Ex_tab_noCH_filt4[Ex_tab_noCH_filt4$delta >= 0.5,]

saveRDS(Ex_tab_noCH_filt4, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Syn/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_05Sept_rect_ASP_syn.rds", compress = T)
