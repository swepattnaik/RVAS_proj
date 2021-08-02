#.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)

Ex_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/MGRB/all_shards_MGRB_combset2020_gnomadAF05.tsv", sep = "\t",
                     header = T, stringsAsFactors = F)

#rm(isks_tsv, mgrb_tsv)
Ex_tab <- Ex_tab[!is.na(Ex_tab$SAMPLE),]
#Ex_tab$SAMPLE <- gsub(".variant.*$", "", Ex_tab$SAMPLE) ##check this in the next run
#remove variants ending with Asterisk; these donot have vep annotation
Ex_tab <- Ex_tab[!is.na(Ex_tab$gene_symbol),]

Ex_tab <- Ex_tab[!is.na(Ex_tab$DP),] ##There are records with NA in DP attribute ; file: sept05_rect

#remove duplicates
dup_samp <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/fin_samp_dup_drop.txt", sep = "", header = T,
                       stringsAsFactors = F)
`%nin%` = Negate(`%in%`)
Ex_tab <- Ex_tab[Ex_tab$SAMPLE %nin% dup_samp$x,]

Ex_tab <- Ex_tab[!(Ex_tab$VAF %in% "NaN"),]
class(Ex_tab$VAF) <- "numeric" ##remove headers from shards in the next run
Ex_tab <- Ex_tab[!is.na(Ex_tab$VAF),] ##some VAFs are NAs in sept05_rect file
Ex_tab1 <- Ex_tab[as.numeric(Ex_tab$VAF) >= 0.35 & as.numeric(Ex_tab$DP) >= 10,]


##Misannotation of MUTYH to HPDL.. change that
#MUTYH hg19 coordinates: chr1: 45,794,835-45,806,142
#chk_gene_ind <- ifelse(Ex_tab$gene_symbol %nin% Ex_tab$clinvar_Gene, 1, 0)
#chk_gene_diff <- unique(Ex_tab$gene_symbol[which(chk_gene_ind == 1)])
#Ex_tab_chk <- Ex_tab[Ex_tab$gene_symbol ]
Ex_tab_clinvar_gene <- Ex_tab1[!is.na(Ex_tab1$clinvar_Gene),]
Ex_tab_clinvar_gene <- Ex_tab_clinvar_gene[Ex_tab_clinvar_gene$VAF >= 0.35 & Ex_tab_clinvar_gene$DP >= 10,]

chk_gene_ind <- ifelse(Ex_tab_clinvar_gene$gene_symbol %nin% Ex_tab_clinvar_gene$clinvar_Gene, 1, 0)
chk_gene_diff_sym <- Ex_tab_clinvar_gene$gene_symbol[which(chk_gene_ind == 1)]
chk_gene_diff_clin <- Ex_tab_clinvar_gene$clinvar_Gene[which(chk_gene_ind == 1)]
gene_map_rescue <- cbind.data.frame(chk_gene_diff_sym, chk_gene_diff_clin)
gene_map_rescue <- unique(gene_map_rescue)
#gene_sub[gene_sub %in% gene_map_rescue_path$clinvar_Gene]
##feature to select: vep_biotype(73), vep_VARIANT_CLASS c(1:9,11,17,22:24,30,32,33,53,73,87)

Ex_tab_Pathogenic <- Ex_tab_clinvar_gene[grep("pathogenic", Ex_tab_clinvar_gene$clinvar_Clinical_Significance, ignore.case = T),]
rm_patho <- c("Conflicting", "protective")
fs_var <- c("indel", "insertion", "sequence_alteration", "SNV")
#grepl(paste(rm_patho,collapse="|"), 
#      Ex_tab_Pathogenic$clinvar_Clinical_Significance)
Ex_tab_Pathogenic_C5 <- Ex_tab_Pathogenic[grep(paste(rm_patho,collapse="|"), 
                                               Ex_tab_Pathogenic$clinvar_Clinical_Significance, invert = T),]
Ex_tab_Pathogenic_C4C5 <- Ex_tab_Pathogenic_C5[Ex_tab_Pathogenic_C5$vep_VARIANT_CLASS %in% fs_var,]

gene_true <- Ex_tab_Pathogenic_C4C5$clinvar_Name
gene_true <- gsub("\\):*.*$", "", gene_true)
gene_true <- gsub("^.*\\(", "", gene_true)
true_id <- ifelse(Ex_tab_Pathogenic_C4C5$gene_symbol %in% gene_true, 1, 0)
Ex_tab_Pathogenic_C4C5$gene_true <- gene_true
gene_map_rescue_path_C45 <- Ex_tab_Pathogenic_C4C5[true_id == 0,]
clinvar_Name_ind <- ifelse(grepl("\\_\\(p", gene_map_rescue_path_C45$clinvar_Name), 1, 0)

gene_map_rescue_path_C45$clinvar_Name_ind <- clinvar_Name_ind
gene_map_rescue_path_C45 <- gene_map_rescue_path_C45[gene_map_rescue_path_C45$clinvar_Name_ind == 1,]
gene_rescue_path_C4C5 <- unique(gene_map_rescue_path_C45[,c(9,120)])

rescue_path_C4C5 <- gene_map_rescue_path_C45[gene_map_rescue_path_C45$vep_consequence %nin% "intron_variant",]
rescue_path_C4C5$vep_consequence <- ifelse(grepl("\\=", rescue_path_C4C5$clinvar_Name), "synonymous_variant", rescue_path_C4C5$vep_consequence)
rescue_path_C4C5 <- rescue_path_C4C5[rescue_path_C4C5$vep_consequence %nin% "synonymous_variant",]
rescue_path_C4C5$vep_consequence <- ifelse(grepl("fs", rescue_path_C4C5$clinvar_Name), "frameshift_variant", rescue_path_C4C5$vep_consequence)
rescue_path_C4C5$vep_consequence <- ifelse(grepl("Ter\\)", rescue_path_C4C5$clinvar_Name), "stop_gained", rescue_path_C4C5$vep_consequence)
rescue_path_C4C5$vep_consequence <- ifelse(grepl("p.Met1", rescue_path_C4C5$clinvar_Name), "start_lost", rescue_path_C4C5$vep_consequence)
rescue_path_C4C5$vep_consequence <- ifelse(rescue_path_C4C5$vep_VARIANT_CLASS %in% "SNV" & 
                                             rescue_path_C4C5$vep_consequence %nin% "stop_gained" &
                                             rescue_path_C4C5$vep_consequence %nin% "start_lost", "missense_variant", rescue_path_C4C5$vep_consequence)
rescue_path_C4C5 <- rescue_path_C4C5[rescue_path_C4C5$vep_biotype %in% "protein_coding",]


##flip the gene_symbol to clinvar_Gene and vep_consequence terms for these rescued variants based on whether they are fs or otherwise.
##populate accurate protein length by searching for protein length from gene_symbol column..
prot_len <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/UCSC_tables_GRCh37_RefSeq_genes_20190716_gene_length_for_canonical_and_longest_transcripts.txt",
                       sep = "\t", header = T, stringsAsFactors = F)
canonical_gene_len <- rescue_path_C4C5$CANONICAL_GENE_LENGTH_IN_AMIMO_ACIDS
rescue_path_C4C5$CANONICAL_GENE_LENGTH_IN_AMIMO_ACIDS <- prot_len[match(rescue_path_C4C5$gene_true, prot_len$gene), 2]
rescue_path_C4C5$CANONICAL_GENE_LENGTH_IN_AMIMO_ACIDS <- ifelse(is.na(rescue_path_C4C5$CANONICAL_GENE_LENGTH_IN_AMIMO_ACIDS),
                                                                canonical_gene_len, rescue_path_C4C5$CANONICAL_GENE_LENGTH_IN_AMIMO_ACIDS)
longest_gene_len <- rescue_path_C4C5$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS
rescue_path_C4C5$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS <- prot_len[match(rescue_path_C4C5$gene_true, prot_len$gene), 3]
rescue_path_C4C5$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS <- ifelse(is.na(rescue_path_C4C5$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS),
                                                              longest_gene_len, rescue_path_C4C5$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS)
rescue_path_C4C5$filt_tab <- paste(rescue_path_C4C5$SAMPLE, rescue_path_C4C5$VARIANT, sep = "_")
rescue_path_C4C5$gene_symbol <- rescue_path_C4C5$gene_true

##Add protein position
vep_comb_pos <- gsub("^*.*p\\.","", rescue_path_C4C5$clinvar_Name)
library(tidyr)
vep_comb_pos <- extract_numeric(vep_comb_pos)
rescue_path_C4C5$vep_comb_pos <- vep_comb_pos
##QC
t1 <- as.numeric(rescue_path_C4C5$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS) - as.numeric(rescue_path_C4C5$vep_comb_pos)
table(t1 < 0)

#rescue_path_C4C5 <- rescue_path_C4C5[,-120] ## remove "clinvar_Name_ind"
##remove rescue_path_C4C5 variants from processed Ex_tab later
#Ex_tab$filt_tab <- paste(Ex_tab$SAMPLE, Ex_tab$VARIANT, sep = "_")
#Ex_tab <- Ex_tab[Ex_tab$filt_tab %nin% rescue_path_C4C5$filt_tab,] # remove mis-annotated variants

##Add reformatted/rescued variants
#Ex_tab <- rbind.data.frame(Ex_tab, rescue_path_C4C5)
#Ex_tab <- Ex_tab[,-120] #remove filt_tab

##remove QC2 failed samples: 55 MGRB, 6 RISC, 7 LIONS, 18 ISKS
QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv", header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_fail <- QC2_dat[QC2_dat$passes_qc2 %in% "FALSE",]$new_sampleid

#DT_gene_tab <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/DT_gene_tab_all.rds")
##ASPREE cancer
load("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/ASPREE_can.RData")
ASPREE_can_samp <- as.character(ASPREE_can$samp[!is.na(as.character(ASPREE_can$samp))])
##45nUP cancer
can269 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/cancer_45UP_269.csv", 
                     sep = ",", header = F, stringsAsFactors = F)

##CHIP filter(CHIP1 is Mark's algorithm, CHIP2 is tweaked version of CHIP1 to detect CHIP in ISKS)
MGRB2_chip <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/discovery_set/MGRB2_chip1.all_categories_nb2_depth_and_vaf_and_position_tallies.tsv",
                         header = T, sep = "\t", stringsAsFactors = F)
MGRB2_chip_pos <- MGRB2_chip[grepl("^CN", MGRB2_chip$algorithm_CN_category),]$sample_id

Ex_tab1$gnomad_AF_NFE <- gsub("\\[", "", Ex_tab1$gnomad_AF_NFE)
Ex_tab1$gnomad_AF_NFE <- gsub("\\]", "", Ex_tab1$gnomad_AF_NFE)
class(Ex_tab1$gnomad_AF_NFE) <- "numeric"
Ex_tab1$gnomad_AF_NFE <- ifelse(is.na(Ex_tab1$gnomad_AF_NFE), 0, Ex_tab1$gnomad_AF_NFE)

##extract genuine protein coding variants using clinvar_Name
#Ex_tab1 <- Ex_tab1[Ex_tab1$vep_consequence %nin% "intron_variant",]
clinvar_Name_ind <- ifelse(grepl("\\_\\(p", Ex_tab1$clinvar_Name), 1, 
                           ifelse(grepl("\\:p",Ex_tab1$vep_hgvsp), 2, 0))
Ex_tab1$clinvar_Name_ind <- clinvar_Name_ind
#clinvar_Name_ind == 1
#c(1:9,11,17,22:24,30,32,33,53,73,80,87,118:122)
##tier1 variants
Ex_tab2 <- Ex_tab1[Ex_tab1$clinvar_Name_ind == 1,]
Ex_tab2$vep_consequence <- ifelse(grepl("\\=", Ex_tab2$clinvar_Name), "synonymous_variant", Ex_tab2$vep_consequence)
Ex_tab2 <- Ex_tab2[Ex_tab2$vep_consequence %nin% "synonymous_variant",]
Ex_tab2$vep_consequence <- ifelse(grepl("fs", Ex_tab2$clinvar_Name), "frameshift_variant", Ex_tab2$vep_consequence)
Ex_tab2$vep_consequence <- ifelse(grepl("Ter\\)", Ex_tab2$clinvar_Name), "stop_gained", Ex_tab2$vep_consequence)
Ex_tab2$vep_consequence <- ifelse(grepl("p.Met1", Ex_tab2$clinvar_Name), "start_lost", Ex_tab2$vep_consequence)
#Ex_tab2$vep_consequence <- ifelse(Ex_tab2$vep_VARIANT_CLASS %in% "SNV" & 
#                                    Ex_tab2$vep_consequence %nin% "stop_gained" &
#                                    Ex_tab2$vep_consequence %nin% "start_lost", "missense_variant", Ex_tab2$vep_consequence)
Ex_tab2$filt_tab <- paste(Ex_tab2$SAMPLE, Ex_tab2$VARIANT, sep = "_")
Ex_tab2 <- Ex_tab2[Ex_tab2$vep_biotype %in% "protein_coding",]
##rectify protein length
canonical_gene_len <- Ex_tab2$CANONICAL_GENE_LENGTH_IN_AMIMO_ACIDS
##extract clinvar_gene
gene_true <- Ex_tab2$clinvar_Name
gene_true <- gsub("\\):*.*$", "", gene_true)
gene_true <- gsub("^.*\\(", "", gene_true)
Ex_tab2$gene_true <- gene_true

##Add protein lengths
canonical_gene_len <- Ex_tab2$CANONICAL_GENE_LENGTH_IN_AMIMO_ACIDS
Ex_tab2$CANONICAL_GENE_LENGTH_IN_AMIMO_ACIDS <- prot_len[match(Ex_tab2$gene_true, prot_len$gene), 2]
Ex_tab2$CANONICAL_GENE_LENGTH_IN_AMIMO_ACIDS <- ifelse(is.na(Ex_tab2$CANONICAL_GENE_LENGTH_IN_AMIMO_ACIDS),
                                                                canonical_gene_len, Ex_tab2$CANONICAL_GENE_LENGTH_IN_AMIMO_ACIDS)
longest_gene_len <- Ex_tab2$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS
Ex_tab2$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS <- prot_len[match(Ex_tab2$gene_true, prot_len$gene), 3]
Ex_tab2$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS <- ifelse(is.na(Ex_tab2$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS),
                                                       canonical_gene_len, Ex_tab2$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS)


vep_comb_pos <- gsub("^*.*p\\.","", Ex_tab2$clinvar_Name)
vep_comb_pos <- gsub("\\_.*$", "", vep_comb_pos)
vep_comb_pos <- gsub("fsTer.*$", "", vep_comb_pos)
library(tidyr)
vep_comb_pos <- extract_numeric(vep_comb_pos)
vep_comb_pos <- ifelse(is.na(vep_comb_pos), Ex_tab2$vep_Protein_position, vep_comb_pos)
Ex_tab2$vep_comb_pos <- vep_comb_pos
t1 <- as.numeric(Ex_tab2$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS) - as.numeric(Ex_tab2$vep_comb_pos)
Ex_tab2$vep_comb_pos <- ifelse(t1 < 0, Ex_tab2$vep_Protein_position, Ex_tab2$vep_comb_pos)
Ex_tab2$vep_comb_pos <- ifelse(is.na(Ex_tab2$vep_comb_pos), vep_comb_pos, Ex_tab2$vep_comb_pos)
#Ex_tab2_gene_mis <- unique(Ex_tab2[Ex_tab2$gene_symbol %nin% Ex_tab2$clinvar_Gene,c(9,32)])
#gene_mis2 <- unique(Ex_tab2[which(t1 < 0),c(9,32)])
#Ex_tab2_gene_mis <- unique(rbind.data.frame(Ex_tab2_gene_mis,gene_mis2))
Ex_tab2$gene_symbol <- Ex_tab2$gene_true
# remove mis-annotated variants
Ex_tab2_comb <- Ex_tab2[Ex_tab2$filt_tab %nin% rescue_path_C4C5$filt_tab,]
##add rescued variants from rescue_path_C4C5
Ex_tab2_comb <- rbind.data.frame(Ex_tab2_comb, rescue_path_C4C5)

#clinvar_Name_ind == 2; tier2 variants; rescue splice variants in coding regions; exclude intronic variants
Ex_tab3 <- Ex_tab1[Ex_tab1$clinvar_Name_ind != 1,]
Ex_tab3 <- Ex_tab3[Ex_tab3$vep_consequence %nin% "synonymous_variant",]
Ex_tab3 <- Ex_tab3[Ex_tab3$vep_biotype %in% "protein_coding",]
Ex_tab3$filt_tab <- paste(Ex_tab3$SAMPLE, Ex_tab3$VARIANT, sep = "_")
Ex_tab3 <- Ex_tab3[Ex_tab3$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS != 0,]
#vep_comb_pos <- ifelse(is.na(Ex_tab3$vep_hgvsp), 0, Ex_tab3$vep_hgvsp)
vep_comb_pos <- ifelse(is.na(Ex_tab3$vep_hgvsp), Ex_tab3$vep_Protein_position, Ex_tab3$vep_hgvsp)
vep_comb_pos <- gsub("^*.*p\\.","", vep_comb_pos)
vep_comb_pos <- gsub("^Ter","", vep_comb_pos)
vep_comb_pos <- gsub("Ter$","", vep_comb_pos)
vep_comb_pos <- gsub("Ter.*$","", vep_comb_pos)
vep_comb_pos <- gsub("\\_.*$","", vep_comb_pos)
vep_comb_pos <- gsub("\\%3D$","", vep_comb_pos)
vep_comb_pos <- extract_numeric(vep_comb_pos)
vep_comb_pos <- ifelse(is.na(vep_comb_pos), 0, vep_comb_pos)
#vep_comb_pos <- ifelse(is.na(vep_comb_pos), Ex_tab3$vep_Protein_position, vep_comb_pos)
Ex_tab3$vep_comb_pos <- vep_comb_pos
##Process Ex_tab3$vep_Protein_position
Ex_tab3$vep_Protein_position <- ifelse(is.na(Ex_tab3$vep_Protein_position), 0, Ex_tab3$vep_Protein_position)
#vep_comb_pos <- as.numeric(gsub("\\-","", vep_comb_pos))

Ex_tab3$vep_comb_pos <- ifelse(Ex_tab3$vep_comb_pos == 0 & Ex_tab3$vep_Protein_position != 0,
                                       Ex_tab3$vep_Protein_position, Ex_tab3$vep_comb_pos)

t1 <- as.numeric(Ex_tab3$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS) - as.numeric(Ex_tab3$vep_comb_pos)
Ex_tab3$vep_comb_pos <- ifelse(t1 < 0, Ex_tab3$vep_Protein_position, Ex_tab3$vep_comb_pos)
Ex_tab3$vep_comb_pos <- ifelse(is.na(Ex_tab3$vep_comb_pos), vep_comb_pos, Ex_tab3$vep_comb_pos)
##combine Ex_tab2 and Ex_tab3 and rescued variants
Ex_tab2_comb <- Ex_tab2_comb[,-122]
Ex_tab23_comb <- rbind.data.frame(Ex_tab2_comb, Ex_tab3)
dim(Ex_tab23_comb)
rm(Ex_tab, Ex_tab2_comb, Ex_tab2, Ex_tab3)
##remove SAMPLE CR57
Ex_tab_filt1 <- Ex_tab23_comb[(Ex_tab23_comb$SAMPLE %nin% "CR57"),]

#don't remove "synonymous_variant" 
toMatch <- c("3_prime_UTR_variant", "5_prime_UTR_variant", "intron_variant", 
             "upstream_gene_variant", "non_coding_transcript_exon_variant",
             "non_coding_transcript_variant", "downstream_gene_variant", "NMD_transcript_variant",
             "intergenic_variant", "mature_miRNA_variant")
#"NMD_transcript_variant"
#Ex_tab_filt2 <- Ex_tab_filt1[!(Ex_tab_filt1$vep_consequence %in% toMatch),]
Ex_tab_filt2 <- Ex_tab_filt1[!(grepl(paste(toMatch,collapse="|"), 
                                     Ex_tab_filt1$vep_consequence)),]

rm(Ex_tab1, Ex_tab2, Ex_tab3, Ex_tab_filt1)
#Ex_tab_filt2$is_MGRB <- ifelse(Ex_tab_filt2$SAMPLE %in% MGRB_all, 1, 0)
Ex_tab_filt2$is_CH <- ifelse(Ex_tab_filt2$SAMPLE %in% MGRB2_chip_pos, 1, 0)
Ex_tab_filt2$is_ASPC <- ifelse(Ex_tab_filt2$SAMPLE %in% ASPREE_can_samp, 1, 0)
Ex_tab_filt2$is_45up <- ifelse(Ex_tab_filt2$SAMPLE %in% can269$V1, 1, 0)

##filter out MGRB samples that have CH
Ex_tab_noCH_filt3 <- Ex_tab_filt2[Ex_tab_filt2$is_ASPC == 0 & Ex_tab_filt2$is_45up == 0,]

rm(Ex_tab_filt2)
##Capture variant present in the last 5% of the length of the protein sequence

dim(Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$vep_consequence %in% "frameshift_variant",])

##remove RNA genes
Ex_tab_noCH_filt3 <- Ex_tab_noCH_filt3[which(!is.na(as.numeric(Ex_tab_noCH_filt3$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS))),]

Ex_tab_noCH_filt3 <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS != 0 
                                       & Ex_tab_noCH_filt3$CANONICAL_GENE_LENGTH_IN_AMIMO_ACIDS != 0,]
Ex_tab_noCH_filt3$last_five_perc <- ifelse((as.numeric(Ex_tab_noCH_filt3$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS) - as.numeric(Ex_tab_noCH_filt3$vep_comb_pos))/
                                             (as.numeric(Ex_tab_noCH_filt3$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS)) <= 0.05, 1,
                                           0)

##David algorithm JUl-22-2019
##Total Exome count
toMatch_var <- c("synonymous_variant", "frameshift_variant", "start_lost", "stop_gained",
                 "missense_variant", "splice_region_variant",
                 "splice_donor_variant", "splice_acceptor_variant",
                 "protein_altering_variant")
Ex_tab_noCH_filt3$vep_con1 <- unlist(lapply(strsplit(Ex_tab_noCH_filt3$vep_consequence, split = "&"), function(x)x[1]))
Ex_tab_noCH_filt3 <- Ex_tab_noCH_filt3[grepl(paste(toMatch_var,collapse="|"), 
                                             Ex_tab_noCH_filt3$vep_con1),]
##Changed on Aug12 after DT meeting
##for sept05_rect dataset
##DT nee suggestions: Oct-2-2019
##If C4 == ClinVar Benign; assign eigenphred score
# toMatch_ben <- c("Benign", "Benign/Likely_benign")
# Ex_tab_noCH_filt3$comb_score <- as.numeric(ifelse(grepl(paste(toMatch_ben,collapse="|"), 
#                                                         Ex_tab_noCH_filt3$clinvar_Clinical_Significance) & 
#                                                     Ex_tab_noCH_filt3$auto_call %in% "C4", Ex_tab_noCH_filt3$EigenPhred,
#                                                   Ex_tab_noCH_filt3$comb_score))
# #frameshift donot have eigenphred scores, hence they are assigned a score of 15 manually
# Ex_tab_noCH_filt3$comb_score <- ifelse(is.na(Ex_tab_noCH_filt3$comb_score), 15, Ex_tab_noCH_filt3$comb_score)
# ##last 5 percent of protein sequence based score penalty
# Cterm_C3_C4 <- c("C3", "C4")
# Ex_tab_noCH_filt3$comb_score <- as.numeric(ifelse(Ex_tab_noCH_filt3$auto_call %in% Cterm_C3_C4 &
#                                                     Ex_tab_noCH_filt3$last_five_perc == 1, (Ex_tab_noCH_filt3$comb_score)/2,
#                                                   Ex_tab_noCH_filt3$comb_score))
# 
# ##remove benign variants
# 
# Ex_tab_noCH_filt3 <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$comb_score > 0,]

##remove samples with multiple frameshift mutations in the same gene
library(tidyr)
library(dplyr)
print("remove GATK artefacts")
## frameshifts; these are handled better in comb_ISKS_MGRB_var_filt.R script using same-gene-same-sample
#Ex_tab_noCH_filt3 %>% filter(gene_symbol %in% "BRIP1" & vep_consequence %in% "frameshift_variant") %>% group_by(SAMPLE, gene_symbol) %>% summarise(n = n()) 
fs_all <- Ex_tab_noCH_filt3[grepl("frameshift_variant", Ex_tab_noCH_filt3$vep_consequence),]
t1 <- fs_all %>% group_by(SAMPLE, gene_symbol) %>% summarise(n = n()) 

t2 <- as.data.frame(t1)

t3 <- t2[t2$n > 1,] #remove these variants by subsetting the frameshift variants from these genes in the corresponding samples before saving
#t3[order(t3$n, decreasing = T),]

Ex_tab_noCH_filt3_fs <- list()
for(j in 1:dim(t3)[1]){
  f1 <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$SAMPLE %in% t3$SAMPLE[j] &
                            Ex_tab_noCH_filt3$gene_symbol %in% t3$gene_symbol[j],]
  f1 <- f1[grep("frameshift_variant",f1$vep_consequence), ]
 # f1 <- f1[order(f1$VAF, f1$comb_score, decreasing = T),]
  f1 <- f1[order(f1$VAF, decreasing = T),]
  f1 <- f1[-1,]
#  print(j)
  # Ex_tab_noCH_filt3_fs[[j]] <- f1[grep("frameshift_variant",f1$vep_consequence), ]
 # if(dim(f1)[1] > 0){
    Ex_tab_noCH_filt3_fs[[j]] <- f1
 # }
#  else { Ex_tab_noCH_filt3_fs[[j]] <- NULL}
}
rm_fs <- do.call("rbind.data.frame", Ex_tab_noCH_filt3_fs)
#c(1:9,11,17,22:24,30,32,33,53,73,80,87,118:122)
## mis-sense
t1_mis <- Ex_tab_noCH_filt3 %>% filter(vep_consequence %in% "missense_variant") %>% group_by(SAMPLE, gene_symbol) %>% summarise(n = n()) 

t2_mis <- as.data.frame(t1_mis)

t3_mis <- t2_mis[t2_mis$n > 3,] #cull all variants found in samples with more than 3 mis sense variants per genes for recessive genes

Ex_tab_noCH_filt3_ms <- list()
for(j in 1:dim(t3_mis)[1]){
  f1 <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$SAMPLE %in% t3_mis$SAMPLE[j] &
                            Ex_tab_noCH_filt3$gene_symbol %in% t3_mis$gene_symbol[j],]
  Ex_tab_noCH_filt3_ms[[j]] <- f1[f1$vep_consequence %in% "missense_variant", ]
  
}

rm_ms <- do.call("rbind.data.frame", Ex_tab_noCH_filt3_ms)
#rm_ms$filt_tab <- paste(rm_ms$SAMPLE, rm_ms$VARIANT, sep = "_")

##for synonymous variants
## mis-sense
t1_syn <- Ex_tab_noCH_filt3[grepl("synonymous_variant",Ex_tab_noCH_filt3$vep_consequence),] %>% group_by(SAMPLE, gene_symbol) %>% summarise(n = n()) 

t2_syn <- as.data.frame(t1_syn)

t3_syn <- t2_syn[t2_syn$n > 3,]

Ex_tab_noCH_filt3_syn <- list()
if(dim(t3_syn)[1] > 0){
for(j in 1:dim(t3_syn)[1]){
  f3 <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$SAMPLE %in% t3_syn$SAMPLE[j] &
                            Ex_tab_noCH_filt3$gene_symbol %in% t3_syn$gene_symbol[j],]
  
  Ex_tab_noCH_filt3_syn[[j]] <- f3[f3$vep_consequence %in% "synonymous_variant", ]
  
}
} else{Ex_tab_noCH_filt3_syn <- NULL }

if(!is.null(Ex_tab_noCH_filt3_syn)){
rm_syn <- do.call("rbind.data.frame", Ex_tab_noCH_filt3_syn)
}else{rm_syn <- NULL}

##combine the sets to remove from 

rm_fin <- rbind.data.frame(rm_fs, rm_ms, rm_syn)
rm_fin$filt_tab <- paste(rm_fin$SAMPLE, rm_fin$VARIANT, sep = "_")
##retain C5 variants, "start_lost", "stop_gained"
#rm_fin <- rm_fin[rm_fin$auto_call %nin% "C5",]
#rm_fin <- rm_fin[rm_fin$vep_consequence %nin% c("start_lost", "stop_gained"),]

Ex_tab_noCH_filt3$filt_tab <- paste(Ex_tab_noCH_filt3$SAMPLE, Ex_tab_noCH_filt3$VARIANT, sep = "_")


Ex_tab_noCH_filt3_fin <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$filt_tab %nin% rm_fin$filt_tab,]

#Tag QC2 failed
##new step added from repset variant filtering in addition to changes in MGRB chip set
Ex_tab_noCH_filt3_fin$QC2 <- ifelse(Ex_tab_noCH_filt3_fin$SAMPLE %in% QC2_dat_fail, "fail", "pass")

#remove benigns and QC2 failed samples
##QC2 failed samples: ("1383", "AABTU", "1", "AACDY")
Ex_tab_noCH_filt3_fin_filt <- Ex_tab_noCH_filt3_fin[Ex_tab_noCH_filt3_fin$QC2 %in% "pass",]

#Ex_tab_noCH_filt3_fin_filt <- Ex_tab_noCH_filt3_fin[as.character(Ex_tab_noCH_filt3_fin$auto_call) %nin% "B",]

Ex_tab_noCH_filt3_fin_filt$set <- c("MGRB_AR_AD")


write.table(Ex_tab_noCH_filt3_fin_filt, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/all_mgrb_combset2020_variants_AR_AD_exome_count_clingene_rnd3_freeze.tsv",
            row.names = F, quote = F, sep = "\t")


