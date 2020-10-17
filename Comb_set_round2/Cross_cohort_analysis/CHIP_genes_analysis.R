##Statistics of CHIP related genes in ISKS, MGRB, EPIT and CAIRNS datasets
##filtering criteria: only remove spurious variants(more than two variants per gene per sample)
##specific gene list: c("TP53", "JAK2", "KRAS", "NRAS", "TET2", "IDH1", "DNMT3A", "KIT",  "NPM1", "IDH2", 
##"RUNX1", "SF3B1", "MYD88", "NOTCH1", "CALR",  "ASXL1", "FLT3", "SRSF2", "WT1", "TERT")
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
library(ggplot2)
library(dplyr)

`%nin%` = Negate(`%in%`)

# isks_tsv <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/all_shards_ISKS_combset2020.tsv", sep = "\t",
#                        header = T, stringsAsFactors = F)
# isks_chip_tsv <- isks_tsv[isks_tsv$gene_symbol %in% c("TP53", "JAK2", "KRAS", "NRAS", "TET2", "IDH1", "DNMT3A", "KIT",  "NPM1", "IDH2", 
#                                                       "RUNX1", "SF3B1", "MYD88", "NOTCH1", "CALR",  "ASXL1", "FLT3", "SRSF2", "WT1", "TERT"),]
# isks_chip_tsv$set <- c("ISKS")
# 
# mgrb_tsv <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/MGRB/all_shards_MGRB_combset2020.tsv", sep = "\t",
#                        header = T, stringsAsFactors = F)
# mgrb_chip_tsv <- mgrb_tsv[mgrb_tsv$gene_symbol %in% c("TP53", "JAK2", "KRAS", "NRAS", "TET2", "IDH1", "DNMT3A", "KIT",  "NPM1", "IDH2", 
#                                                       "RUNX1", "SF3B1", "MYD88", "NOTCH1", "CALR",  "ASXL1", "FLT3", "SRSF2", "WT1", "TERT"),]
# mgrb_chip_tsv$set <- c("MGRB")
# 
# epit_tsv <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/all_shards_EPIT_combset2020.tsv", sep = "\t",
#                                    header = T, stringsAsFactors = F)
# epit_chip_tsv <- epit_tsv[epit_tsv$gene_symbol %in% c("TP53", "JAK2", "KRAS", "NRAS", "TET2", "IDH1", "DNMT3A", "KIT",  "NPM1", "IDH2", 
#                                                       "RUNX1", "SF3B1", "MYD88", "NOTCH1", "CALR",  "ASXL1", "FLT3", "SRSF2", "WT1", "TERT"),]
# epit_chip_tsv$set <- c("EPIT")
# 
# cairns_tsv <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CAIRNS/all_shards_CAIRNS.tsv", sep = "\t",
#                          header = T, stringsAsFactors = F)
# cairns_chip_tsv <- cairns_tsv[cairns_tsv$gene_symbol %in% c("TP53", "JAK2", "KRAS", "NRAS", "TET2", "IDH1", "DNMT3A", "KIT",  "NPM1", "IDH2", 
#                                                       "RUNX1", "SF3B1", "MYD88", "NOTCH1", "CALR",  "ASXL1", "FLT3", "SRSF2", "WT1", "TERT"),]
# 
# cairns_chip_tsv$set <- c("CAIRNS")
# 
# 
# ##filter MGRB samples from analysis
# ##remove QC2 failed samples: 55 MGRB, 6 RISC, 7 LIONS, 18 ISKS
# QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv", header = T, sep = "\t", stringsAsFactors = F)
# QC2_dat_fail <- QC2_dat[QC2_dat$passes_qc2 %in% "FALSE",]$new_sampleid
# ##ASPREE cancer
# load("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/ASPREE_can.RData")
# ASPREE_can_samp <- as.character(ASPREE_can$samp[!is.na(as.character(ASPREE_can$samp))])
# ##45nUP cancer
# can269 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/cancer_45UP_269.csv", 
#                      sep = ",", header = F, stringsAsFactors = F)
# 
# ##CHIP filter(CHIP1 is Mark's algorithm, CHIP2 is tweaked version of CHIP1 to detect CHIP in ISKS)
# MGRB2_chip <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/discovery_set/MGRB2_chip1.all_categories_nb2_depth_and_vaf_and_position_tallies.tsv",
#                          header = T, sep = "\t", stringsAsFactors = F)
# MGRB2_chip_pos <- MGRB2_chip[grepl("^CN", MGRB2_chip$algorithm_CN_category),]$sample_id
# 
# mgrb_chip_tsv <- mgrb_chip_tsv[mgrb_chip_tsv$SAMPLE %nin% QC2_dat_fail, ]
# mgrb_chip_tsv <- mgrb_chip_tsv[mgrb_chip_tsv$SAMPLE %nin% ASPREE_can_samp, ]
# mgrb_chip_tsv <- mgrb_chip_tsv[mgrb_chip_tsv$SAMPLE %nin% can269$V1, ]
# mgrb_chip_tsv <- mgrb_chip_tsv[mgrb_chip_tsv$SAMPLE %nin% MGRB2_chip_pos, ]
# 
# comb_chip_all <- rbind.data.frame(isks_chip_tsv, mgrb_chip_tsv, epit_chip_tsv, cairns_chip_tsv)
# write.table(comb_chip_all, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/Cross_cohort_analysis/chip_genes_cohort.tsv",
#             sep = "\t", quote = F, row.names = F)

##start here
##Add C3/C4/C5 and EigenScores
comb_chip_all <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/Cross_cohort_analysis/chip_genes_cohort.tsv",
                            sep = "\t", header = T, stringsAsFactors = F)
comb_chip_all <- comb_chip_all[as.numeric(comb_chip_all$VAF) >= 0.05 & as.numeric(comb_chip_all$DP) >= 10,]
Ex_samp_id <- unique(comb_chip_all$SAMPLE)
#remove duplicates
dup_samp <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/fin_samp_dup_drop.txt", sep = "", header = T,
                       stringsAsFactors = F)
`%nin%` = Negate(`%in%`)
Ex_samp_id <- Ex_samp_id[Ex_samp_id %nin% dup_samp$x]
comb_chip_all <- comb_chip_all[comb_chip_all$SAMPLE %in% Ex_samp_id,]

comb_chip_all$gnomad_AF_NFE <- gsub("\\[", "", comb_chip_all$gnomad_AF_NFE)
comb_chip_all$gnomad_AF_NFE <- gsub("\\]", "", comb_chip_all$gnomad_AF_NFE)
class(comb_chip_all$gnomad_AF_NFE) <- "numeric"
comb_chip_all$gnomad_AF_NFE <- ifelse(is.na(comb_chip_all$gnomad_AF_NFE), 0, comb_chip_all$gnomad_AF_NFE)
#TERT variants
#comb_chip_all_TERT <- comb_chip_all[comb_chip_all$gene_symbol %in% "TERT",]

##filter out non-coding variants; splice variants were not removed
toMatch <- c("3_prime_UTR_variant", "5_prime_UTR_variant", "intron_variant", 
             "synonymous_variant", "upstream_gene_variant", "non_coding_transcript_exon_variant",
             "non_coding_transcript_variant", "downstream_gene_variant", "NMD_transcript_variant",
             "intergenic_variant", "mature_miRNA_variant")
#"NMD_transcript_variant"
Ex_tab_noCH_filt3 <- comb_chip_all[!(grepl(paste(toMatch,collapse="|"), 
                                      comb_chip_all$vep_consequence)),]

##Capture variant present in the last 5% of the length of the protein sequence

Ex_tab_noCH_filt3$vep_Protein_position <- ifelse(is.na(Ex_tab_noCH_filt3$vep_Protein_position), 0, Ex_tab_noCH_filt3$vep_Protein_position)
class(Ex_tab_noCH_filt3$vep_Protein_position) <- "numeric"
Ex_tab_noCH_filt3$vep_Protein_position <- ifelse(is.na(Ex_tab_noCH_filt3$vep_Protein_position), 0, Ex_tab_noCH_filt3$vep_Protein_position)
vep_hgvsp <- ifelse(is.na(Ex_tab_noCH_filt3$vep_hgvsp), 0, Ex_tab_noCH_filt3$vep_hgvsp)
vep_hgvsp <- gsub("^*.*p\\.","", vep_hgvsp)
vep_hgvsp <- gsub("Ter*.*$","", vep_hgvsp)
library(tidyr)
vep_hgvsp <- extract_numeric(vep_hgvsp)
#combine position from vep_Protein_position and vep_hgvsp
vep_comb_pos <- ifelse(vep_hgvsp == 0, Ex_tab_noCH_filt3$vep_Protein_position, vep_hgvsp)
Ex_tab_noCH_filt3$vep_comb_pos <- vep_comb_pos

dim(Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$vep_consequence %in% "frameshift_variant",])

vep_hgvsp <- ifelse(is.na(Ex_tab_noCH_filt3$vep_hgvsp), 0, Ex_tab_noCH_filt3$vep_hgvsp)
vep_hgvsp <- gsub("^*.*p\\.","", vep_hgvsp)
vep_hgvsp <- gsub("Ter*.*$","", vep_hgvsp)
library(tidyr)
vep_hgvsp <- extract_numeric(vep_hgvsp)
#combine position from vep_Protein_position and vep_hgvsp
vep_comb_pos <- ifelse(vep_hgvsp == 0, Ex_tab_noCH_filt3$vep_Protein_position, vep_hgvsp)
Ex_tab_noCH_filt3$vep_comb_pos <- vep_comb_pos

dim(Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$vep_consequence %in% "frameshift_variant",])

Ex_tab_noCH_filt3$vep_comb_pos <- ifelse(is.na(Ex_tab_noCH_filt3$vep_comb_pos), 0, Ex_tab_noCH_filt3$vep_comb_pos)
##remove RNA genes
Ex_tab_noCH_filt3 <- Ex_tab_noCH_filt3[which(!is.na(as.numeric(Ex_tab_noCH_filt3$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS))),]

Ex_tab_noCH_filt3 <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS != 0 
                                       & Ex_tab_noCH_filt3$CANONICAL_GENE_LENGTH_IN_AMIMO_ACIDS != 0,]
Ex_tab_noCH_filt3$last_five_perc <- ifelse((as.numeric(Ex_tab_noCH_filt3$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS) - Ex_tab_noCH_filt3$vep_comb_pos)/
                                             (as.numeric(Ex_tab_noCH_filt3$LONGEST_GENE_LENGTH_IN_AMINO_ACIDS)) <= 0.05, 1,
                                           0)

##David algorithm JUl-22-2019
toMatch1 <- c("Likely_pathogenic", "Pathogenic", "Pathogenic/Likely_pathogenic", 
              "Pathogenic/Likely_pathogenic/drug_response", "Pathogenic/Likely_pathogenic/other",
              "Pathogenic/Affects", "Pathogenic/risk_factor", "Pathogenic/Likely_pathogenic/risk_factor",
              "risk_factor")


##modify scoring function since the annotation from the new VEP are slightly different
toMatch2 <- c("frameshift_variant", "start_lost", "stop_gained")

sco_fun5 <- function(filt3_inp){
  ##construct C5,C4, C3 and score them as sco_fun4  
  ##All clinvar_sig variants
  filt3_inp$vep_con1 <- unlist(lapply(strsplit(filt3_inp$vep_consequence, split = "&"), function(x)x[1]))
  metric1 <-  ifelse(grepl(paste(toMatch1,collapse="|"), 
                           filt3_inp$clinvar_Clinical_Significance) , "C5", 
                     ifelse(grepl(paste(toMatch2,collapse="|"), 
                                  filt3_inp$vep_con1), "C4", 
                            ifelse(filt3_inp$vep_con1 %in% 
                                     c("missense_variant", "splice_region_variant",
                                       "splice_donor_variant", "splice_acceptor_variant",
                                       "protein_altering_variant")
                                   & filt3_inp$gnomad_AF_NFE == 0, "C3", "B")))
  ##lapply(strsplit(var_type, split = "&"), function(x)x[1]) think about this
  
  metric2 <- as.numeric(ifelse(metric1 %in% "C5", 45, 
                               ifelse(metric1 %in% "C4", 45, 
                                      ifelse(metric1 %in% "C3", filt3_inp$EigenPhred, 0))))
  
  metric2[is.na(metric2)] <- 0
  
  ##changed on Aug-11-2019: Adjust this outside the function later
  #metric3 <- as.numeric(ifelse(metric1 %in% "C4" & filt3_inp$last_exon == 1, -10, 0))
  #metric3 <- as.numeric(ifelse(metric1 %in% "C4" & filt3_inp$last_exon == 1, -25, 0))
  
  #met <- metric2 + metric3
  
  #met_df <- cbind.data.frame("auto_call" = metric1, "comb_score" = met )
  met_df <- cbind.data.frame("auto_call" = metric1, "comb_score" = metric2 )
  return(met_df)
}


Ex_tab_noCH_filt3 <- cbind.data.frame(Ex_tab_noCH_filt3, sco_fun5(Ex_tab_noCH_filt3))

##Changed on Aug12 after DT meeting
##for sept05_rect dataset
##DT nee suggestions: Oct-2-2019
##If C4 == ClinVar Benign; assign eigenphred score
toMatch_ben <- c("Benign", "Benign/Likely_benign")
Ex_tab_noCH_filt3$comb_score <- as.numeric(ifelse(grepl(paste(toMatch_ben,collapse="|"), 
                                                        Ex_tab_noCH_filt3$clinvar_Clinical_Significance) & 
                                                    Ex_tab_noCH_filt3$auto_call %in% "C4", Ex_tab_noCH_filt3$EigenPhred,
                                                  Ex_tab_noCH_filt3$comb_score))
#frameshift donot have eigenphred scores, hence they are assigned a score of 15 manually
Ex_tab_noCH_filt3$comb_score <- ifelse(is.na(Ex_tab_noCH_filt3$comb_score), 15, Ex_tab_noCH_filt3$comb_score)
##last 5 percent of protein sequence based score penalty
Cterm_C3_C4 <- c("C3", "C4")
Ex_tab_noCH_filt3$comb_score <- as.numeric(ifelse(Ex_tab_noCH_filt3$auto_call %in% Cterm_C3_C4 &
                                                    Ex_tab_noCH_filt3$last_five_perc == 1, (Ex_tab_noCH_filt3$comb_score)/2,
                                                  Ex_tab_noCH_filt3$comb_score))

##remove benign variants

#Ex_tab_noCH_filt3 <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$auto_call %nin% "B",]
Ex_tab_noCH_filt3 <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$comb_score > 0,]
##filter : remove spurious variants(more than two variants per gene per sample)
## mis-sense
t1_mis <- Ex_tab_noCH_filt3 %>% filter(vep_consequence %in% "missense_variant") %>% group_by(SAMPLE, gene_symbol) %>% summarise(n = n()) 

t2_mis <- as.data.frame(t1_mis)

t3_mis <- t2_mis[t2_mis$n > 2,] #cull sample variants with more than 2 mis sense variants per genes

Ex_tab_noCH_filt3_ms <- list()
for(j in 1:dim(t3_mis)[1]){
  f1 <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$SAMPLE %in% t3_mis$SAMPLE[j] &
                            Ex_tab_noCH_filt3$gene_symbol %in% t3_mis$gene_symbol[j],]
  Ex_tab_noCH_filt3_ms[[j]] <- f1[f1$vep_consequence %in% "missense_variant", ]
  
}

rm_ms <- do.call("rbind.data.frame", Ex_tab_noCH_filt3_ms)

## frameshifts
t0_fs <- Ex_tab_noCH_filt3[grep("frameshift_variant", Ex_tab_noCH_filt3$vep_consequence),]

#t1 <- Ex_tab_noCH_filt3 %>% filter(vep_consequence %in% "frameshift_variant") %>% group_by(SAMPLE, gene_symbol) %>% summarise(n = n()) 
t1 <- t0_fs %>% group_by(SAMPLE, gene_symbol) %>% summarise(n = n()) 
  
t2 <- as.data.frame(t1)

t3 <- t2[t2$n > 1,] #remove these variants by subsetting the frameshift variants from these genes in the corresponding samples before saving
#t3[order(t3$n, decreasing = T),]

Ex_tab_noCH_filt3_fs <- list()
for(j in 1:dim(t3)[1]){
  f1 <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$SAMPLE %in% t3$SAMPLE[j] &
                            Ex_tab_noCH_filt3$gene_symbol %in% t3$gene_symbol[j],]
  Ex_tab_noCH_filt3_fs[[j]] <- f1[grep("frameshift_variant",f1$vep_consequence), ]
  
}
rm_fs <- do.call("rbind.data.frame", Ex_tab_noCH_filt3_fs)

##combine the sets to remove  

rm_fin <- rbind.data.frame(rm_fs, rm_ms)
rm_fin$filt_tab <- paste(rm_fin$SAMPLE, rm_fin$VARIANT, sep = "_")
rm_fin <- rm_fin[rm_fin$auto_call %nin% "C5",]
Ex_tab_noCH_filt3$filt_tab <- paste(Ex_tab_noCH_filt3$SAMPLE, Ex_tab_noCH_filt3$VARIANT, sep = "_")

Ex_tab_noCH_filt3_fin_filt <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$filt_tab %nin% rm_fin$filt_tab,]

##remove variants from samples with more than 2 variants per gene
t_msfs <- Ex_tab_noCH_filt3_fin_filt %>% group_by(SAMPLE, gene_symbol) %>% summarise(n = n())
tdf_msfs <- as.data.frame(t_msfs[t_msfs$n > 3, ])
Ex_tab_noCH_filt3_msfs <- list()
if(dim(tdf_msfs)[1] > 0){
for(j in 1:dim(tdf_msfs)[1]){
  Ex_tab_noCH_filt3_msfs[[j]] <- Ex_tab_noCH_filt3_fin_filt[Ex_tab_noCH_filt3_fin_filt$SAMPLE %in% tdf_msfs$SAMPLE[j] &
                                                              Ex_tab_noCH_filt3_fin_filt$gene_symbol %in% tdf_msfs$gene_symbol[j],]
  # Ex_tab_noCH_filt3_msfs[[j]] <- f1[f1$vep_consequence %in% "missense_variant", ]
}
rm_ms_fs <- do.call("rbind.data.frame", Ex_tab_noCH_filt3_msfs)
rm_ms_fs$filt_tab <- paste(rm_ms_fs$SAMPLE, rm_ms_fs$VARIANT, sep = "_")
rm_ms_fs <- rm_ms_fs[rm_ms_fs$auto_call %nin% "C5",]
Ex_tab_noCH_filt3_fin_filt <- Ex_tab_noCH_filt3_fin_filt[Ex_tab_noCH_filt3_fin_filt$filt_tab %nin% rm_ms_fs$filt_tab,]
 } else {
Ex_tab_noCH_filt3_fin_filt <- Ex_tab_noCH_filt3_fin_filt
 }

write.table(Ex_tab_noCH_filt3_fin_filt, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/Cross_cohort_analysis/CHIP_variants_ISKS_EPIT_CAIRNS_MGRB.tsv",
            sep = "\t", quote = F, row.names = F)

##Include promoter variants for TERT 
Ex_tab_trunc <- Ex_tab_noCH_filt3_fin_filt[,c(1:121)] 
upstream_match <- c("5_prime_UTR_variant", "upstream_gene_variant")
TERT_var_upstream <- comb_chip_all[comb_chip_all$gene_symbol %in% "TERT" & grepl(paste(upstream_match,collapse="|"), 
                                                                                  comb_chip_all$vep_consequence),]
Ex_tab_inc_TERT_prom <- rbind.data.frame(Ex_tab_trunc, TERT_var_upstream)
write.table(Ex_tab_inc_TERT_prom, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/Cross_cohort_analysis/CHIP_variants_w_TERT_prom_ISKS_EPIT_CAIRNS_MGRB.tsv",
            sep = "\t", quote = F, row.names = F)

library(ggplot2)
#t1_df <- Ex_tab_noCH_filt3_fin_filt %>% filter(vep_consequence %in% "missense_variant") %>% 
t1_df <- Ex_tab_noCH_filt3_fin_filt %>% filter(auto_call %nin% "C3") %>% 
  group_by(set, gene_symbol) %>% summarise(n = n()) 

t1_df <- Ex_tab_noCH_filt3_fin_filt %>% filter(auto_call %nin% "C3") %>% 
  group_by(set, gene_symbol) %>% summarise(Eigen_scores = sum(comb_score))
#t1_df <- Ex_tab_noCH_filt3_fin_filt %>% group_by(set, gene_symbol) %>% summarise(Eigen_scores = sum(comb_score))
t2_df <- as.data.frame(t1_df)
#ggplot(t2_df, aes(set, gene_symbol, fill= n)) + geom_tile() 
ggplot(t2_df, aes(gene_symbol, Eigen_scores)) + geom_point(aes(colour = factor(set))) + ggtitle("Cumulative EigenPhred Scores across Cohorts")
ggplot(t2_df, aes(gene_symbol, Eigen_scores)) + geom_point(aes(colour = factor(set))) + ggtitle("Cumulative EigenPhred Scores across Cohorts minus C3")
# ##make C3/C4/C5
# ##remove intronic variants
# comb_chip_all_fin_ns <- comb_chip_all_fin[comb_chip_all_fin$vep_consequence %nin% "intron_variant", ]
# comb_chip_all_fin_ns <- comb_chip_all_fin_ns[comb_chip_all_fin_ns$vep_consequence %nin% "synonymous_variant", ]
# comb_chip_all_fin_ns <- comb_chip_all_fin_ns[comb_chip_all_fin_ns$vep_consequence %nin% "3_prime_UTR_variant", ]
# comb_chip_all_fin_ns <- comb_chip_all_fin_ns[comb_chip_all_fin_ns$vep_consequence %nin% "non_coding_transcript_exon_variant", ]
# comb_chip_all_fin_ns <- comb_chip_all_fin_ns[comb_chip_all_fin_ns$vep_consequence %nin% "downstream_gene_variant", ]
# comb_chip_all_fin_TERT <- comb_chip_all_fin_ns[comb_chip_all_fin_ns$gene_symbol %in% "TERT",]
# comb_chip_all_fin_ns <- comb_chip_all_fin_ns[comb_chip_all_fin_ns$vep_consequence %nin% "5_prime_UTR_variant", ]
# comb_chip_all_fin_ns <- comb_chip_all_fin_ns[comb_chip_all_fin_ns$vep_consequence %nin% "upstream_gene_variant", ]
# comb_chip_all_fin_ns_inc <- rbind.data.frame(comb_chip_all_fin_ns, comb_chip_all_fin_TERT)
# 
# ##gene_cohort distribution
# 
# t1 <- comb_chip_all_fin_ns_inc %>% filter(vep_consequence %in% "missense_variant") %>% group_by(set, gene_symbol) %>% summarise(n = n()) 
# 
# t2 <- as.data.frame(t1)

#library(ggplot2)
#ggplot(t2, aes(set, gene_symbol, fill= n)) + geom_tile() + theme_ipsum()
