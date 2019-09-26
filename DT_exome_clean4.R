##DT genes
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

##Aug6 remove VAF NaN's(multiallelic split issue; error in code for HAIL computation)
##Aug7 remove NMD variants
##Rcommendation for SKAT; use a hard filter of comb_score > 5 compared to comb_score >= 5
##Aug12 changed scoring function, increased penalty for last exon variants (-25)

##test shard chromosome 1 for testing new feature
#tshard <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test_shard.tsv", header = T, sep = "\t",
 #                    stringsAsFactors = F)


##same as DT_call_analysis.R (refer it for source file details)
##Add filter using VAF > 0.3
##Use scores for C5[30], C4[30], C3 [10-20] DT_gene_tab_ISKS_MGRB_nCH_fil1$DavidThomas_call
##weights assigned: C5 -> 30; C4 -> 30; C3 -> 15; stop_gained -> 30; clin_var(pathogenic) <- 35; indels <- 10
#Ex_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/comb_isks_risc_exome_latest2.tsv", sep = "\t",
#           header = T, stringsAsFactors = F)

#Ex_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/comb_isks_latest1Aug.tsv", sep = "\t",
#                                header = T, stringsAsFactors = F)

#Ex_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/comb_isks_latest4Aug.tsv", sep = "\t",
#                     header = T, stringsAsFactors = F)

#Ex_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/comb_isks_risc_latest12Aug.tsv", sep = "\t",
#                     header = T, stringsAsFactors = F)

Ex_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/comb_isks_risc_latest05Sept_rect1.tsv", sep = "\t",
                                          header = T, stringsAsFactors = F)
Ex_tab <- Ex_tab[!is.na(Ex_tab$SAMPLE),]
Ex_tab$SAMPLE <- gsub(".variant.*$", "", Ex_tab$SAMPLE) ##check this in the next run
#remove variants ending with Asterisk; these donot have vep annotation
Ex_tab <- Ex_tab[!is.na(Ex_tab$gene_symbol),]

Ex_tab <- Ex_tab[!is.na(Ex_tab$DP),] ##There are records with NA in DP attribute ; file: sept05_rect

#Ex_tab <- Ex_tab[is.na(Ex_tab$VAF),] ##VAF's have NaNs as DP = 0. Its a division by zero for VAF

#Ex_tab <- Ex_tab[!is.na(Ex_tab$VAF),]
#Ex_tab <- Ex_tab[Ex_tab$DP != 0,]
#table(Ex_tab[Ex_tab$VAF == "NaN",]$gene_symbol %in% c("TP53", "TINF2", "POT1", "SMARCAL1"))
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
#Ex_tab_filt1 <- Ex_tab_filt1[!(Ex_tab_filt1$SAMPLE %in% "CR57"),]
Ex_tab_filt1 <- Ex_tab1[!(Ex_tab1$SAMPLE %in% "CR57"),]

toMatch <- c("3_prime_UTR_variant", "5_prime_UTR_variant", "intron_variant", 
             "synonymous_variant", "upstream_gene_variant", "non_coding_transcript_exon_variant",
             "non_coding_transcript_variant", "downstream_gene_variant", "NMD_transcript_variant",
             "intergenic_variant", "mature_miRNA_variant")
#"NMD_transcript_variant"
#Ex_tab_filt2 <- Ex_tab_filt1[!(Ex_tab_filt1$vep_consequence %in% toMatch),]
Ex_tab_filt2 <- Ex_tab_filt1[!(grepl(paste(toMatch,collapse="|"), 
                                    Ex_tab_filt1$vep_consequence)),]


Ex_tab_filt2$is_MGRB <- ifelse(Ex_tab_filt2$SAMPLE %in% MGRB_all, 1, 0)

Ex_tab_filt2$is_ASPC <- ifelse(Ex_tab_filt2$SAMPLE %in% ASPREE_can$samp, 1, 0)
#length(unique(Ex_tab_filt2[Ex_tab_filt2$is_CH == 1 & Ex_tab_filt2$is_MGRB == 1,]$SAMPLE))


#saveRDS(Ex_tab_filt2, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/comb_isks_risc_exome_latest_filt2.rds", 
 #       compress = T)
#saveRDS(Ex_tab_filt2, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/comb_isks_latest1Aug_filt2.rds", 
#        compress = T)


##filter out MGRB samples that have CH
#Ex_tab_filt2 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/comb_isks_latest30Jul_filt2.rds")
#Ex_tab_filt2 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/comb_isks_latest1Aug_filt2.rds")
Ex_tab_noCH_filt3 <- Ex_tab_filt2[!(Ex_tab_filt2$is_CH == 0 & Ex_tab_filt2$is_MGRB == 1),]
Ex_tab_noCH_filt3 <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$is_ASPC == 0,]
#dim(Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$gene_symbol %in% "TP53",])

##filter somatic snps
#DT_gene_tab_ISKS_MGRB_nCH_fil3 <- DT_gene_tab_ISKS_MGRB_nCH_fil3[is.na(DT_gene_tab_ISKS_MGRB_nCH_fil3$`Somatic?`),]

##last 2 exon check; this input contains vep_exon filed while specifies the exon number relative to the last exon number.
#last2 <- unlist(lapply(strsplit(TP53_ex$vep_exon, split = "/"), function(x)ifelse(as.numeric(x[2]) - as.numeric(x[1]) <= 1, 1, 0)))
#last2[is.na(last2)] <- 0

#Ex_tab_noCH_filt3$last_exon <- unlist(lapply(strsplit(Ex_tab_noCH_filt3$vep_exon, split = "/"), function(x)ifelse(as.numeric(x[2]) - as.numeric(x[1]) <= 1, 1, 0)))
##capture only last exon: based on DT meeting Aug-12-2019

# Ex_tab_noCH_filt3$last_exon <- unlist(lapply(strsplit(Ex_tab_noCH_filt3$vep_exon, split = "/"), function(x)
#   ifelse(as.numeric(x[2]) <= 5, 0,
#          ifelse(as.numeric(x[2]) - as.numeric(x[1]) < 1, 1, 0)))) ##only last exon
# Ex_tab_noCH_filt3$last_exon[is.na(Ex_tab_noCH_filt3$last_exon)] <- 0


##Capture variant present in the last 5% of the length of the protein sequence

Ex_tab_noCH_filt3$vep_Protein_position <- ifelse(is.na(Ex_tab_noCH_filt3$vep_Protein_position), 0, Ex_tab_noCH_filt3$vep_Protein_position)
class(Ex_tab_noCH_filt3$vep_Protein_position) <- "numeric"
Ex_tab_noCH_filt3$vep_Protein_position <- ifelse(is.na(Ex_tab_noCH_filt3$vep_Protein_position), 0, Ex_tab_noCH_filt3$vep_Protein_position)

#vep_hgvsp <- as.numeric(unlist(regmatches(vep_hgvsp, gregexpr("[[:digit:]]+", vep_hgvsp))))
#vep_hgvsp <- as.numeric(gsub("[^\\d]", "", vep_hgvsp, perl=TRUE))
#vep_hgvsp1 <- regmatches(vep_hgvsp, gregexpr("[[:digit:]]+", vep_hgvsp))
#vep_hgvsp2 <- as.numeric(unlist(lapply(vep_hgvsp1, function(x)x[1])))
#vep_comb_pos <- ifelse(vep_hgvsp2 == 0, Ex_tab_noCH_filt3$vep_Protein_position, vep_hgvsp2) #works but takes long

vep_hgvsp <- ifelse(is.na(Ex_tab_noCH_filt3$vep_hgvsp), 0, Ex_tab_noCH_filt3$vep_hgvsp)
vep_hgvsp <- gsub("^*.*p\\.","", vep_hgvsp)
vep_hgvsp <- gsub("Ter*.*$","", vep_hgvsp)
library(tidyr)
vep_hgvsp <- extract_numeric(vep_hgvsp)
#combine position from vep_Protein_position and vep_hgvsp
vep_comb_pos <- ifelse(vep_hgvsp == 0, Ex_tab_noCH_filt3$vep_Protein_position, vep_hgvsp)
Ex_tab_noCH_filt3$vep_comb_pos <- vep_comb_pos

dim(Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$vep_consequence %in% "frameshift_variant",])
#check redundant records
#dim(Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$vep_consequence %in% "frameshift_variant" & Ex_tab_noCH_filt3$vep_comb_pos == 0,])
#Ex_tab_noCH_filt3 <- Ex_tab_noCH_filt3[!is.na(Ex_tab_noCH_filt3$SAMPLE),]
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
              "Pathogenic/Affects", "Pathogenic/risk_factor", "Pathogenic/Likely_pathogenic/risk_factor")


##modify scoring function since the annotation from the new VEP are slightly different
toMatch2 <- c("frameshift_variant", "start_lost", "stop_gained",
              "splice_donor_variant", "splice_acceptor_variant")

##sco_fun5 changed on Aug-11-2019 to further underweight last exon variants (aftermath of POT1 variants in MGRB)
##sco_fun5 changed on Sept-05-2019 to further underweight variants that are C4 and present in the last
## 5% of the C-terminal tail of the protein.(aftermath of POT1 variants in MGRB)
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



##Experiments

##scoring function5
#DT_gene_tab_ISKS_MGRB_nCH_fil3 <- filt3_inp
Ex_tab_noCH_filt3 <- cbind.data.frame(Ex_tab_noCH_filt3, sco_fun5(Ex_tab_noCH_filt3))
Cterm_C3_C4 <- c("C3", "C4")
##Changed on Aug12 after DT meeting
##for sept05_rect dataset
Ex_tab_noCH_filt3$comb_score <- as.numeric(ifelse(Ex_tab_noCH_filt3$auto_call %in% Cterm_C3_C4 &
                                                    Ex_tab_noCH_filt3$last_five_perc == 1, (Ex_tab_noCH_filt3$comb_score)/2,
                                                  Ex_tab_noCH_filt3$comb_score))

#saveRDS(Ex_tab_noCH_filt3, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_noNAs_auto.rds", compress = T)
#saveRDS(Ex_tab_noCH_filt3, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_noNAs_auto_4Aug.rds", compress = T)
#saveRDS(Ex_tab_noCH_filt3, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_noNAs_auto_6Aug.rds", compress = T)
#saveRDS(Ex_tab_noCH_filt3, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_auto_7Aug.rds", compress = T)
#saveRDS(Ex_tab_noCH_filt3, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_12Aug.rds", compress = T)
#saveRDS(Ex_tab_noCH_filt3, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_05Sept_rect.rds", compress = T)
#saveRDS(Ex_tab_noCH_filt3, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_05Sept_splice.rds", compress = T)
saveRDS(Ex_tab_noCH_filt3, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_05Sept_splice_ASP.rds", compress = T)


######################SKAT analysis###################

# ##Additive model
# Ex_samp_id <- unique(Ex_tab_noCH_filt3$SAMPLE)
# samp_vec <- list()
# for(k in 1:dim(Ex_tab_noCH_filt3)[1]){
#   sam_gene_gt <- Ex_tab_noCH_filt3[grepl(Ex_tab_noCH_filt3$VARIANT[k],
#                                          Ex_tab_noCH_filt3$VARIANT),][,c(1,3,9)]
#   #sam_gene_gt$add_mod <- ifelse(sam_gene_gt$GT %in% "0/1", 1, ifelse(sam_gene_gt$GT %in% "0/0", 0, 2))
#   sam_gene_gt$add_mod <- as.numeric(sam_gene_gt$GT)
#   sam10 <- ifelse(Ex_samp_id %in% sam_gene_gt$SAMPLE, 1, 0)
#   sam10[which(sam10 != 0)] <- sam_gene_gt$add_mod ##additive model
#   samp_vec[[k]] <- sam10
#   # samp_vec[[k]] <- ifelse(DT_samp_id %in% sam_gene$SAMPLE, 1, 0)
#   # gene_vec[[k]] <- ifelse(genes %in% sam_gene$gene_symbol, 1, 0)
# }
# 
# samp_vec_mat <- do.call("rbind", samp_vec)
# colnames(samp_vec_mat) <- Ex_samp_id
# ##Account for last_exon column
# #samp_vec_mat_df <- cbind.data.frame(DT_gene_tab_ISKS_MGRB_nCH_fil3[,c(1:2,6,14:15,178)], samp_vec_mat)
# #for auto call
# samp_vec_mat_df <- cbind.data.frame(Ex_tab_noCH_filt3[,c(2:3,9,11,121)], samp_vec_mat)
# samp_vec_mat_df <- samp_vec_mat_df[!is.na(samp_vec_mat_df$SAMPLE),]
# dim(samp_vec_mat_df[is.na(samp_vec_mat_df$comb_score),]) ##QC step
# 
# samp_vec_mat_df[is.na(samp_vec_mat_df$comb_score),] <- 0
# 
# saveRDS(samp_vec_mat_df, file = "~/RVAS/shard_sub_tier3/DT_sheet/SKAT_inp_filt_var_noCH_C5eqC4_noNA_Exome_gt.rds", compress = T)
# 
# ##Split Ex_tab_noCH_filt3 into smaller chunks to parallelise generation of samp_vec_mat_df
# ##split into 50000 rows per data frame
# Ex_tab_noCH_filt3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_noNAs_auto.rds")
# print_table <- function(df, df_name){
#   saveRDS(df, file = paste("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Genomat/", 
#                            df_name,".rds", sep = ""), compress = T)
# }
# 
# ##tried to parallelise (abandoned for now); preliminary files generated successfully
# df_sp <- split(Ex_tab_noCH_filt3, (seq(nrow(Ex_tab_noCH_filt3))-1) %/% 50000)
# saveRDS(df_sp, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/list_Ex_tab_noCH_filt3.rds", compress = T)
# #lapply(split(Ex_tab_noCH_filt3, (seq(nrow(Ex_tab_noCH_filt3))-1) %/% 50000), function(x)print_table(x))
# mapply(print_table, df_sp, paste0(names(df_sp), "sp")) ##works
# 

#########QC
Ex_tab_noCH_filt3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_05Sept_splice_ASP.rds")
var_check <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/var_control_set_Aug12.tsv", sep = "\t",
           header = T)

table(var_check$VARIANT %in% Ex_tab_noCH_filt3$VARIANT)

##centrosome genes
aug12_var_cent <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/aug_12_vars.tsv", sep = "\t",
                             header = T)
table(aug12_var_cent$VARIANT %in% Ex_tab_noCH_filt3$VARIANT)

#c("RANGAP1", "CENPC", "CENPF", "NDEL1", "MIS18A")

# DT_gene_tab_ISKS_MGRB_nCH <- DT_gene_tab_ISKS_MGRB_nCH[DT_gene_tab_ISKS_MGRB_nCH$gnomad_AF_NFE <= 0.0001,]
# DT_gene_tab_ISKS_MGRB_nCH <- DT_gene_tab_ISKS_MGRB_nCH[DT_gene_tab_ISKS_MGRB_nCH$VAF >= 0.25 & DT_gene_tab_ISKS_MGRB_nCH$VAF < 1,]
# 
# table(unique(DT_gene_tab_ISKS_MGRB_nCH$gene_symbol) %in% unique(Ex_tab_noCH_filt3$gene_symbol))
# table(unique(DT_gene_tab_ISKS_MGRB_nCH$VARIANT) %in% unique(Ex_tab_noCH_filt3$VARIANT))
# 
# 
# TP53_ex <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$gene_symbol %in% "TP53",]
# TP53_ex$last_exon <- unlist(lapply(strsplit(TP53_ex$vep_exon, split = "/"), function(x)ifelse(as.numeric(x[2]) - as.numeric(x[1]) <= 1, 1, 0)))
# TP53_DT <- DT_gene_tab_ISKS_MGRB_nCH[DT_gene_tab_ISKS_MGRB_nCH$gene_symbol %in% "TP53",]
# table(unique(TP53_DT$VARIANT) %in% unique(TP53_ex$VARIANT))

