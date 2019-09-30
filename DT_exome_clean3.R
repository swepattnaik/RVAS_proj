##Uses scoring function without splice acceptor and donor as C4
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

##Aug6 remove VAF NaN's(multiallelic split issue; error in code for HAIL computation)
##Aug7 remove NMD variants
##Rcommendation for SKAT; use a hard filter of comb_score > 5 compared to comb_score >= 5
##Aug12 changed scoring function, increased penalty for last exon variants (-25)
##same as DT_call_analysis.R (refer it for source file details)
##Add filter using VAF > 0.3
##Use scores for C5 = C4 > C3 = Eigen
##weights assigned: C5 -> 45; C4 -> 45; C3 -> Eigen; last 5% protein sequence: score/2
##Scoring function: C5 := Clinvar Pathogenic; C4 := c("frameshift_variant", "start_lost", "stop_gained") 



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

##filter out MGRB samples that have CH

Ex_tab_noCH_filt3 <- Ex_tab_filt2[!(Ex_tab_filt2$is_CH == 0 & Ex_tab_filt2$is_MGRB == 1),]
Ex_tab_noCH_filt3 <- Ex_tab_noCH_filt3[Ex_tab_noCH_filt3$is_ASPC == 0,]


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
toMatch2 <- c("frameshift_variant", "start_lost", "stop_gained")

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
saveRDS(Ex_tab_noCH_filt3, file = "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/Exome_filt3_nCH_C5eqC4_nonmds_iskrisc_05Sept_rect_ASP.rds", compress = T)


#########QC
#var_check <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/var_control_set_Aug12.tsv", sep = "\t",
#           header = T)

#table(var_check$VARIANT %in% Ex_tab_noCH_filt3$VARIANT)