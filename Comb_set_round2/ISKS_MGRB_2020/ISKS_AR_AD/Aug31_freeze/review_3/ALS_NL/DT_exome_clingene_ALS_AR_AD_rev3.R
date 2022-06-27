#.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)


##HAIL pipeline output with gnomadAF < 0.05
Ex_tab <- read.delim("~/RVAS/Jan_NL_Harwig/ALS_dat/comb_genes_var.tsv", sep = "\t",
                     header = T, stringsAsFactors = F)

Ex_tab <- Ex_tab[!is.na(Ex_tab$sample),]
#Ex_tab$SAMPLE <- gsub(".variant.*$", "", Ex_tab$SAMPLE) ##check this in the next run
#remove variants ending with Asterisk; these donot have vep annotation
#Ex_tab <- Ex_tab[!is.na(Ex_tab$SYMBOL),]

Ex_tab <- Ex_tab[!is.na(Ex_tab$DP),] ##There are records with NA in DP attribute ; file: sept05_rect


Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")
#Shelterin <- c("POT1", "TINF2", "TERF1", "SMARCAL1", "STAG3")
##new extended CEP_HAUS_core
#CEP_HAUS_core_extn = read.delim("~/RVAS/GIST_32/CEP_HAUS_C45_genes.txt", sep = "", header = T, stringsAsFactors = F)
#CEP_HAUS_core_extn2 = CEP_HAUS_core_extn$x[grep("CEP|HAUS", CEP_HAUS_core_extn$x)]
CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
sarcoma_genes <- c("TP53", "SDHA", "SDHB", "SDHD")
Breast_cancer_genes <- c("BRCA1", "BRCA2", "PALB2")
genes = c(Shelterin, CEP_HAUS_core, sarcoma_genes, Breast_cancer_genes)

Ex_tab <- Ex_tab[Ex_tab$SYMBOL %in% genes,]
Ex_tab <- Ex_tab[!is.na(Ex_tab$SYMBOL),]
saveRDS(Ex_tab, "~/RVAS/Jan_NL_Harwig/ALS_dat/comb_genes_var.rds", compress = T)
toMatch <- c("3_prime_UTR_variant", "5_prime_UTR_variant", "intron_variant", 
             "synonymous_variant", "upstream_gene_variant", "non_coding_transcript_exon_variant",
             "non_coding_transcript_variant", "downstream_gene_variant", "NMD_transcript_variant",
             "intergenic_variant", "mature_miRNA_variant")
#"NMD_transcript_variant"
#Ex_tab_filt2 <- Ex_tab_filt1[!(Ex_tab_filt1$vep_consequence %in% toMatch),]
Ex_tab1 <- Ex_tab[!(grepl(paste(toMatch,collapse="|"), 
                                     Ex_tab$Consequence)),]


dim(Ex_tab1[Ex_tab1$Consequence %in% "frameshift_variant",])

##David algorithm JUl-22-2019
toMatch1 <- c("^likely_pathogenic", "^pathogenic", "^pathogenic/likely_pathogenic", 
              "Pathogenic/Likely_pathogenic/drug_response", "Pathogenic/Likely_pathogenic/other",
              "Pathogenic/Affects", "Pathogenic/risk_factor", "Pathogenic/Likely_pathogenic/risk_factor",
              "^risk_factor")


##modify scoring function since the annotation from the new VEP are slightly different
toMatch2 <- c("frameshift_variant", "start_lost", "stop_gained")
toMatch2 <- c("frameshift_variant", "start_lost", "stop_gained", "splice_donor_variant", "splice_acceptor_variant")
##sco_fun5 changed on Aug-11-2019 to further underweight last exon variants (aftermath of POT1 variants in MGRB)
##sco_fun5 changed on Sept-05-2019 to further underweight variants that are C4 and present in the last
## 5% of the C-terminal tail of the protein.(aftermath of POT1 variants in MGRB)
sco_fun5 <- function(filt3_inp){
  ##construct C5,C4, C3 and score them as sco_fun4  
  ##All clinvar_sig variants
  filt3_inp$vep_con1 <- unlist(lapply(strsplit(filt3_inp$Consequence, split = ","), function(x)x[1]))
  metric1 <-  ifelse(grepl(paste(toMatch1,collapse="|"), 
                           filt3_inp$CLIN_SIG) , "C5", 
                     ifelse(grepl(paste(toMatch2,collapse="|"), 
                                  filt3_inp$vep_con1), "C4", 
                            ifelse(filt3_inp$vep_con1 %in% 
                                     c("missense_variant", "splice_region_variant",
                                       "splice_donor_variant", "splice_acceptor_variant",
                                       "protein_altering_variant")
                                   & filt3_inp$gnomAD_NFE_AF == 0, "C3", "B")))
  ##lapply(strsplit(var_type, split = "&"), function(x)x[1]) think about this
  
  # metric2 <- as.numeric(ifelse(metric1 %in% "C5", 45, 
  #                              ifelse(metric1 %in% "C4", 45, 
  #                                     ifelse(metric1 %in% "C3", filt3_inp$EigenPhred, 0))))
  # 
  # metric2[is.na(metric2)] <- 0
  
  ##changed on Aug-11-2019: Adjust this outside the function later
  #metric3 <- as.numeric(ifelse(metric1 %in% "C4" & filt3_inp$last_exon == 1, -10, 0))
  #metric3 <- as.numeric(ifelse(metric1 %in% "C4" & filt3_inp$last_exon == 1, -25, 0))
  
  #met <- metric2 + metric3
  
  #met_df <- cbind.data.frame("auto_call" = metric1, "comb_score" = met )
  #met_df <- cbind.data.frame("auto_call" = metric1, "comb_score" = metric2 )
  #return(met_df)
  return(metric1)
}



##Experiments
Ex_tab2 = Ex_tab1
##scoring function5
#DT_gene_tab_ISKS_MGRB_nCH_fil3 <- filt3_inp


Ex_tab2$auto_call <- sco_fun5(Ex_tab2)
Ex_tab3_C45 = Ex_tab2[Ex_tab2$auto_call %in% c("C4", "C5"),]


##not filtered for intra_cohort frequency
write.table(Ex_tab2, "~/RVAS/Jan_NL_Harwig/ALS_dat/ALS_Hartwig_C345_cpx_genes_var.tsv", sep = "\t",
            quote = F, row.names = F)
write.table(Ex_tab3_C45, "~/RVAS/Jan_NL_Harwig/ALS_dat/ALS_Hartwig_C45_cpx_genes_var.tsv", sep = "\t",
            quote = F, row.names = F)

Ex_tab3_C45 = read.delim("~/RVAS/Jan_NL_Harwig/ALS_dat/ALS_Hartwig_C45_cpx_genes_var.tsv", sep = "\t",
                         header = T, stringsAsFactors = F)
##remove samples with multiple frameshift mutations in the same gene
library(tidyr)
library(dplyr)

## frameshifts; these are handled better in comb_ISKS_MGRB_var_filt.R script using same-gene-same-sample
#Ex_tab_noCH_filt3 %>% filter(gene_symbol %in% "BRIP1" & vep_consequence %in% "frameshift_variant") %>% group_by(SAMPLE, gene_symbol) %>% summarise(n = n()) 

fs_all <- Ex_tab3_C45[grepl("frameshift_variant", Ex_tab3_C45$Consequence),]
t1 <- fs_all %>% group_by(sample, ID, SYMBOL) %>% summarise(n = n()) 

t2 <- as.data.frame(t1)

t3 <- t2[t2$n > 1,] #remove these variants by subsetting the frameshift variants from these genes in the corresponding samples before saving
#t3[order(t3$n, decreasing = T),]
##no frameshift artefacts found in ALS data

##remove high frequency variants based on intra-cohort frequency or any variant that is present 
#in greater than 3 samples

var_coh_maf = Ex_tab3_C45 %>% group_by(ID, SYMBOL) %>% summarise(n = n())

var_coh_maf = as.data.frame(var_coh_maf)
var_coh_maf = var_coh_maf[order(var_coh_maf$n, decreasing = T),]
var_coh_maf_filt = var_coh_maf[var_coh_maf$n < 4,]
var_coh_maf_filt$subtract_n = ifelse(var_coh_maf_filt$ID %in% t3$ID, 1, 0)
var_coh_maf_filt$final_n = var_coh_maf_filt$n - var_coh_maf_filt$subtract_n


var_coh_maf_filt[var_coh_maf_filt$SYMBOL %in% Shelterin,]
var_coh_maf_filt[var_coh_maf_filt$SYMBOL %in% CEP_HAUS_core,]
var_coh_maf_filt[var_coh_maf_filt$SYMBOL %in% sarcoma_genes,]
var_coh_maf_filt[var_coh_maf_filt$SYMBOL %in% Breast_cancer_genes,]
##fisher's test:replication set note Hartwig (n = 278); ALS (n = 2892)
##Centrosome
inp = c(3, 275, 4, 2888)
sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
colnames(sim_mat) <- c("case", "cont")
rownames(sim_mat) <- c("hits", "no_hits")
ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.95)

##Shelterin
inp = c(6, 272, 12, 2880)
sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
colnames(sim_mat) <- c("case", "cont")
rownames(sim_mat) <- c("hits", "no_hits")
ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.95)

##Sarcoma
inp = c(4, 274, 4, 2888)
sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
colnames(sim_mat) <- c("case", "cont")
rownames(sim_mat) <- c("hits", "no_hits")
ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.95)

##Check with Breast cancer genes

##Check for combined TCGA + Hartwig

