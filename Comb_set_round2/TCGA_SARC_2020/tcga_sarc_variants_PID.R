.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)


Ex_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/TCGA_SARC_AR_AD/all_tcga_sarc_combset2020_variants_AR_AD_all_fields_clingene_rnd3.tsv", sep = "\t", header = T, 
                     stringsAsFactors = F)
Ex_tab$gnomad_AF <- gsub("\\[", "", Ex_tab$gnomad_AF)
Ex_tab$gnomad_AF <- gsub("\\]", "", Ex_tab$gnomad_AF)
class(Ex_tab$gnomad_AF) <- "numeric"
Ex_tab$gnomad_AF_NFE <- gsub("\\[", "", Ex_tab$gnomad_AF_NFE)
Ex_tab$gnomad_AF_NFE <- gsub("\\]", "", Ex_tab$gnomad_AF_NFE)
class(Ex_tab$gnomad_AF_NFE) <- "numeric"
Ex_tab$gnomad_AF <- ifelse(is.na(Ex_tab$gnomad_AF), 0, Ex_tab$gnomad_AF)
Ex_tab$swegen_AF <- ifelse(is.na(Ex_tab$swegen_AF), 0, Ex_tab$swegen_AF)
Ex_tab$gnomad_AF_NFE <- ifelse(is.na(Ex_tab$gnomad_AF_NFE), 0, Ex_tab$gnomad_AF_NFE)
##GQ filter

Ex_tab <- Ex_tab[Ex_tab$GQ >= 80,]

##AD
##C4/C5; gnomad_AF_NFE <=  0.0002
#Ex_tab_C4C5 <- Ex_tab[Ex_tab$auto_call %in% c("C4", "C5") & Ex_tab$gnomad_AF < 0.01,]
Ex_tab_C4C5 <- Ex_tab[Ex_tab$auto_call %in% c("C4", "C5") & Ex_tab$gnomad_AF_NFE <= 0.0002,]
##Add intra_cohort MAF filter cohort_MAF < 0.001 ##Added later cohort_MAF <= 3/2*253
Ex_tab_C4C5 <- Ex_tab_C4C5[Ex_tab_C4C5$cohort_MAF <= 0.006 & Ex_tab_C4C5$swegen_AF <= 0.001,]
##C3; gnomad_AF == 0
#Ex_tab_C3 <- Ex_tab[Ex_tab$auto_call %in% "C3" & Ex_tab$gnomad_AF <= 0.01 & Ex_tab$gnomad_AF_NFE == 0,]
Ex_tab_C3 <- Ex_tab[Ex_tab$auto_call %in% "C3" & Ex_tab$gnomad_AF_NFE  == 0 & Ex_tab$swegen_AF == 0,]
##cohort_MAF <= 0.001
Ex_tab_C3 <- Ex_tab_C3[Ex_tab_C3$cohort_MAF <= 0.006 & Ex_tab_C3$comb_score >= 5.6,]



##For SKAT(use the same set of variants used above, compound het filtered, gnomad_AF_NFE filtered)
Ex_tab_C3C4C5_AD_filt <- unique(rbind.data.frame(Ex_tab_C4C5, Ex_tab_C3))
rm_X_GT2_AD <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$GT == 2 & grepl("^X", Ex_tab_C3C4C5_AD_filt$VARIANT), ]
Ex_tab_C3C4C5_AD_filt <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$filt_tab %nin% rm_X_GT2_AD$filt_tab,]

##overlap with PID genes
#PIDgenes
#PIDgenes
top_SKAT_cgc <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/DT_skatbin_cgc.txt",
                           sep = "", header = F, stringsAsFactors = F)
cgc_genes <- read.delim("~/RVAS/cancer_gene_census_hg37.csv", sep = ",", header = T, stringsAsFactors = F)
#remove fusion genes from cgc list
#cgc_genes$Gene.Symbol[(grep("fusion", cgc_genes$Role.in.Cancer))]
sheltrin_complex <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "ACD")
Telo_extension <- c("TIMELESS", "TIPIN", "FANCM", "BRCA1", "BLM") ##replication stress is a source of telomere recombination
Sheltrin_comp_extn = c("ACD", "POT1", "TERF1", "TERF2", "TERF2IP", "TINF2", "ATM", 
                       "BAG3", "BLM", "BRCA1", "CALD1", "CLK3", "DCLRE1B", "FANCD2", 
                       "FBXO4", "HSPA4", "KIAA1191", "MRE11A", "NBN", "PINX1", "PRKDC", 
                       "RAD50", "SLX4", "STUB1", "TNKS", "TNKS2", "U2AF2", "UCHL1", 
                       "WRN", "XRCC5", "XRCC6")
gene_sub <- unique(c(top_SKAT_cgc$V1, cgc_genes$Gene.Symbol, sheltrin_complex, Telo_extension, Sheltrin_comp_extn))

table(Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$gene_symbol %in% gene_sub,]$auto_call)

##NF1, WRN, MTHFR were found

#Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$gene_symbol %in% c("POT1", "TINF2", "TERF1","TP53", "NF1", "TERF2IP", "ACD", "NBN", "CCT2", "CCT6A", "TCP1"),c(1:9, 127:128)]
##clueGO modules
mods <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/clueGO_mod_gene.txt",
                   sep = "", header = F, stringsAsFactors = F)
mod_genes_filt <- Ex_tab_C3C4C5_AD_filt[Ex_tab_C3C4C5_AD_filt$gene_symbol %in% mods$V1, c(1:3,9,11,127:128)]
write.table(mod_genes_filt, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/TCGA_SARC_AR_AD/mod_genes_filt.tsv",
            sep = "\t", row.names = F, quote = F)

Ex_tab_C3 <- Ex_tab[Ex_tab$auto_call %in% "C3" & Ex_tab$gnomad_AF_NFE  == 0 & Ex_tab$swegen_AF == 0,]
##cohort_MAF <= 0.001
Ex_tab_C3 <- Ex_tab_C3[Ex_tab_C3$cohort_MAF <= 0.006,]

##For SKAT(use the same set of variants used above, compound het filtered, gnomad_AF_NFE filtered)
Ex_tab_C4C5_AD_C3unfilt <- unique(rbind.data.frame(Ex_tab_C4C5, Ex_tab_C3))
rm_X_GT2_AD_unfilt <- Ex_tab_C4C5_AD_C3unfilt[Ex_tab_C4C5_AD_C3unfilt$GT == 2 & grepl("^X", Ex_tab_C4C5_AD_C3unfilt$VARIANT), ]
Ex_tab_C4C5_AD_C3unfilt <- Ex_tab_C4C5_AD_C3unfilt[Ex_tab_C4C5_AD_C3unfilt$filt_tab %nin% rm_X_GT2_AD_unfilt$filt_tab,]

mod_genes_unfilt <- Ex_tab_C4C5_AD_C3unfilt[Ex_tab_C4C5_AD_C3unfilt$gene_symbol %in% mods$V1, c(1:3,9,11,127:128)]
write.table(mod_genes_unfilt, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/TCGA_SARC_AR_AD/mod_genes_C3_unfilt.tsv",
            sep = "\t", row.names = F, quote = F)

##make table

gene_tab <- function(gene_sym, var_file){
  var_file1 <- var_file[var_file$gene_symbol %in% gene_sym,]
  sarc <- dim(var_file1)[1]
  case_wt <- sum(var_file1$comb_score)
  case_call_auto <- paste(apply(as.data.frame(table(as.character(var_file1$auto_call))) , 1 , paste , collapse = ":" ), collapse = ",")
  case_vep_var <- paste(apply(as.data.frame(table(as.character(var_file1$vep_consequence))) , 1 , paste , collapse = ":" ), collapse = ",")
  res <- cbind.data.frame("TCGA_SARC" = sarc, "gene" = gene_sym, 
                          "case_wt" = case_wt,  "case_call_auto" = case_call_auto, 
                          "case_vep_var" = case_vep_var)
  return(res)
}

#gene_tab("TP53", mod_genes_filt)

library(doParallel)
library(doMC)
registerDoMC(30)

res20 <- list()


#for(k in 1:length(unique(fil_tab_noCH$gene_symbol))){
res_list <- list()
#  genes <-  as.character(Exome_pc123_srt_SKAT[[k]][,1])
genes <- unique(mods$V1)
system.time(res_list <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
{gene_tab(genes[i], mod_genes_filt)})
res_filt <- do.call("rbind.data.frame", res_list)
write.table(res_filt, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/TCGA_SARC_AR_AD/filt_clueGO_mod_tab.tsv",
            sep = "\t", row.names = F, quote = F)

res_list1 <- list()

system.time(res_list1 <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
{gene_tab(genes[i], mod_genes_unfilt)})
res_unfilt <- do.call("rbind.data.frame", res_list1)
write.table(res_unfilt, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/TCGA_SARC_AR_AD/unfilt_clueGO_mod_tab.tsv",
            sep = "\t", row.names = F, quote = F)
