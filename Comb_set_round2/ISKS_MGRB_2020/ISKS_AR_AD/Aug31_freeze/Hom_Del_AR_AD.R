.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)
library(dplyr)
library(tidyr)
print("Reading input files from MGRB and ISKS")
Ex_tab_mgrb <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/MGRB/all_mgrb_combset2020_variants_AR_AD_all_fields_clingene_rnd3.tsv", sep = "\t", header = T, stringsAsFactors = F)
Ex_tab_isks <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/all_isks_combset2020_variants_AR_AD_all_fields_clingene_rnd3_freeze.tsv", sep = "\t", header = T, stringsAsFactors = F)
Ex_tab_mgrb$set <- gsub("ISKS", "MGRB", Ex_tab_mgrb$set)

Ex_tab <- rbind.data.frame(Ex_tab_mgrb, Ex_tab_isks)
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

##remove duplicates
print("removing duplicates")
fil_tab_noCH <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
Ex_tab <- Ex_tab[Ex_tab$SAMPLE %in% fil_tab_noCH$SAMPLE, ]
Ex_tab <- Ex_tab[!is.na(Ex_tab$SAMPLE),]

##Extract all deletions in the data involving PID genes

#PIDgenes
mirabello_genes <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/mirabello_genes.txt",
                              sep = "", header = F, stringsAsFactors = F)
##remove trailing and leading white spaces
mira_genes <- mirabello_genes %>% 
  mutate(V1 = str_trim(mirabello_genes$V1, side = "both"))

##mitotic check point genes
mito_chk_point = read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/mitotic_checkpoint/mitocheck_point.txt",
                            sep = "", header = F, skip = 1, stringsAsFactors = F)

##centrosome maturation complex
centrosome_mat <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/Centrosome maturation.xlsx",
                             sheet = 1)
centrosome_mat <- as.data.frame(centrosome_mat)
centrosome_mat <- centrosome_mat[centrosome_mat$MoleculeType %nin% "Chemical Compounds",]

##cep_haus C3C4C5
cep_haus_all <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CEP_HAUS_C345.txt",
                           sep = "", header = F, stringsAsFactors = F)

top_SKAT_cgc <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/DT_skatbin_cgc.txt",
                           sep = "", header = F, stringsAsFactors = F)
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

##Add mito_chk_point
print("Add mito_chk_point genes")
gene_sub <- unique(c(top_SKAT_cgc$V1, cgc_genes$Gene.Symbol, sheltrin_complex, Telo_extension, Sheltrin_comp_extn, 
                     mito_chk_point$V1, mira_genes$V1, centrosome_mat$X__1, cep_haus_all$V1))

Ex_tab_PID <- Ex_tab[Ex_tab$gene_symbol %in% gene_sub, ]
##deletions
Ex_tab_PID_del <- Ex_tab_PID[Ex_tab_PID$vep_VARIANT_CLASS %in% c("deletion", "indel"),]

##remove FS within 5 bp of other FS; also check for compound hets
all_genes <- unique(Ex_tab_PID_del$gene_symbol)
#extract hits for removal; all variants within 5 bps of the fs variant
sam_sel_ns_rm <- list()
sam_sel_ns_rm_ls <- list()
for(i in 1:length(all_genes)){
  gene_df <- Ex_tab_PID_del[Ex_tab_PID_del$gene_symbol %in% all_genes[i],]
  sam_freq <- as.data.frame(table(gene_df$SAMPLE))
  sam_sel <- as.character(sam_freq[sam_freq$Freq >1,1])
  if(length(sam_sel) > 0){
    for(k in 1:length(sam_sel)){
      sam_sel_df <- gene_df[gene_df$SAMPLE %in% sam_sel[k],]
      #  pos_fs <- sam_sel_df[sam_sel_df$vep_consequence %in% "frameshift_variant",]$pos
      pos_fs <- sam_sel_df[grepl("frameshift_variant", sam_sel_df$vep_consequence),]$pos
      #   pos_nfs <- sam_sel_df[sam_sel_df$vep_consequence %nin% "frameshift_variant",]$pos
      pos_nfs <- sam_sel_df[sam_sel_df$pos != pos_fs,]$pos
      sel_var <- ifelse(abs(pos_fs - pos_nfs) <= 5, 1, 0) 
      pos_nfs <- pos_nfs[which(sel_var == 1)]
      sam_sel_ns_rm[[k]] <-  sam_sel_df[sam_sel_df$pos %in% pos_nfs,]
    }
    sam_sel_ns_rm_ls[[i]] <- do.call("rbind.data.frame",sam_sel_ns_rm)
  }
  else{sam_sel_ns_rm_ls[[i]] <- NULL}
}
sam_rm_het_df <- do.call("rbind.data.frame",sam_sel_ns_rm_ls)
sam_rm_het_df <- unique(sam_rm_het_df)

##Ex_tab compound het filter: removes all variants within 5 bps of the FS variant
Ex_tab_PID_del <- Ex_tab_PID_del[Ex_tab_PID_del$filt_tab %nin% sam_rm_het_df$filt_tab,]
##compound het filtering
s1 <- Ex_tab_PID_del %>% group_by(SAMPLE, gene_symbol) %>% summarise(n = n()) 
s2 <- as.data.frame(s1)
s3 <- s2[s2$n > 1,]
PID_cpd_het_pos_fin <- Ex_tab_PID_del[Ex_tab_PID_del$SAMPLE %in% s3$SAMPLE,]
##No compound het deletions; as expected.


##Monoallelic deletions
Ex_tab_PID_mon_del <- Ex_tab_PID_del[Ex_tab_PID_del$vep_VARIANT_CLASS %in% "deletion" & Ex_tab_PID_del$GT == 1,]
##Homozygous deletions
Ex_tab_PID_hom_del <- Ex_tab_PID_del[Ex_tab_PID_del$vep_VARIANT_CLASS %in% "deletion" & Ex_tab_PID_del$GT == 2,]
Ex_tab_PID_hom_del_isks <- Ex_tab_PID_hom_del[Ex_tab_PID_hom_del$set %in% "ISKS_AR_AD",]
View(Ex_tab_PID_hom_del_isks[,c(1:4,9,11,23:24,87)])
