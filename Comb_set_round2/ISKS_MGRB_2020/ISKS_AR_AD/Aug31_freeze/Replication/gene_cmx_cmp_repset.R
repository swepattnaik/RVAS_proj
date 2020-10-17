##gene_complexes
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

`%nin%` = Negate(`%in%`)
##all genes combined:positive control
# gene_all <- c("TP53", "NF1", "EXT1", "EXT2", "BRCA2", "ERCC2", "SDHA", "SDHB", "SDHD")
# gene_all <- c("TP53", "NF1")
# ##Sheltrin complex
# gene_sheltrin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "ACD")
# ##Extended Sheltrin complex
# Sheltrin_comp_extn = c("ACD", "POT1", "TERF1", "TERF2", "TERF2IP", "TINF2", "ATM", 
#                        "BAG3", "BLM", "BRCA1", "CALD1", "CLK3", "DCLRE1B", "FANCD2", 
#                        "FBXO4", "HSPA4", "KIAA1191", "MRE11A", "NBN", "PINX1", "PRKDC", 
#                        "RAD50", "SLX4", "STUB1", "TNKS", "TNKS2", "U2AF2", "UCHL1", 
#                        "WRN", "XRCC5", "XRCC6")
# ##HAUS5
# cmp_cep <- c("CEP63", "CEP72", "HAUS4", "HAUS5")
# ##GRIN2A
# cmp_grin <- c("DLG1", "GRIN2A", "GRIA1", "CAM2KB", "CAM2KD")
# ##ZBTB16
# novel <- c("ZBTB16", "ATG7", "ASB10", "ASB8", "FBXO7", "TRIM4")

Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "NBN", "SMARCAL1", "STAG3")

Glutamate_NMDA <- c("DLG1", "GRIN2A", "GRIA1", "CAM2KB", "CAM2KD")

ZBTB16_complex <- c("ZBTB16", "ATG7", "ASB10", "ASB8", "FBXO7")

CEP_HAUS_complex <- c("CEP63", "CEP72", "HAUS4", "HAUS5")

CENP_complex <- c("CENPC", "CENPF", "BUB1", "BUB1B", "AURKB", "RANGAP1", "RACGAP1", "MAD2L2")

TCP1_complex <- c("TCP1", "CCT2", "CCT6A")

PRPF4_complex <- c("PRPF4", "DHX9", "DHX16", "SF3A1", "CPSF3", "PPIL3")

Sarcoma_genes <- c("TP53", "NF1", "BRCA2", "ERCC2", "SDHA", "SDHB", "SDHD")

Breast_cancer_genes <- c("BRCA1", "BRCA2", "PALB2", "RAD51C", "CDH1")

cpx_genes <- unique(c(Shelterin, Glutamate_NMDA, ZBTB16_complex, CEP_HAUS_complex, CENP_complex, TCP1_complex,
                      PRPF4_complex, Breast_cancer_genes, Sarcoma_genes))

cpx_list <- list(Shelterin, Glutamate_NMDA, ZBTB16_complex, CEP_HAUS_complex, CENP_complex, TCP1_complex,
                 PRPF4_complex, Breast_cancer_genes, Sarcoma_genes)
names(cpx_list) <- c("Shelterin", "Glutamate_NMDA", "ZBTB16_complex", "CEP_HAUS", "CENP_complex", 
                     "TCP1_complex", "PRPF4_complex", "BRCA_genes", "SARC_genes")


##Odds ratio using fisher's exact test(use ppi_final_res... file for fisher's test)

#Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_para_pc1234_SKAT_Enriched_ISKS_2020_Aug31.rds")
#isks_mgrb_genes <- Exome_pc123_srt_SKAT_case_enr_nCH[[4]]
#isks_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Exome_para_tab_ISKS_MGRB_minusC3_Aug31.rds")

isks_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/isks_mgrb_all.rds")
isks_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/isks_mgrb_all_noC3.rds")
# ##Epithelial data
# epit_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/SKAT/Exome_para_MGRB_EPIT_2020.rds")
# epit_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/SKAT/Exome_para_MGRB_EPIT_2020_minusC3.rds")
# 
# ##ISKS_complex vs MGRB
# isks_complex_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Complex/Exome_para_ISKS_Complex_Aug31.rds")
# isks_complex_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Complex/Exome_para_ISKS_Complex_minusC3_Aug31.rds")
# 
# ##ISKS_TAS vs MGRB
# isks_tas_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/TAS/Exome_para_ISKS_TAS_Aug31.rds")
# isks_tas_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/TAS/Exome_para_ISKS_TAS_minusC3_Aug31.rds")
# ##ISKS_SIMPLE vs MGRB
# isks_simple_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Simple/Exome_para_ISKS_SIMPLE_2020_Aug31.rds")
# isks_simple_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Simple/Exome_para_ISKS_SIMPLE_2020_minusC3_Aug31.rds")
# 
# ##Combined Sarcoma replication set (not subsetted with PID genes)
# sarc_comb_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/Exome_para_CombSarc_Aug31.rds")
# sarc_comb_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/Exome_para_CombSarc_minusC3Aug31.rds")

##Augmented replication sets (the variants here have been trimmed based on the BED coordinates of exome kit)
## (subsetted with PID genes)
repset_563 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/repset_563.rds")
repset_563_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/repset_563_noC3.rds")

repset_663 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/repset_663.rds")
repset_663_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/repset_663_noC3.rds")

repset_763 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/repset_763.rds")
repset_763_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/repset_763_noC3.rds")

repset_963 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/repset_963.rds")
repset_963_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/repset_963_noC3.rds")

##ISKS complemented have the exome coordinates of discovery set NOT replication set
##ISKS noC3 complemented
isks_min100_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/isks_min100_noC3.rds")
isks_min200_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/isks_min200_noC3.rds")
isks_min400_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/isks_min400_noC3.rds")

##ISKS C345 complemented
isks_min100 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/isks_min100.rds")
isks_min200 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/isks_min200.rds")
isks_min400 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/isks_min400.rds")

########

cpx_OR_fisher <- function(ppi_res,case_coh_size, cont_coh_size, coh){
  ft_df <- list()
  for(i in 1:length(cpx_list)){
    ppi_res_tab <- ppi_res[ppi_res$gene %in% cpx_list[[i]],]
   # ppi_res_tab[,2] <- ifelse(ppi_res_tab[,2] == 0, 1, ppi_res_tab[,2])
    inp <- c(sum(ppi_res_tab[,1]), case_coh_size - sum(ppi_res_tab[,1]) , 
             sum(ppi_res_tab[,2]), cont_coh_size - sum(ppi_res_tab[,2]))
    sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
    colnames(sim_mat) <- c("case", "cont")
    rownames(sim_mat) <- c("hits", "no_hits")
    #ft <- fisher.test(sim_mat, alternative = "greater")
    ft <- fisher.test(sim_mat)
    ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.99)
    ft_df[[i]] <- cbind.data.frame("gene" = names(cpx_list)[i] ,"Cases" = sum(ppi_res_tab[,1]),
                                   "Controls" = sum(ppi_res_tab[,2]),
                                   "Fish_pval" = ft$p.value,"CI_lower" = ft$conf.int[1],
                                   "CI_upper" = ft$conf.int[2],
                                   "OR_Fish" = ft$estimate, "Coh" = coh)
  }
  return(ft_df)
}



cpx_ISKS_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_mgrb_genes, 1644, 3205, "ISKSvsMGRB"))
cpx_isks_min100 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_min100, 1544, 3205, "ISKS100vsMGRB"))
cpx_isks_min200 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_min200, 1444, 3205, "ISKS200vsMGRB"))
cpx_isks_min400 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_min400, 1244, 3205, "ISKS400vsMGRB"))

cpx_repset_563 <- do.call("rbind.data.frame", cpx_OR_fisher(repset_563, 563, 3205, "rep563vsMGRB"))
cpx_repset_663 <- do.call("rbind.data.frame", cpx_OR_fisher(repset_663, 663, 3205, "rep663vsMGRB"))
cpx_repset_763 <- do.call("rbind.data.frame", cpx_OR_fisher(repset_763, 763, 3205, "rep763vsMGRB"))
cpx_repset_963 <- do.call("rbind.data.frame", cpx_OR_fisher(repset_963, 963, 3205, "rep963vsMGRB"))

# cpx_ISKS_complex_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_complex_mgrb_genes, 832, 3205, "COMPLEXvsMGRB"))
# cpx_ISKS_tas_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_tas_mgrb_genes, 354, 3205, "TASvsMGRB"))
# cpx_ISKS_simple_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_simple_mgrb_genes, 370, 3205, "SIMPLEvsMGRB"))
# cpx_EPIT_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(epit_mgrb_genes, 842, 3205, "EPITvsMGRB"))
# #cpx_CAIRN_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(cairns_mgrb_genes, 413, 3205, "CAIRNvsMGRB"))
# cpx_REPSET_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(sarc_comb_mgrb_genes, 563, 3205, "REPvsMGRB"))
df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_OR_df, cpx_isks_min100, cpx_isks_min200, cpx_isks_min400,
                                 cpx_repset_563, cpx_repset_663, cpx_repset_763, cpx_repset_963)
df_comb_rnd2$Coh <- as.factor(gsub("vsMGRB", "", df_comb_rnd2$Coh))
df_comb_rnd2_filt <- df_comb_rnd2[df_comb_rnd2$Fish_pval < 0.05, ]

##ISKS_MGRB, EPIT_MGRB and CAIRNS_MGRB without C3.
##no_C3
cpx_ISKS_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_mgrb_genes_noC3, 1644, 3205, "ISKSvsMGRB"))
cpx_isks_min100_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_min100_noC3, 1544, 3205, "ISKS100vsMGRB"))
cpx_isks_min200_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_min200_noC3, 1444, 3205, "ISKS200vsMGRB"))
cpx_isks_min400_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_min400_noC3, 1244, 3205, "ISKS400vsMGRB"))

cpx_repset_563_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(repset_563_noC3, 563, 3205, "rep563vsMGRB"))
cpx_repset_663_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(repset_663_noC3, 663, 3205, "rep663vsMGRB"))
cpx_repset_763_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(repset_763_noC3, 763, 3205, "rep763vsMGRB"))
cpx_repset_963_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(repset_963_noC3, 963, 3205, "rep963vsMGRB"))


# cpx_ISKS_complex_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_complex_mgrb_genes_noC3, 832, 3205, "COMPLEXvsMGRB"))
# cpx_ISKS_tas_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_tas_mgrb_genes_noC3, 354, 3205, "TASvsMGRB"))
# cpx_ISKS_simple_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_simple_mgrb_genes_noC3, 370, 3205, "SIMPLEvsMGRB"))
# cpx_EPIT_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(epit_mgrb_genes_noC3, 842, 3205, "EPITvsMGRB"))
# cpx_REPSET_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(sarc_comb_mgrb_genes_noC3, 563, 3205, "REPvsMGRB"))
# #cpx_CAIRN_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(Cairns_noC3_fin, 413, 3205, "CAIRNvsMGRB"))
#noC3_df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_OR_df_noC3, cpx_ISKS_complex_OR_df_noC3, 
#                                      cpx_ISKS_tas_OR_df_noC3,cpx_ISKS_simple_OR_df_noC3,
#                                      cpx_EPIT_OR_df_noC3, cpx_CAIRN_OR_df_noC3)
##without Cairns
noC3_df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_OR_df_noC3, cpx_isks_min100_noC3, 
                                      cpx_isks_min200_noC3,cpx_isks_min400_noC3,
                                      cpx_repset_563_noC3, cpx_repset_663_noC3,
                                      cpx_repset_763_noC3, cpx_repset_963_noC3)
noC3_df_comb_rnd2$pval_adj <- p.adjust(noC3_df_comb_rnd2$Fish_pval, 
                                       n = length(noC3_df_comb_rnd2$Fish_pval), method = "bonferroni")

noC3_df_comb_rnd2$Coh <- as.factor(gsub("vsMGRB", "", noC3_df_comb_rnd2$Coh))

noC3_df_comb_rnd2$OR_Fish <- ifelse(noC3_df_comb_rnd2$OR_Fish == 0, 1, noC3_df_comb_rnd2$OR_Fish)
noC3_df_comb_rnd2_filt <- noC3_df_comb_rnd2[noC3_df_comb_rnd2$Fish_pval < 0.05 & noC3_df_comb_rnd2$CI_lower > 1,]
#noC3_df_comb_rnd2_filt <- noC3_df_comb_rnd2[noC3_df_comb_rnd2$pval_adj < 0.05 & noC3_df_comb_rnd2$CI_lower > 1,]
##Forest plot
#define colours for dots and bars
library(ggplot2)
forest_custplot <- function(df, repset = NULL){
#  dotCOLS = c("#a6d8f0","#f9b282", "#78f542", "#f5c6e5", "#ffe3bf", "#cfabff")
#  barCOLS = c("#008fd5","#de6b35", "#7a9406", "#fc0373", "#f79514", "#6c0ee8")
  df$Coh <- droplevels(df$Coh)
  ##with replication set (combined sarcomas)
  if(repset == 1){
  dotCOLS = c("#a6d8f0","#f9b282", "#78f542", "#f5c6e5","#cfabff", "#232533", "#e3c45d", "#05b362")
  barCOLS = c("#008fd5","#de6b35", "#7a9406", "#fc0373","#6c0ee8", "#a3ceff", "#704006", "#035730")
  }
  else{
    df <- df[df$Coh %nin% "REP",]
    ##without Cairns
    dotCOLS = c("#a6d8f0","#f9b282", "#78f542", "#f5c6e5","#cfabff")
    barCOLS = c("#008fd5","#de6b35", "#7a9406", "#fc0373","#6c0ee8")
  }
  p <- ggplot(df, aes(x=gene, y=log2(OR_Fish), ymin=log2(CI_lower), ymax=log2(CI_upper),col=Coh,fill=Coh)) + 
    #specify position here
    geom_linerange(size=1,position=position_dodge(width = 0.5)) +
    geom_hline(yintercept=0, lty=2) +
  #  geom_hline(yintercept=1, lty=2) +
    #specify position here too
    geom_point(size=5, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
    scale_fill_manual(values=barCOLS)+
    scale_color_manual(values=dotCOLS)+
    scale_x_discrete(name="Protein modules") + 
    scale_y_continuous(name="Log Odds ratio") +
    coord_flip() +
    theme_minimal() + theme(legend.position="bottom") + theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size=8))
  return(p)
}

#t1_df <- noC3_df_comb_rnd2[noC3_df_comb_rnd2$gene %nin% c("Sheltrin"),]
#forest_custplot(t1_df)
p1 <- forest_custplot(noC3_df_comb_rnd2, repset = 1) + ggtitle("C4_C5") + guides(fill=guide_legend(nrow=2,byrow=TRUE))
p2 <-forest_custplot(noC3_df_comb_rnd2_filt, repset = 1) + ggtitle("C4_C5_filtered") + guides(fill=guide_legend(nrow=2,byrow=TRUE))

p3 <-forest_custplot(df_comb_rnd2, repset = 1) + ggtitle("C3_C4_C5") + guides(fill=guide_legend(nrow=2,byrow=TRUE))
p4 <-forest_custplot(df_comb_rnd2_filt, repset = 1) + ggtitle("C3_C4_C5_filt") + guides(fill=guide_legend(nrow=2,byrow=TRUE))

source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")
#multiplot(p1, p2, p3, p4, cols = 2)
#multiplot(p4, p2, cols = 2)
multiplot(p1, p2, cols = 2)
multiplot(p3, p4, cols = 2)
