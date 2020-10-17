##gene_complexes
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
library(readxl)

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


##Complexes updated on Sep.27.2020
##Complexes updated on Sep.29.2020; remove NBN from Shelterin

#Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "NBN", "SMARCAL1", "STAG3", "TIMELESS")

Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")

Glutamate_NMDA <- c("DLG1", "GRIN2A", "GRIA1", "CAM2KB", "CAM2KD")

ZBTB16_complex <- c("ZBTB16", "ATG7", "ASB10", "ASB8", "FBXO7")

#CEP_HAUS_complex <- c("CEP63", "CEP72", "HAUS4", "HAUS5")
CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
CEP_HAUS_extn1 <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1", "CEP57", "PCM1")
#CEP_HAUS_extn2 <- c("MACC1", "SMPD2", "CEP63", "CEP72", "CEP89", "HAUS4", "HAUS5", "MZT1", "PCM1")


ppi_res_fil_final <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/ppi_res_fil_final_SKATbin_Aug31.tsv", sep = "\t", header = T)
CEP_interactors <- as.character(ppi_res_fil_final[ppi_res_fil_final$gene %in% CEP_HAUS_core,]$int)
CEP_interactors1 <- unique(unlist(lapply(CEP_interactors, function(x) strsplit(x, ","))))
CEP_HAUS_SKATint_complex <- CEP_interactors1

##CEP_HAUS_Biogrid_complex

library(igraph)
can_net <- read.delim("~/VDLab_scripts/BioGrid/biogrid_db_all_subnet.sif", header = T, sep = " ", stringsAsFactor = F)
#can_net <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Biogrid_latest/Biog_net_hs.sif", header = T, sep = "\t", stringsAsFactor = F)
can_net_graph <- igraph::graph.data.frame(can_net, directed = F)
can_net_graph1 <- igraph::simplify(can_net_graph, remove.loops=T, remove.multiple = T)
can_net1 <- igraph::as_data_frame(can_net_graph1, what = "edges")
prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% CEP_HAUS_core
| as.character(can_net1$to) %in% CEP_HAUS_core, ])
CEP_HAUS_Biogrid_complex <- unique(c(prot_np_all$from, prot_np_all$to))

CENP_complex <- c("CENPC", "CENPF", "BUB1", "BUB1B", "AURKB", "RANGAP1", "RACGAP1", "MAD2L2")

TCP1_complex <- c("TCP1", "CCT2", "CCT6A")

PRPF4_complex <- c("PRPF4", "DHX9", "DHX16", "SF3A1", "CPSF3", "PPIL3")

#Sarcoma_genes <- c("TP53", "NF1", "BRCA2", "ERCC2", "SDHA", "SDHB", "SDHD")

#Sarcoma_genes <- c("TP53", "NF1", "SDHA", "SDHB", "SDHD")

Sarcoma_genes <- c("TP53", "SDHA", "SDHB", "SDHD")

TP53_control <- c("TP53") #positive control

MPNST_pos <- c("NF1", "LZTR1", "SDHA", "SDHB", "SDHD")

Breast_cancer_genes <- c("BRCA1", "BRCA2", "PALB2", "RAD51C", "CDH1")

Centrosome_genes <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/Centrosome maturation.xlsx",
                               sheet = 1)
Centrosome_genes <- as.data.frame(Centrosome_genes)
Centrosome_genes <- Centrosome_genes[Centrosome_genes$MoleculeType %nin% "Chemical Compounds",]

Centrosome_genes <- Centrosome_genes[,4]

mito_cell_cycle <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/Mitotic_cell_cycle_genes_HSA_69278.txt", header = T, stringsAsFactors = F)
mito_HSA_69278 <- unique(as.character(mito_cell_cycle$X))

##mitotic check point genes (GO term)
mito_chk_point = read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/mitotic_checkpoint/mitocheck_point.txt",
                            sep = "", header = F, skip = 1, stringsAsFactors = F)
mito_chk_point <- as.character(mito_chk_point$V1)

##mito overlap
mito_43_overlap <- mito_chk_point[mito_chk_point %in% mito_HSA_69278]

#Centriole replication (GO Term)
cent_rep <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/centriole_replication_genes_GO.xlsx",
                       sheet = 1)
cent_rep <- as.data.frame(cent_rep)
cent_rep_genes <- toupper(unique(cent_rep$Symbol)) #40 genes

cpx_genes <- unique(c(Shelterin, Glutamate_NMDA, ZBTB16_complex, CEP_HAUS_core, CEP_HAUS_extn1, CEP_HAUS_SKATint_complex, 
                      CEP_HAUS_Biogrid_complex, CENP_complex, TCP1_complex,PRPF4_complex, Breast_cancer_genes, Sarcoma_genes, 
                      Centrosome_genes, mito_HSA_69278, mito_chk_point, mito_43_overlap, cent_rep_genes, TP53_control, MPNST_pos))

cpx_list <- list(Shelterin, Glutamate_NMDA, ZBTB16_complex, CEP_HAUS_core, CEP_HAUS_extn1, CEP_HAUS_SKATint_complex, 
                 CEP_HAUS_Biogrid_complex, CENP_complex, TCP1_complex,PRPF4_complex, Breast_cancer_genes, Sarcoma_genes, 
                 Centrosome_genes, mito_HSA_69278, mito_chk_point, mito_43_overlap, cent_rep_genes, TP53_control, MPNST_pos)
names(cpx_list) <- c("Shelterin", "Glutamate_NMDA", "ZBTB16_complex", "CEP_HAUS_core", "CEP_HAUS_extn1", 
                     "CEP_HAUS_SKAT", "CEP_HAUS_Biogrid", "CENP_complex", "TCP1_complex", "PRPF4_complex", "BRCA_genes", 
                     "SARC_genes", "Centrosome_genes", "mito_HSA_69278", "mito_chk_GO", "mito_43", "cent_rep_genes",
                     "TP53_control", "MPNST_pos")


##Odds ratio using fisher's exact test(use ppi_final_res... file for fisher's test)

#Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_para_pc1234_SKAT_Enriched_ISKS_2020_Aug31.rds")
#isks_mgrb_genes <- Exome_pc123_srt_SKAT_case_enr_nCH[[4]]
isks_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_para_tab_ISKS_MGRB_C345_Aug31.rds")
isks_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Exome_para_tab_ISKS_MGRB_minusC3_Aug31.rds")

##Epithelial data
epit_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/SKAT/Exome_para_MGRB_EPIT_2020.rds")
epit_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/SKAT/Exome_para_MGRB_EPIT_2020_minusC3.rds")

##Epithelial data (minus melanoma)
nmela_epit_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/SKAT/Exome_para_EPIT_nmela_2020.rds")
nmela_epit_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/SKAT/Exome_para_EPIT_nmela_2020_noC3.rds")

##Epithelial data (melanoma)
mela_epit_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/SKAT/Exome_para_EPIT_mela_2020.rds")
mela_epit_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/SKAT/Exome_para_EPIT_mela_2020_noC3.rds")


##ISKS_complex vs MGRB
isks_complex_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Complex/Exome_para_ISKS_Complex_Aug31.rds")
isks_complex_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Complex/Exome_para_ISKS_Complex_minusC3_Aug31.rds")

##ISKS_TAS vs MGRB
isks_tas_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/TAS/Exome_para_ISKS_TAS_Aug31.rds")
isks_tas_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/TAS/Exome_para_ISKS_TAS_minusC3_Aug31.rds")
##ISKS_SIMPLE vs MGRB
isks_simple_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Simple/Exome_para_ISKS_SIMPLE_2020_Aug31.rds")
isks_simple_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Simple/Exome_para_ISKS_SIMPLE_2020_minusC3_Aug31.rds")

##ISKS_MPNST vs MGRB
isks_mpnst_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/MPNST/Exome_para_ISKS_MPNST_sep30.rds")
isks_mpnst_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/MPNST/Exome_para_ISKS_MPNST_sep30_noC3.rds")


##Combined Sarcoma replication set
sarc_comb_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/Exome_para_CombSarc_Aug31.rds")
sarc_comb_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/Exome_para_CombSarc_minusC3Aug31.rds")

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
    #ft <- fisher.test(sim_mat)
    ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.99)
    ft_df[[i]] <- cbind.data.frame("gene" = names(cpx_list)[i] ,"Cases" = sum(ppi_res_tab[,1]),
                                   "Controls" = sum(ppi_res_tab[,2]),
                                   "Fish_pval" = ft$p.value,"CI_lower" = ft$conf.int[1],
                                   "CI_upper" = ft$conf.int[2],
                                   "OR_Fish" = ft$estimate, "case_coh_size" = case_coh_size,
                                     "Coh" = coh)
  }
  return(ft_df)
}
##removed 7 oesphageal samples (Sep30)
cpx_ISKS_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_mgrb_genes, 1644, 3205, "ISKSvsMGRB"))
cpx_ISKS_complex_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_complex_mgrb_genes, 837, 3205, "COMPLEXvsMGRB"))
cpx_ISKS_tas_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_tas_mgrb_genes, 362, 3205, "TASvsMGRB"))
cpx_ISKS_simple_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_simple_mgrb_genes, 371, 3205, "SIMPLEvsMGRB"))
cpx_ISKS_mpnst_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_mpnst_mgrb_genes, 157, 3205, "MPNSTvsMGRB"))
cpx_EPIT_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(epit_mgrb_genes, 835, 3205, "EPITvsMGRB"))
cpx_NMELA_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(nmela_epit_mgrb_genes, 660, 3205, "NMELAvsMGRB"))
cpx_MELA_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(mela_epit_mgrb_genes, 175, 3205, "MELAvsMGRB"))
#cpx_CAIRN_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(cairns_mgrb_genes, 413, 3205, "CAIRNvsMGRB"))
cpx_REPSET_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(sarc_comb_mgrb_genes, 563, 3205, "REPvsMGRB"))
df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_OR_df, cpx_ISKS_complex_OR_df, cpx_ISKS_tas_OR_df, 
                                 cpx_ISKS_simple_OR_df, cpx_ISKS_mpnst_OR_df, 
                                 cpx_EPIT_OR_df,cpx_NMELA_OR_df,
                                 cpx_MELA_OR_df, cpx_REPSET_OR_df)
df_comb_rnd2$Coh <- as.factor(gsub("vsMGRB", "", df_comb_rnd2$Coh))
df_comb_rnd2_filt <- df_comb_rnd2[df_comb_rnd2$Fish_pval < 0.05, ]

##ISKS_MGRB, EPIT_MGRB and CAIRNS_MGRB without C3.
##no_C3
cpx_ISKS_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_mgrb_genes_noC3, 1644, 3205, "ISKSvsMGRB"))
cpx_ISKS_complex_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_complex_mgrb_genes_noC3, 837, 3205, "COMPLEXvsMGRB"))
cpx_ISKS_tas_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_tas_mgrb_genes_noC3, 362, 3205, "TASvsMGRB"))
cpx_ISKS_simple_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_simple_mgrb_genes_noC3, 371, 3205, "SIMPLEvsMGRB"))
cpx_ISKS_mpnst_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_mpnst_mgrb_genes_noC3, 157, 3205, "MPNSTvsMGRB"))
cpx_EPIT_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(epit_mgrb_genes_noC3, 835, 3205, "EPITvsMGRB"))
cpx_NMELA_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(nmela_epit_mgrb_genes_noC3, 660, 3205, "NMELAvsMGRB"))
cpx_MELA_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(mela_epit_mgrb_genes_noC3, 175, 3205, "MELAvsMGRB"))

cpx_REPSET_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(sarc_comb_mgrb_genes_noC3, 563, 3205, "REPvsMGRB"))
#cpx_CAIRN_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(Cairns_noC3_fin, 413, 3205, "CAIRNvsMGRB"))
#noC3_df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_OR_df_noC3, cpx_ISKS_complex_OR_df_noC3, 
#                                      cpx_ISKS_tas_OR_df_noC3,cpx_ISKS_simple_OR_df_noC3,
#                                      cpx_EPIT_OR_df_noC3, cpx_CAIRN_OR_df_noC3)
##without Cairns
noC3_df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_OR_df_noC3, cpx_ISKS_complex_OR_df_noC3, 
                                      cpx_ISKS_tas_OR_df_noC3, cpx_ISKS_simple_OR_df_noC3, 
                                      cpx_ISKS_mpnst_OR_df_noC3, cpx_EPIT_OR_df_noC3, 
                                      cpx_NMELA_OR_df_noC3, cpx_MELA_OR_df_noC3,
                                      cpx_REPSET_OR_df_noC3)
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
  dotCOLS = c("#a6d8f0","#f9b282", "#78f542", "#f5c6e5","#cfabff", "#232533", "#bdb071", "#f52a6b", "#2a7ff5")
  barCOLS = c("#008fd5","#de6b35", "#7a9406", "#fc0373","#6c0ee8", "#a3ceff", "#717bbd", "#8f5467", "#f5822a")
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

##intercohort comparison
tab_noC3 <- noC3_df_comb_rnd2
tab_C345 <- df_comb_rnd2

inter_coh_fisher <- function(tab_dat, cpx_name, case_coh, cont_coh){
  cohs <- c(case_coh, cont_coh)
  df_sel <- tab_dat[tab_dat$gene %in% cpx_name & tab_dat$Coh %in% c(case_coh, cont_coh),c(2,8:9)]
  df_sel <- df_sel[match(cohs, df_sel$Coh),]
  inp <- c(df_sel[1,1], df_sel[1,2] - df_sel[1,1] , 
           df_sel[2,1], df_sel[2,2] - df_sel[2,1])
  inp[3] <- ifelse(inp[3] == 0, 1, inp[3])
  sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
  colnames(sim_mat) <- c("case", "cont")
  rownames(sim_mat) <- c("hits", "no_hits")
  #ft <- fisher.test(sim_mat, alternative = "greater")
  #ft <- fisher.test(sim_mat)
  ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.99)
  ft_df <- cbind.data.frame("gene" = cpx_name ,"Cases" = df_sel[1,1],
                                 "Controls" = df_sel[2,1], "Case_Coh" = df_sel[1,2],
                                 "Cont_Coh" = df_sel[2,2], "Fish_pval" = ft$p.value,
                                 "CI_lower" = ft$conf.int[1],
                                 "CI_upper" = ft$conf.int[2],
                                 "OR_Fish" = ft$estimate, "Coh" = paste(cohs, collapse = "vs"))
  return(ft_df)
}

cont_cohs <- c("EPIT", "MELA", "NMELA")
case_cohs <- c("ISKS", "COMPLEX", "SIMPLE", "TAS", "MPNST")


multi_inter_coh <- function(tab, cases, controls, cpx_name){
all_comp <- list()
all_comp_sub <- list()
for(m in 1:length(controls)){
  for(n in 1:length(cases)){
    all_comp_sub[[n]] <- inter_coh_fisher(tab_dat = tab, cpx_name = cpx_name,  case_coh = cases[n], 
                                          cont_coh = controls[m])
  }
  all_comp[[m]] <- do.call("rbind.data.frame", all_comp_sub)
}
all_comp_df <- do.call("rbind.data.frame", all_comp)
return(all_comp_df)
}
##C45
Sheltrin_cmp <- multi_inter_coh(tab_noC3, case_cohs, cont_cohs, "Shelterin")
CEP_HAUS_core_cmp <- multi_inter_coh(tab_noC3, case_cohs, cont_cohs, "CEP_HAUS_core")
CEP_HAUS_extn1_cmp <- multi_inter_coh(tab_noC3, case_cohs, cont_cohs, "CEP_HAUS_extn1")
BRCA_cmp <- multi_inter_coh(tab_noC3, case_cohs, cont_cohs, "BRCA_genes")
TP53_cmp <- multi_inter_coh(tab_noC3, case_cohs, cont_cohs, "TP53_control")
fin_comb_tab <- rbind.data.frame(Sheltrin_cmp, CEP_HAUS_core_cmp, CEP_HAUS_extn1_cmp, BRCA_cmp, TP53_cmp)
fin_comb_tab$sig <- ifelse(fin_comb_tab$Fish_pval > 0.05, "ns", ifelse(fin_comb_tab$Fish_pval > 0.01 & fin_comb_tab$Fish_pval <= 0.05, "*",
                           ifelse(fin_comb_tab$Fish_pval > 0.001 & fin_comb_tab$Fish_pval <= 0.01, "**",
                           ifelse(fin_comb_tab$Fish_pval > 0.0001 & fin_comb_tab$Fish_pval <= 0.001, "***", "****"))))
#fin_comb_tab1 <- fin_comb_tab
#fin_comb_tab1[fin_comb_tab1$Fish_pval > 0.1,]$OR_Fish <- 0.1 
##heatmap

##Finalised Oct.8 (use OR values for heatmap and use p-values to represent significance)
library(viridis)
# ggplot2(fin_comb_tab, aes(Coh, gene, fill= OR_Fish)) + 
#   geom_tile() + scale_fill_viridis(discrete=FALSE) + theme(axis.text.x = element_text(angle = 45))

ggplot(fin_comb_tab, aes(Coh, gene, fill= log(OR_Fish))) + geom_tile() + scale_fill_distiller(palette = "RdBu") + 
  theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.title.x=element_blank()) + 
  theme(axis.text.x = element_text(margin = margin(t = 20))) + 
  theme(legend.position="bottom") + geom_text(aes(label = sig))

# ggplot(fin_comb_tab1, aes(Coh, gene, fill= log2(OR_Fish))) + geom_tile() + scale_fill_distiller(palette = "RdBu") + 
#   theme(axis.text.x = element_text(angle = 45)) + 
#   theme(axis.title.x=element_blank()) + 
#   theme(axis.text.x = element_text(margin = margin(t = 20))) + 
#   theme(legend.position="bottom") 

ggplot(fin_comb_tab, aes(Coh, gene, fill= log(OR_Fish))) + geom_tile() + scale_fill_distiller(palette = "RdBu") + 
  theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.title.x=element_blank()) + 
  theme(axis.text.x = element_text(margin = margin(t = 20))) + 
  theme(legend.position="bottom") + geom_text(aes(label = sig))
  
################heatmap (ISKS vs MGRB)
tab_noC3$Coh <- factor(tab_noC3$Coh, levels = unique(tab_noC3$Coh))
tab_noC3$sig <- ifelse(tab_noC3$Fish_pval > 0.05, "ns", ifelse(tab_noC3$Fish_pval > 0.01 & tab_noC3$Fish_pval <= 0.05, "*",
                                                               ifelse(tab_noC3$Fish_pval > 0.001 & tab_noC3$Fish_pval <= 0.01, "**",
                                                                      ifelse(tab_noC3$Fish_pval > 0.0001 & tab_noC3$Fish_pval <= 0.001, "***", "****"))))
p1 <- ggplot(tab_noC3, aes(Coh, gene, fill= log(OR_Fish))) + geom_tile() + theme_minimal()
p1 + scale_fill_distiller(palette = "RdBu") + 
  theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.title.x=element_blank()) + 
  theme(axis.text.x = element_text(margin = margin(t = 20))) + 
  theme(legend.position="bottom") + geom_text(aes(label = sig))

##Add dendrogram to heatmap
library("ggplot2")
library("ggdendro")
library("reshape2")
library("grid")
##C45 subset complex and cohorts
tab_noC3_sub <- tab_noC3[tab_noC3$Coh %nin% c("TAS", "SIMPLE", "EPIT", "REP") &
                           tab_noC3$gene %in% c("Shelterin", "CEP_HAUS_core", 
                                                "SARC_genes", "BRCA_genes", 
                                                "TP53_control", "MPNST_pos"),]
tab_noC3_sub$Coh <- droplevels(tab_noC3_sub$Coh)
tab_noC3_sub$gene <- droplevels(tab_noC3_sub$gene)
tab_noC3_sub$Coh <- gsub("MPNST", "MPNST+GIST",tab_noC3_sub$Coh)
tab_noC3_sub$Coh <- factor(gsub("NMELA", "EPIT",tab_noC3_sub$Coh))
tab_noC3_sub$Coh <- factor(x = tab_noC3_sub$Coh,
                            levels =c("ISKS", "COMPLEX", "MPNST+GIST",
                                      "EPIT", "MELA") , 
                            ordered = TRUE)
df_dendro <- tab_noC3_sub[,c(1,9,7)]
df_dendro$OR_Fish <- ifelse(is.infinite(df_dendro$OR_Fish), 0, df_dendro$OR_Fish) 
df_dendro.scaled <- df_dendro
df_dendro.scaled[, 3] <- scale(df_dendro.scaled[, 3])
df_dendro_mat <- acast(df_dendro.scaled, gene ~ Coh)
dend_inp <- as.dendrogram(hclust(d = dist(x = df_dendro_mat))) 

# Create dendro
dendro.plot <- ggdendrogram(data = dend_inp, rotate = TRUE) + 
  theme(axis.text.y = element_text(size = 9))

# Preview the plot
print(dendro.plot)

##order heatmap
dend_inp.order <- order.dendrogram(dend_inp)
tab_noC3_sub$gene <- factor(x = tab_noC3_sub$gene,
                               levels = tab_noC3_sub$gene[dend_inp.order], 
                               ordered = TRUE)
##C45
p1 <- ggplot(tab_noC3_sub, aes(Coh, gene, fill= log(OR_Fish))) + geom_tile() + theme_minimal()
p1_order <- p1 + scale_fill_distiller(palette = "RdBu") + 
  theme(axis.text.x = element_text(angle = 0)) + 
  theme(axis.title.x=element_blank()) + 
  theme(axis.text.x = element_text(margin = margin(t = 5))) + 
  theme(legend.position="top") + geom_text(aes(label = sig)) + 
  theme(axis.text.y = element_blank()) + theme(axis.title.y=element_blank())
grid.newpage()
print(p1_order, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.85, y = 0.4, width = 0.18, height = 0.82))
dev.off()

##figure heatmap

################C345
Sheltrin_cmp <- multi_inter_coh(tab_C345, case_cohs, cont_cohs, "Shelterin")
CEP_HAUS_core_cmp <- multi_inter_coh(tab_C345, case_cohs, cont_cohs, "CEP_HAUS_core")
CEP_HAUS_extn1_cmp <- multi_inter_coh(tab_C345, case_cohs, cont_cohs, "CEP_HAUS_extn1")
BRCA_cmp <- multi_inter_coh(tab_C345, case_cohs, cont_cohs, "BRCA_genes")
TP53_cmp <- multi_inter_coh(tab_C345, case_cohs, cont_cohs, "TP53_control")
fin_comb_tab_C345 <- rbind.data.frame(Sheltrin_cmp, CEP_HAUS_core_cmp, CEP_HAUS_extn1_cmp, BRCA_cmp, TP53_cmp)
fin_comb_tab_C345$sig <- ifelse(fin_comb_tab_C345$Fish_pval > 0.05, "ns", ifelse(fin_comb_tab_C345$Fish_pval > 0.01 & fin_comb_tab_C345$Fish_pval <= 0.05, "*",
                                                                       ifelse(fin_comb_tab_C345$Fish_pval > 0.001 & fin_comb_tab_C345$Fish_pval <= 0.01, "**",
                                                                              ifelse(fin_comb_tab_C345$Fish_pval > 0.0001 & fin_comb_tab_C345$Fish_pval <= 0.001, "***", "****"))))
#fin_comb_tab1 <- fin_comb_tab
#fin_comb_tab1[fin_comb_tab1$Fish_pval > 0.1,]$OR_Fish <- 0.1 
##heatmap

##Finalised Oct.8 (use OR values for heatmap and use p-values to represent significance)
library(viridis)

ggplot(fin_comb_tab_C345, aes(Coh, gene, fill= log(OR_Fish))) + geom_tile() + scale_fill_distiller(palette = "RdBu") + 
  theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.title.x=element_blank()) + 
  theme(axis.text.x = element_text(margin = margin(t = 20))) + 
  theme(legend.position="bottom") + geom_text(aes(label = sig))

###vs MGRB
tab_C345$Coh <- factor(tab_C345$Coh, levels = unique(tab_C345$Coh))
tab_C345$sig <- ifelse(tab_C345$Fish_pval > 0.05, "ns", ifelse(tab_C345$Fish_pval > 0.01 & tab_C345$Fish_pval <= 0.05, "*",
                                                               ifelse(tab_C345$Fish_pval > 0.001 & tab_C345$Fish_pval <= 0.01, "**",
                                                                      ifelse(tab_C345$Fish_pval > 0.0001 & tab_C345$Fish_pval <= 0.001, "***", "****"))))
df_dendro <- tab_C345[,c(1,9,7)]
df_dendro$OR_Fish <- ifelse(is.infinite(df_dendro$OR_Fish), 0, df_dendro$OR_Fish) 
df_dendro.scaled <- df_dendro
df_dendro.scaled[, 3] <- scale(df_dendro.scaled[, 3])
df_dendro_mat <- acast(df_dendro.scaled, gene ~ Coh)
dend_inp <- as.dendrogram(hclust(d = dist(x = df_dendro_mat))) 

# Create dendro
dendro.plot <- ggdendrogram(data = dend_inp, rotate = TRUE) + 
  theme(axis.text.y = element_text(size = 10))

# Preview the plot
print(dendro.plot)

##order heatmap
dend_inp.order <- order.dendrogram(dend_inp)
tab_C345$gene <- factor(x = tab_C345$gene,
                        levels = tab_C345$gene[dend_inp.order], 
                        ordered = TRUE)
##C345
p1 <- ggplot(tab_C345, aes(Coh, gene, fill= log(OR_Fish))) + geom_tile() + theme_minimal()
p1_order <- p1 + scale_fill_distiller(palette = "RdBu") + 
  theme(axis.text.x = element_text(angle = 0)) + 
  theme(axis.title.x=element_blank()) + 
  theme(axis.text.x = element_text(margin = margin(t = 10))) + 
  theme(legend.position="top") + geom_text(aes(label = sig)) + 
  theme(axis.text.y = element_blank()) + theme(axis.title.y=element_blank())
grid.newpage()
print(p1_order, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.9, y = 0.44, width = 0.18, height = 0.94))
dev.off()


##Bubble plots
