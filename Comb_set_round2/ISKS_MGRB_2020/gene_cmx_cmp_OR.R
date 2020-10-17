##gene_complexes
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
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

#CEP_HAUS_complex <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "PCM1")

CEP_HAUS_complex <- c("CEP63", "CEP72", "HAUS4", "HAUS5")

CENP_complex <- c("CENPC", "CENPF", "BUB1", "BUB1B", "AURKB", "RANGAP1", "RACGAP1", "MAD2L2")

TCP1_complex <- c("TCP1", "CCT2", "CCT6A")

PRPF4_complex <- c("PRPF4", "DHX9", "DHX16", "SF3A1", "CPSF3", "PPIL3")

Sarcoma_genes <- c("TP53", "NF1", "BRCA2", "ERCC2", "SDHA", "SDHB", "SDHD")

Breast_cancer_genes <- c("BRCA1", "BRCA2", "PALB2", "RAD51C", "CDH1")

cust_genes <- c("OS9", "NR1I3", "LINGO1", "LINGO4", "L3MBTL3")

Glutamate_NMDA_comb <- c(Glutamate_NMDA, cust_genes)

pos_TP53 <- c("TP53")

##Odds ratio using fisher's exact test(use ppi_final_res... file for fisher's test)

Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Exome_para_pc123_SKAT_Enriched_ISKS_2020.rds")
isks_mgrb_genes <- Exome_pc123_srt_SKAT_case_enr_nCH[[4]]

##Epithelial data
Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/round2/Exome_para_pc123_SKAT_Enriched_EPIT_2020_no_fbrca.rds")
epit_mgrb_genes <- Exome_pc123_srt_SKAT_case_enr_nCH[[4]]

##Cairns data
Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CAIRNS/Exome_para_pc123_SKAT_Enriched_CAIRNS_2020.rds")
cairns_mgrb_genes <- Exome_pc123_srt_SKAT_case_enr_nCH[[4]]

##ISKS_complex vs MGRB
isks_complex_mgrb <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Exome_para_pc123_SKAT_Enriched_ISKS_2020_complex.rds")
isks_complex_mgrb_genes <- isks_complex_mgrb[[4]]

##ISKS_TAS vs MGRB
isks_tas_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Exome_para_ISKS_TAS_2020.rds")

##ISKS_SIMPLE vs MGRB
isks_simple_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Exome_para_ISKS_SIMPLE_2020.rds")
########

#isks_epit_cairn <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_EPIT_cmp/round2/isks_ppi_epit_no_fbrca_cairns_comb_genes_top.tsv", 
#                              sep = "\t", header = T, stringsAsFactors = F)
##use only enriched(SKATBin <0.1) genes for sheltrin extn complex, for sarc genes use TP53, NF1.
#Sheltrin_comp_extn_enr <- isks_epit_cairn[isks_epit_cairn$gene %in% Sheltrin_comp_extn & isks_epit_cairn$pval_SKATbin < 0.1,]$gene

##Exact test OR and Pvalues
#gene_all <- c("TP53", "NF1")
#cpx_list <- list(gene_all,gene_sheltrin, Sheltrin_comp_extn, Sheltrin_comp_extn_enr, cmp_cep, cmp_grin, novel)
#names(cpx_list) <- c("Sarc_genes", "Sheltrin", "Sheltrin_extn", "Sheltrin_enr", "HAUS5_cmp", "GRIN_cmp", "ZBTB16_cmp")

cpx_list <- list(Shelterin, Glutamate_NMDA, ZBTB16_complex, CEP_HAUS_complex, CENP_complex, TCP1_complex,
                 PRPF4_complex, Breast_cancer_genes, Sarcoma_genes, cust_genes, Glutamate_NMDA_comb, pos_TP53)
names(cpx_list) <- c("Shelterin", "Glutamate_NMDA", "ZBTB16_complex", "CEP_HAUS", "CENP_complex", 
                     "TCP1_complex", "PRPF4_complex", "BRCA_genes", "SARC_genes", "cust_genes", "Glutamate_NMDA_comb" ,"pos_TP53")

cpx_OR_fisher <- function(ppi_res,case_coh_size, cont_coh_size, coh){
  ft_df <- list()
  for(i in 1:length(cpx_list)){
    ppi_res_tab <- ppi_res[ppi_res$gene %in% cpx_list[[i]],]
    ppi_res_tab[,2] <- ifelse(ppi_res_tab[,2] == 0, 1, ppi_res_tab[,2])
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

cpx_ISKS_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_mgrb_genes, 1646, 3205, "ISKSvsMGRB"))
cpx_ISKS_complex_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_complex_mgrb_genes, 832, 3205, "COMPLEXvsMGRB"))
cpx_ISKS_tas_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_tas_mgrb_genes, 354, 3205, "TASvsMGRB"))
cpx_ISKS_simple_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_simple_mgrb_genes, 370, 3205, "SIMPLEvsMGRB"))
cpx_EPIT_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(epit_mgrb_genes, 842, 3205, "EPITvsMGRB"))
cpx_CAIRN_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(cairns_mgrb_genes, 413, 3205, "CAIRNvsMGRB"))
df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_OR_df, cpx_ISKS_complex_OR_df, cpx_EPIT_OR_df, cpx_CAIRN_OR_df)
df_comb_rnd2_filt <- df_comb_rnd2[df_comb_rnd2$Fish_pval < 0.05, ]

##make_tab for ISKS_MGRB, EPIT_MGRB and CAIRNS_MGRB without C3.
##CAIRNS
#Cairns_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CAIRNS/Exome_para_pc123_SKAT_Enriched_CAIRNS_2020_minusC3.rds")
#Cairns_noC3_fin <- Cairns_noC3[[4]]
isks_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Exome_para_pc123_SKAT_Enriched_ISKS_2020_minusC3.rds")
isks_noC3_fin <- isks_noC3[[4]]
isks_cpx_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Exome_para_pc123_SKAT_Enriched_ISKScomplex_2020_minusC3.rds")
isks_cpx_noC3_fin <- isks_cpx_noC3[[4]]
epit_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT/round2/Exome_para_pc123_SKAT_Enriched_EPIT_2020_no_fbrca_minusC3.rds")
epit_noC3_fin <- epit_noC3[[4]]
isks_tas_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Exome_para_ISKS_TAS_2020_minusC3.rds")
isks_simple_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Exome_para_ISKS_SIMPLE_2020_minusC3.rds")
##no_C3
cpx_ISKS_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_noC3_fin, 1646, 3205, "ISKSvsMGRB"))
cpx_ISKS_complex_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_cpx_noC3_fin, 832, 3205, "COMPLEXvsMGRB"))
cpx_ISKS_tas_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_tas_noC3, 354, 3205, "TASvsMGRB"))
cpx_ISKS_simple_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_simple_noC3, 370, 3205, "SIMPLEvsMGRB"))
cpx_EPIT_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(epit_noC3_fin, 842, 3205, "EPITvsMGRB"))
#cpx_CAIRN_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(Cairns_noC3_fin, 413, 3205, "CAIRNvsMGRB"))
#noC3_df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_OR_df_noC3, cpx_ISKS_complex_OR_df_noC3, 
#                                      cpx_ISKS_tas_OR_df_noC3,cpx_ISKS_simple_OR_df_noC3,
#                                      cpx_EPIT_OR_df_noC3, cpx_CAIRN_OR_df_noC3)
##without Cairns
noC3_df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_OR_df_noC3, cpx_ISKS_complex_OR_df_noC3, 
                                      cpx_ISKS_tas_OR_df_noC3,cpx_ISKS_simple_OR_df_noC3,
                                      cpx_EPIT_OR_df_noC3)
noC3_df_comb_rnd2$pval_adj <- p.adjust(noC3_df_comb_rnd2$Fish_pval, 
                                       n = length(noC3_df_comb_rnd2$Fish_pval), method = "bonferroni")

noC3_df_comb_rnd2$Coh <- as.factor(gsub("vsMGRB", "", noC3_df_comb_rnd2$Coh))

noC3_df_comb_rnd2$OR_Fish <- ifelse(noC3_df_comb_rnd2$OR_Fish == 0, 1, noC3_df_comb_rnd2$OR_Fish)
noC3_df_comb_rnd2_filt <- noC3_df_comb_rnd2[noC3_df_comb_rnd2$Fish_pval < 0.05 & noC3_df_comb_rnd2$CI_lower > 1,]
#noC3_df_comb_rnd2_filt <- noC3_df_comb_rnd2[noC3_df_comb_rnd2$pval_adj < 0.05 & noC3_df_comb_rnd2$CI_lower > 1,]
##Forest plot
#define colours for dots and bars
library(ggplot2)
forest_custplot <- function(df){
#  dotCOLS = c("#a6d8f0","#f9b282", "#78f542", "#f5c6e5", "#ffe3bf", "#cfabff")
#  barCOLS = c("#008fd5","#de6b35", "#7a9406", "#fc0373", "#f79514", "#6c0ee8")
  ##without Cairns
  dotCOLS = c("#a6d8f0","#f9b282", "#78f542", "#f5c6e5","#cfabff")
  barCOLS = c("#008fd5","#de6b35", "#7a9406", "#fc0373","#6c0ee8")
  
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

`%nin%` = Negate(`%in%`)
#t1_df <- noC3_df_comb_rnd2[noC3_df_comb_rnd2$gene %nin% c("Sheltrin"),]
#forest_custplot(t1_df)
p1 <- forest_custplot(noC3_df_comb_rnd2) + ggtitle("C4_C5") + guides(fill=guide_legend(nrow=2,byrow=TRUE))
p2 <-forest_custplot(noC3_df_comb_rnd2_filt) + ggtitle("C4_C5_filtered") + guides(fill=guide_legend(nrow=2,byrow=TRUE))

p3 <-forest_custplot(df_comb_rnd2)
p4 <-forest_custplot(df_comb_rnd2_filt) + ggtitle("C3_C4_C5")

source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")
#multiplot(p1, p2, p3, p4, cols = 2)
#multiplot(p4, p2, cols = 2)
multiplot(p1, p2, cols = 2)
##ISKS_complex, ISKS_simple, ISKS_TAS; use fil_tab, PID file (genomic class)
##use make_tab to generate corresponding groups



##for talk : Maya cust_genes
#df_comb_rnd2_fin <- df_comb_rnd2[df_comb_rnd2$gene %in% c("Sarc_genes", "Sheltrin",
#                                                          "HAUS5_cmp", "GRIN_cmp", "ZBTB16_cmp"),]
df_comb_rnd2_fin_filt <- noC3_df_comb_rnd2_filt[noC3_df_comb_rnd2_filt$gene %in% c("Shelterin", "Glutamate_NMDA", "ZBTB16_complex", "CEP_HAUS", "CENP_complex", 
                                                                                   "TCP1_complex", "PRPF4_complex", "BRCA_genes", "SARC_genes", "Glutamate_NMDA_comb", "cust_genes"),]
# p5 <-forest_custplot(df_comb_rnd2_fin)
p6 <-forest_custplot(df_comb_rnd2_fin_filt)
multiplot(p5, p6, cols = 2)

#########non coding OR######
##non-coding obtained from Exome data
isks_tas_nc <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/nc_Exome_para_ISKS_TAS_2020_C2.rds")
isks_simple_nc <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/nc_Exome_para_ISKS_SIMPLE_2020_C2.rds")
isks_complex_nc <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/nc_Exome_para_ISKScomplex_2020_C2.rds")
isks_nc <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/nc_Exome_para_ISKS_2020_C2.rds")
##no_C3
cpx_ISKS_tas_OR_df_nc <- do.call("rbind.data.frame", cpx_OR_fisher(isks_tas_nc, 354, 3205, "TASvsMGRB"))
cpx_ISKS_simple_OR_df_nc <- do.call("rbind.data.frame", cpx_OR_fisher(isks_simple_nc, 370, 3205, "SIMPLEvsMGRB"))
cpx_ISKS_complex_OR_df_nc <- do.call("rbind.data.frame", cpx_OR_fisher(isks_complex_nc, 832, 3205, "COMPLEXvsMGRB"))
cpx_ISKS_OR_df_nc <- do.call("rbind.data.frame", cpx_OR_fisher(isks_nc, 1646, 3205, "ISKSvsMGRB"))
##combine results
nc_df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_tas_OR_df_nc, cpx_ISKS_simple_OR_df_nc, 
                                    cpx_ISKS_complex_OR_df_nc,cpx_ISKS_OR_df_nc)

##library(ggplot2)
forest_custplot_nc <- function(df){
  #  dotCOLS = c("#a6d8f0","#f9b282", "#78f542", "#f5c6e5", "#ffe3bf", "#cfabff")
  #  barCOLS = c("#008fd5","#de6b35", "#7a9406", "#fc0373", "#f79514", "#6c0ee8")
  ##without Cairns
  dotCOLS = c("#a6d8f0","#f9b282", "#78f542", "#f5c6e5")
  barCOLS = c("#008fd5","#de6b35", "#7a9406", "#fc0373")
  
 # p <- ggplot(df, aes(x=gene, y=log2(OR_Fish), ymin=log2(CI_lower), ymax=log2(CI_upper),col=Coh,fill=Coh)) + 
  p <- ggplot(df, aes(x=gene, y=OR_Fish, ymin=CI_lower, ymax=CI_upper,col=Coh,fill=Coh)) +
    #specify position here
    geom_linerange(size=1,position=position_dodge(width = 0.5)) +
    geom_hline(yintercept=1, lty=2) +
    #  geom_hline(yintercept=1, lty=2) +
    #specify position here too
    geom_point(size=5, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
    scale_fill_manual(values=barCOLS)+
    scale_color_manual(values=dotCOLS)+
    scale_x_discrete(name="Protein modules") + 
    scale_y_continuous(name="Odds ratio") +
    coord_flip() +
    theme_minimal() + theme(legend.position="bottom") + theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size=8))
  return(p)
}

p_nc <- forest_custplot_nc(nc_df_comb_rnd2) + ggtitle("Non_coding variants") + guides(fill=guide_legend(nrow=2,byrow=TRUE))

##non-coding obtained from regulatory elements

