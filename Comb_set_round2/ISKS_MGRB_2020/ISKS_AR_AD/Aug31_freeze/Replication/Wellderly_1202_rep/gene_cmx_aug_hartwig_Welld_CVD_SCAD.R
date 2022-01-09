##SCAD_88+CVD_117+Wellderly_1202 filtering and analysis

##gene_complexes
#.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
##New analysis with Wellderly data as control instead of MGRB

`%nin%` = Negate(`%in%`)

Glutamate_NMDA <- c("DLG1", "GRIN2A", "GRIA1", "CAM2KB", "CAM2KD")

ZBTB16_complex <- c("ZBTB16", "ATG7", "ASB10", "ASB8", "FBXO7")

CEP_HAUS_complex <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")

CENP_complex <- c("CENPC", "CENPF", "BUB1", "BUB1B", "AURKB", "RANGAP1", "RACGAP1", "MAD2L2")

TCP1_complex <- c("TCP1", "CCT2", "CCT6A")

PRPF4_complex <- c("PRPF4", "DHX9", "DHX16", "SF3A1", "CPSF3", "PPIL3")

#Sarcoma_genes <- c("TP53", "NF1", "BRCA2", "ERCC2", "SDHA", "SDHB", "SDHD")

#Breast_cancer_genes <- c("BRCA1", "BRCA2", "PALB2", "RAD51C", "CDH1")

##make consistent with Fig.2
Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")
Sarcoma_genes <- c("TP53", "SDHA", "SDHB", "SDHD")
Breast_cancer_genes <- c("BRCA1", "BRCA2", "PALB2")

#####

TP53 <- c("TP53")

cpx_genes <- unique(c(Shelterin, Glutamate_NMDA, ZBTB16_complex, CEP_HAUS_complex, CENP_complex, TCP1_complex,
                      PRPF4_complex, Breast_cancer_genes, Sarcoma_genes, TP53))

cpx_list <- list(Shelterin, Glutamate_NMDA, ZBTB16_complex, CEP_HAUS_complex, CENP_complex, TCP1_complex,
                 PRPF4_complex, Breast_cancer_genes, Sarcoma_genes, TP53)
names(cpx_list) <- c("Shelterin", "Glutamate_NMDA", "ZBTB16_complex", "CEP_HAUS", "CENP_complex", 
                     "TCP1_complex", "PRPF4_complex", "BRCA_genes", "SARC_genes", "TP53")


##Odds ratio using fisher's exact test(use ppi_final_res... file for fisher's test)


#isks_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/isks_mgrb_all.rds")
isks_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/isks_mgrb_all_noC3.rds")


##Wellderly variants
Wellder_var <- read.delim("~/RVAS/Wellderly_2012/Wellderly_vep_pois_filt4_C45_GQX40.tsv", sep = "\t",
                          header = T, stringsAsFactors = F)
Wellder_var_gene_tab <- as.data.frame(table(Wellder_var$SYMBOL))
colnames(Wellder_var_gene_tab) <- c("gene", "count")

##use Hart_var and sarcoma_comb variants from previous script gene_cmx_POIS_aug_hartwig.R
aug_all_hart_genes = readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/repset_563_Hartwig_POIS_noC3_genes_tab.rds")
aug_all_hart_genes$well_var = Wellder_var_gene_tab[match(aug_all_hart_genes$gene, Wellder_var_gene_tab$gene), 2]
aug_all_hart_genes$well_var = ifelse(is.na(aug_all_hart_genes$well_var), 0 , aug_all_hart_genes$well_var)
aug_all_hart_genes_welld <- aug_all_hart_genes[,c(1,5,3,4)]
saveRDS(aug_all_hart_genes_welld, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/repset_563_Hartwig_POIS_noC3_vs_Wellderly_genes_tab.rds", compress = T)
##Augment Wellderly with CVD + SCAD
CVD_var <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CVD_DIANE/round2/CVD_AR_AD/SKAT/Exome_para_pc123_SKAT_Enriched_CVD_2021_noC3.rds")

aug_all_hart_genes_welld$cvd_var <- CVD_var[match(aug_all_hart_genes_welld$gene, CVD_var$gene),1]
aug_all_hart_genes_welld$cvd_var <- ifelse(is.na(aug_all_hart_genes_welld$cvd_var), 0, aug_all_hart_genes_welld$cvd_var)

SCAD_var <- read.delim("~/RVAS/SCAD_88/SCAD88_C345_variants.tsv", sep = "\t", header = T,
                       stringsAsFactors = F)
SCAD_var_mod = SCAD_var[,c(11,29,80)]
SCAD_var_mod$SCAD = 1
SCAD_var_mod$MGRB = CVD_var[match(SCAD_var_mod$SYMBOL, CVD_var$gene), 2]
SCAD_var_mod = SCAD_var_mod[!is.na(SCAD_var_mod$MGRB),]
SCAD_var_mod_noC3 = SCAD_var_mod[SCAD_var_mod$auto_call %nin% "C3",]

##consolidate variants

aug_all_hart_genes_welld$scad_var <- SCAD_var_mod_noC3[match(aug_all_hart_genes_welld$gene, SCAD_var_mod_noC3$SYMBOL),4]
  
aug_all_hart_genes_welld$scad_var <-  ifelse(is.na(aug_all_hart_genes_welld$scad_var), 0, aug_all_hart_genes_welld$scad_var)

aug_all_hart_genes_welld$well_var <- aug_all_hart_genes_welld$well_var + aug_all_hart_genes_welld$cvd_var + aug_all_hart_genes_welld$scad_var

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
    ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.95)
    ft_df[[i]] <- cbind.data.frame("gene" = names(cpx_list)[i] ,"Cases" = sum(ppi_res_tab[,1]),
                                   "Controls" = sum(ppi_res_tab[,2]),
                                   "Fish_pval" = ft$p.value,"CI_lower" = ft$conf.int[1],
                                   "CI_upper" = ft$conf.int[2],
                                   "OR_Fish" = ft$estimate, 
                                   "case_coh_size" = case_coh_size, 
                                   "cont_coh_size" = cont_coh_size,
                                   "Coh" = coh)
  }
  return(ft_df)
}

cpx_ISKS_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_mgrb_genes_noC3, 1644, 3205, "ISKSvsMGRB"))
#cpx_repset_563 <- do.call("rbind.data.frame", cpx_OR_fisher(repset_563, 563, 3205, "rep563vsMGRB")) ##contains C3
##563(repset) + 276(Hartwig)
##Wellderly
#cpx_repset_839_POIS_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(aug_all_hart_genes_welld, 839, 1202, "rep839vsWellderly"))
##Wellderly(n=1202) + CVD(n=117) 
#cpx_repset_839_POIS_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(aug_all_hart_genes_welld, 839, 1319, "rep839vsWellderly"))

#Wellderly(n=1202) + CVD(n=117) + SCAD(n=88)

cpx_repset_839_POIS_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(aug_all_hart_genes_welld, 839, 1407, "rep839vsWellderly_CVD_SCAD"))


noC3_df_comb_hart <- rbind.data.frame(cpx_ISKS_OR_df_noC3, cpx_repset_839_POIS_noC3)
#noC3_df_comb_hart_sp <- rbind.data.frame(cpx_ISKS_OR_df_noC3, cpx_repset_839_POIS_noC3_sp)
##Plots
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
  # p <- ggplot(df, aes(x=gene, y=log2(OR_Fish), ymin=log2(CI_lower), ymax=log2(CI_upper),col=Coh,fill=Coh)) +
  p <- ggplot(df, aes(x=gene, y=OR_Fish, ymin=CI_lower, ymax=CI_upper,col=Coh,fill=Coh)) +
    #specify position here
    geom_linerange(size=1,position=position_dodge(width = 0.5)) +
    geom_hline(yintercept=0, lty=2) +
    #  geom_hline(yintercept=1, lty=2) +
    #specify position here too
    geom_point(size=5, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
    scale_fill_manual(values=barCOLS)+
    scale_color_manual(values=dotCOLS)+
    scale_x_discrete(name="Protein modules") + 
    scale_y_continuous(trans = "log", name="Log Odds ratio") +
    coord_flip() +
    theme_minimal() + theme(legend.position="bottom") + theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size=8))
  return(p)
}

noC3_df_comb_hart <- noC3_df_comb_hart[noC3_df_comb_hart$gene %in% c("BRCA_genes","Shelterin", "CEP_HAUS", "SARC_genes", "TP53"),]
p1 <- forest_custplot(noC3_df_comb_hart, repset = 1) + ggtitle("C4_C5") + guides(fill=guide_legend(nrow=2,byrow=TRUE))
# noC3_df_comb_hart_sp <- noC3_df_comb_hart_sp[noC3_df_comb_hart_sp$gene %in% c("BRCA_genes","Shelterin", "CEP_HAUS", "SARC_genes"),]
# p2 <-forest_custplot(noC3_df_comb_hart_sp, repset = 1) + ggtitle("C4_C5_splice_filtered") + guides(fill=guide_legend(nrow=2,byrow=TRUE))


#source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")
#multiplot(p1, p2, cols = 2)



################
##Forestplot for publication
library(forestplot)

#aug_all_hart_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/repset_563_Hartwig_POIS_noC3_genes_tab.rds")
cpx_ISKS_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_mgrb_genes_noC3, 1644, 3205, "ISKSvsMGRB"))
##563(repset) + 276(Hartwig)
#cpx_repset_839_POIS_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(aug_all_hart_genes_welld, 839, 1202, "rep839vsWellderly"))
##563(repset) + 276(Hartwig) ; 1202(Wellderly) + 117 (CVD)
#cpx_repset_839_POIS_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(aug_all_hart_genes_welld, 839, 1319, "rep839vsWellderly"))
#563(repset) + 276(Hartwig) ; Wellderly(n=1202) + CVD(n=117) + SCAD(n=88)
cpx_repset_839_POIS_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(aug_all_hart_genes_welld, 839, 1407, "rep839vsWellderly_CVD_SCAD"))


comb_hart_list <- list()
comb_hart_list$ISKS <-  cpx_ISKS_OR_df_noC3
comb_hart_list$Repset <-  cpx_repset_839_POIS_noC3
#comb_hart_list <- lapply(comb_hart_list, function(x)x[x[,1] %in% c("BRCA_genes","Shelterin", "CEP_HAUS", "SARC_genes", "TP53"),])
##remove TP53
comb_hart_list <- lapply(comb_hart_list, function(x)x[x[,1] %in% c("BRCA_genes","Shelterin", "CEP_HAUS", "SARC_genes"),])


comb_hart_list <- lapply(comb_hart_list, function(x) 
{
  xlog <- apply(x[,5:7], 2, log2)
  x[,4] <- formatC(x[,4], format = "e", digits = 2)
  x[,5:7] <- xlog
  x[,1] <- gsub("CEP_HAUS", "Centrosome", as.character(x[,1]))
  x[,1] <- gsub("BRCA_genes", "HBOC", as.character(x[,1]))
  rownames(x) <- as.character(x[,1])
  return(x)
})


# clrs <- fpColors(box = "royalblue",line = "darkblue", summary = "royalblue")
# forestplot(as.character(comb_hart_list$Repset[,1]), 
#            boxsize = 0.2,
#            rbind(comb_hart_list$Repset[,c(7,5,6)]),
#            clip = c(-3, 10),
#            col = clrs,
#            xlab = "log2 OR")


##Paired forestplot
library(forestplot)
forestplot::forestplot(as.character(comb_hart_list$Repset[,1]), 
                       legend = c("Discovery", "Replication"),
                       fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
                       txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex = 1),
                                        ticks = gpar(fontfamily = "", cex = 1),
                                        xlab  = gpar(fontfamily = "", cex = 1),
                                        legend = gpar(fontfamily = "", cex = 1)),
                       boxsize = .15, # We set the box size to better visualize the type
                       line.margin = .2, # We need to add this to avoid crowding
                       mean = cbind(comb_hart_list$ISKS[, "OR_Fish"], comb_hart_list$Repset[, "OR_Fish"]),
                       lower = cbind(comb_hart_list$ISKS[, "CI_lower"], comb_hart_list$Repset[, "CI_lower"]),
                       upper = cbind(comb_hart_list$ISKS[, "CI_upper"], comb_hart_list$Repset[, "CI_upper"]),
                       clip = c(-2, 10),
                       col = fpColors(box = c("blue", "darkred"), line = c("darkblue", "darkred")),
                       xlab = "log2 OR (95% CI)",
                       vertices = TRUE)


##add OR and p-values

comb_hart_list_all <- do.call("rbind.data.frame", comb_hart_list)
colnames(comb_hart_list_all)[7] <- "log2_OR"

#write.table(comb_hart_list_all, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/ISKS_repset_563_Hartwig_POIS_noC3_geneset_OR_synced_Wellderly_CVD_95CI.tsv",
#            row.names = F, quote = F, sep = "\t")
write.table(comb_hart_list_all, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/ISKS_repset_563_Hartwig_POIS_noC3_geneset_OR_synced_Wellderly_CVD_SCAD_95CI.tsv",
            row.names = F, quote = F, sep = "\t")
##done......

##Below are test runs
##Paired forestplot
# svg("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/augmented_sarc/Fig3A_synced_wellderly.svg", height = 8, width = 8)
# forestplot(as.character(comb_hart_list$Repset[,1]), 
#            legend = c("Discovery", "Replication"),
#            fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
#            txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex = 1),
#                             ticks = gpar(fontfamily = "", cex = 1),
#                             xlab  = gpar(fontfamily = "", cex = 1),
#                             legend = gpar(fontfamily = "", cex = 1)),
#            boxsize = .15, # We set the box size to better visualize the type
#            line.margin = .2, # We need to add this to avoid crowding
#            mean = cbind(comb_hart_list$ISKS[, "OR_Fish"], comb_hart_list$Repset[, "OR_Fish"]),
#            lower = cbind(comb_hart_list$ISKS[, "CI_lower"], comb_hart_list$Repset[, "CI_lower"]),
#            upper = cbind(comb_hart_list$ISKS[, "CI_upper"], comb_hart_list$Repset[, "CI_upper"]),
#            clip = c(-2, 10),
#            col = fpColors(box = c("blue", "darkred"), line = c("darkblue", "darkred")),
#            xlab = "log2 OR",
#            vertices = TRUE)
# dev.off()
#legend_args = fpLegend(pos = list(x = .85, y = 0.25), 
#gp = gpar(col = "#CCCCCC", fill = "#F9F9F9")),

library(grid)
grid.newpage()
borderWidth <- unit(4, "pt")
width <- unit(convertX(unit(1, "npc") - borderWidth, unitTo = "npc", valueOnly = TRUE)/2, "npc")
pushViewport(viewport(layout = grid.layout(nrow = 1, 
                                           ncol = 3, 
                                           widths = unit.c(width,
                                                           borderWidth,
                                                           width))
)
)
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 1))

clrs <- fpColors(box = "royalblue",line = "darkblue", summary = "royalblue")
forestplot(as.character(comb_hart_list$ISKS[,1]),
           title = "Discovery",
           boxsize = 0.2,
           rbind(comb_hart_list$ISKS[,c(7,5,6)]),
           clip = c(-3, 10),
           col = clrs,
           xlab = "log2 OR",
           new_page = FALSE)

upViewport()
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 2))
grid.rect(gp = gpar(fill = "#dddddd", col = "#eeeeee"))
upViewport()
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 3))

#clrs <- fpColors(box = "royalred",line = "darkred", summary = "royalred")
forestplot(as.character(comb_hart_list$Repset[,1]),
           title = "Replication",
           boxsize = 0.2,
           rbind(comb_hart_list$Repset[,c(7,5,6)]),
           clip = c(-3, 10),
           col = clrs,
           xlab = "log2 OR",
           new_page = FALSE)
upViewport(2)


##Alternate approach
#par(mfrow=c(2,1))
#for adding all values
#cbind(comb_hart_list$ISKS[1:4], round(comb_hart_list$ISKS[5],3),
#      round(comb_hart_list$ISKS[6],3), round(comb_hart_list$ISKS[7],3))
forestplot(comb_hart_list$Repset[1:4],
           boxsize = 0.125,
           txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex = 1),
                            ticks = gpar(fontfamily = "", cex = 1),
                            xlab  = gpar(fontfamily = "", cex = 1),
                            legend = gpar(fontfamily = "", cex = 1)),
           rbind(comb_hart_list$ISKS[,c(7,5,6)]),
           new_page = TRUE,
           clip = c(-3, 10), 
           col = fpColors(box = "darkred",
                          line = "darkred",
                          summary = "darkred",
                          hrz_lines = "#444444"),
           xlab = "log2 OR",
           vertices = TRUE)

forestplot(comb_hart_list$Repset[1:4], 
           boxsize = 0.125,
           txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex = 1),
                            ticks = gpar(fontfamily = "", cex = 1),
                            xlab  = gpar(fontfamily = "", cex = 1),
                            legend = gpar(fontfamily = "", cex = 1)),
           rbind(comb_hart_list$Repset[,c(7,5,6)]),new_page = TRUE,
           clip = c(-3, 10), 
           col = fpColors(box = "royalblue",
                          line = "darkblue",
                          summary = "royalblue",
                          hrz_lines = "#444444"),
           xlab = "log2 OR",
           vertices = TRUE)

# tab_rep0 <- comb_hart_list$Repset[2:4]
# #colnames(tab_rep0) <- NULL
# tab_rep <- list(c(NA, rownames(comb_hart_list$Repset)), tab_rep0)
# forestplot(tab_rep, 
#            rbind(comb_hart_list$Repset[,c(7,5,6)]),new_page = TRUE,
#            clip = c(-3, 10), 
#            col = fpColors(box = "royalblue",
#                           line = "darkblue",
#                           summary = "royalblue",
#                           hrz_lines = "#444444"),
#            vertices = TRUE)

#####
# library(tidyverse)
# library(ggforestplot)
# df_linear <-
#   ggforestplot::df_linear_associations %>%
#   dplyr::arrange(name) %>%
#   dplyr::filter(dplyr::row_number() <= 30)
# 
# noC3_df_comb_hart_fp <- noC3_df_comb_hart
# noC3_df_comb_hart_fp$se = (noC3_df_comb_hart_fp$CI_upper - noC3_df_comb_hart_fp$CI_lower)/3.92 ##95% CI
# colnames(noC3_df_comb_hart_fp)[4] = "pvalue"
# noC3_df_comb_hart_fp$gene <- as.character(noC3_df_comb_hart_fp$gene)
# noC3_df_comb_hart_fp$Coh <- as.character(noC3_df_comb_hart_fp$Coh)
# df_hart_rep <- noC3_df_comb_hart_fp %>%
#   dplyr::arrange(gene) 
# df_hart_rep1 <- tibble(df_hart_rep)
# ggforestplot::forestplot(
#   df = df_hart_rep1,
#   name = gene,
#   estimate = OR_Fish,
#   se = se,
#   logodds = TRUE,
#   colour = Coh,
#   psignif = 1,
#   ci = 0.95,
#   title = "Validation of functional modules",
#   xlab = "logOdds ratio"
# )