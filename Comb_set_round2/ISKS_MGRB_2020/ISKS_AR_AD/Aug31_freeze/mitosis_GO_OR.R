.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

`%nin%` = Negate(`%in%`)

library(org.Hs.eg.db)
library(GO.db)

#mitotic G2 DNA damage checkpoint

Go_terms <- c("mitotic cytokinesis checkpoint", "mitotic DNA integrity checkpoint", "mitotic G1/S transition checkpoint",
              "mitotic G2/M transition checkpoint", "signal transduction involved in mitotic cell cycle checkpoint",
              "mitotic spindle checkpoint")

##HAUS complex

#A protein complex that localizes to interphase centrosomes and to mitotic spindle tubules and 
#regulates mitotic spindle assembly and centrosome integrity;

##cell cycle process
Go_terms_new <- c("attachment of mitotic spindle microtubules to kinetochore","mitotic centrosome separation",
                  "mitotic spindle assembly", "mitotic spindle midzone assembly", "mitotic spindle midzone assembly",
                  "centrosome duplication", "centrosome duplication", "centrosome separation", 
                  "establishment of mitotic spindle localization", "regulation of centrosome cycle",
                  "regulation of mitotic centrosome separation")

get_goterm_list <- function(go_terms){
go_genes <- list()
for(i in 1:length(go_terms)){
  g1 <- GOID( GOTERM[ Term(GOTERM) == go_terms[i]])
  allegs = get(g1, org.Hs.egGO2ALLEGS)
  go_genes[[i]] = unique(unlist(mget(allegs,org.Hs.egSYMBOL)))
  print(i)
  }
names(go_genes) <- go_terms
return(go_genes)
}

#go_geneset <- get_goterm_list(Go_terms)
go_geneset <- get_goterm_list(Go_terms_new)

##Odds ratio using fisher's exact test(use ppi_final_res... file for fisher's test)

#Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_para_pc1234_SKAT_Enriched_ISKS_2020_Aug31.rds")
#isks_mgrb_genes <- Exome_pc123_srt_SKAT_case_enr_nCH[[4]]
isks_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_para_tab_ISKS_MGRB_C345_Aug31.rds")
isks_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Exome_para_tab_ISKS_MGRB_minusC3_Aug31.rds")

##Epithelial data
epit_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/SKAT/Exome_para_MGRB_EPIT_2020.rds")
epit_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/EPIT_AR_AD/SKAT/Exome_para_MGRB_EPIT_2020_minusC3.rds")

##ISKS_complex vs MGRB
isks_complex_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Complex/Exome_para_ISKS_Complex_Aug31.rds")
isks_complex_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Complex/Exome_para_ISKS_Complex_minusC3_Aug31.rds")

##ISKS_TAS vs MGRB
isks_tas_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/TAS/Exome_para_ISKS_TAS_Aug31.rds")
isks_tas_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/TAS/Exome_para_ISKS_TAS_minusC3_Aug31.rds")
##ISKS_SIMPLE vs MGRB
isks_simple_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Simple/Exome_para_ISKS_SIMPLE_2020_Aug31.rds")
isks_simple_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Simple/Exome_para_ISKS_SIMPLE_2020_minusC3_Aug31.rds")

##Combined Sarcoma replication set
sarc_comb_mgrb_genes <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/Exome_para_CombSarc_Aug31.rds")
sarc_comb_mgrb_genes_noC3 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/Comb_TCGA_BERG_NoSarc/Exome_para_CombSarc_minusC3Aug31.rds")

########

cpx_OR_fisher <- function(ppi_res,case_coh_size, cont_coh_size, coh, go_input){
  ft_df <- list()
  for(i in 1:length(go_input)){
    ppi_res_tab <- ppi_res[ppi_res$gene %in% go_input[[i]],]
    # ppi_res_tab[,2] <- ifelse(ppi_res_tab[,2] == 0, 1, ppi_res_tab[,2])
    inp <- c(sum(ppi_res_tab[,1]), case_coh_size - sum(ppi_res_tab[,1]) , 
             sum(ppi_res_tab[,2]), cont_coh_size - sum(ppi_res_tab[,2]))
    sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
    colnames(sim_mat) <- c("case", "cont")
    rownames(sim_mat) <- c("hits", "no_hits")
    #ft <- fisher.test(sim_mat, alternative = "greater")
    #ft <- fisher.test(sim_mat)
    ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.99)
    ft_df[[i]] <- cbind.data.frame("gene" = names(go_input)[i] ,"Cases" = sum(ppi_res_tab[,1]),
                                   "Controls" = sum(ppi_res_tab[,2]),
                                   "Fish_pval" = ft$p.value,"CI_lower" = ft$conf.int[1],
                                   "CI_upper" = ft$conf.int[2],
                                   "OR_Fish" = ft$estimate, "Coh" = coh)
  }
  return(ft_df)
}

cpx_ISKS_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_mgrb_genes, 1644, 3205, "ISKSvsMGRB", go_geneset))
cpx_ISKS_complex_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_complex_mgrb_genes, 837, 3205, "COMPLEXvsMGRB", go_geneset))
cpx_ISKS_tas_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_tas_mgrb_genes, 362, 3205, "TASvsMGRB", go_geneset))
cpx_ISKS_simple_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(isks_simple_mgrb_genes, 371, 3205, "SIMPLEvsMGRB", go_geneset))
cpx_EPIT_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(epit_mgrb_genes, 842, 3205, "EPITvsMGRB", go_geneset))
#cpx_CAIRN_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(cairns_mgrb_genes, 413, 3205, "CAIRNvsMGRB", go_geneset))
cpx_REPSET_OR_df <- do.call("rbind.data.frame", cpx_OR_fisher(sarc_comb_mgrb_genes, 563, 3205, "REPvsMGRB", go_geneset))
df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_OR_df, cpx_ISKS_complex_OR_df, cpx_ISKS_tas_OR_df, 
                                 cpx_ISKS_simple_OR_df, cpx_EPIT_OR_df, cpx_REPSET_OR_df)
df_comb_rnd2$Coh <- as.factor(gsub("vsMGRB", "", df_comb_rnd2$Coh))
df_comb_rnd2_filt <- df_comb_rnd2[df_comb_rnd2$Fish_pval < 0.05, ]

##ISKS_MGRB, EPIT_MGRB and CAIRNS_MGRB without C3.
##no_C3
cpx_ISKS_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_mgrb_genes_noC3, 1644, 3205, "ISKSvsMGRB", go_geneset))
cpx_ISKS_complex_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_complex_mgrb_genes_noC3, 837, 3205, "COMPLEXvsMGRB", go_geneset))
cpx_ISKS_tas_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_tas_mgrb_genes_noC3, 362, 3205, "TASvsMGRB", go_geneset))
cpx_ISKS_simple_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(isks_simple_mgrb_genes_noC3, 371, 3205, "SIMPLEvsMGRB", go_geneset))
cpx_EPIT_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(epit_mgrb_genes_noC3, 842, 3205, "EPITvsMGRB", go_geneset))
cpx_REPSET_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(sarc_comb_mgrb_genes_noC3, 563, 3205, "REPvsMGRB", go_geneset))
#cpx_CAIRN_OR_df_noC3 <- do.call("rbind.data.frame", cpx_OR_fisher(Cairns_noC3_fin, 413, 3205, "CAIRNvsMGRB", go_geneset))
#noC3_df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_OR_df_noC3, cpx_ISKS_complex_OR_df_noC3, 
#                                      cpx_ISKS_tas_OR_df_noC3,cpx_ISKS_simple_OR_df_noC3,
#                                      cpx_EPIT_OR_df_noC3, cpx_CAIRN_OR_df_noC3)
##without Cairns
noC3_df_comb_rnd2 <- rbind.data.frame(cpx_ISKS_OR_df_noC3, cpx_ISKS_complex_OR_df_noC3, 
                                      cpx_ISKS_tas_OR_df_noC3,cpx_ISKS_simple_OR_df_noC3,
                                      cpx_EPIT_OR_df_noC3,cpx_REPSET_OR_df_noC3)
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
    dotCOLS = c("#a6d8f0","#f9b282", "#78f542", "#f5c6e5","#cfabff", "#232533")
    barCOLS = c("#008fd5","#de6b35", "#7a9406", "#fc0373","#6c0ee8", "#a3ceff")
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
#multiplot(p1, p2, cols = 2)
#multiplot(p3, p4, cols = 2)
#multiplot(p1, p3, cols = 2)
multiplot(p2, p4, cols = 2)

