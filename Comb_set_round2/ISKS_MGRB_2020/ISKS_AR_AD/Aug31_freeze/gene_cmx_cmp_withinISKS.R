.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
library(readxl)

`%nin%` = Negate(`%in%`)

##Complexes updated on Sep.27.2020

#Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "NBN", "SMARCAL1", "STAG3", "TIMELESS")

Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")

Glutamate_NMDA <- c("DLG1", "GRIN2A", "GRIA1", "CAM2KB", "CAM2KD")

ZBTB16_complex <- c("ZBTB16", "ATG7", "ASB10", "ASB8", "FBXO7")

#CEP_HAUS_complex <- c("CEP63", "CEP72", "HAUS4", "HAUS5")
CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
CEP_HAUS_extn1 <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1", "CEP57", "PCM1")
#CEP_HAUS_extn2 <- c("MACC1", "SMPD2", "CEP63", "CEP72", "CEP89", "HAUS4", "HAUS5", "MZT1", "PCM1")
MPNST_pos <- c("NF1", "LZTR1", "SDHA", "SDHB", "SDHD")
TP53_cont <- c("TP53") ##negative control for MPNST 

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

Sarcoma_genes <- c("TP53", "NF1", "SDHA", "SDHB", "SDHD")

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

cpx_genes <- unique(c(Shelterin, Glutamate_NMDA, ZBTB16_complex, CEP_HAUS_core, CEP_HAUS_extn1, MPNST_pos, TP53_cont, CEP_HAUS_SKATint_complex, 
                      CEP_HAUS_Biogrid_complex, CENP_complex, TCP1_complex,PRPF4_complex, Breast_cancer_genes, Sarcoma_genes, 
                      Centrosome_genes, mito_HSA_69278, mito_chk_point, mito_43_overlap, cent_rep_genes))

cpx_list <- list(Shelterin, Glutamate_NMDA, ZBTB16_complex, CEP_HAUS_core, CEP_HAUS_extn1, MPNST_pos, TP53_cont, CEP_HAUS_SKATint_complex, 
                 CEP_HAUS_Biogrid_complex, CENP_complex, TCP1_complex,PRPF4_complex, Breast_cancer_genes, Sarcoma_genes, 
                 Centrosome_genes, mito_HSA_69278, mito_chk_point, mito_43_overlap, cent_rep_genes)
names(cpx_list) <- c("Shelterin", "Glutamate_NMDA", "ZBTB16_complex", "CEP_HAUS_core", "CEP_HAUS_extn1", "MPNST_pos", "TP53_cont",
                     "CEP_HAUS_SKAT", "CEP_HAUS_Biogrid", "CENP_complex", "TCP1_complex", "PRPF4_complex", "BRCA_genes", 
                     "SARC_genes", "Centrosome_genes", "mito_HSA_69278", "mito_chk_GO", "mito_43", "cent_rep_genes")


###gene complex within ISKS sarcoma subtype enrichment
PID_file <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/PID_combset2020_CGC_skatBin_repstress_potint_mito_chkpt_centrosome_predNFE_clueGO_Sep192020_AD_addC4C5.tsv",
                       sep = "\t", header = T, stringsAsFactors = F)
##MPNST genes
#CEP_HAUS_extn2 <- c("MACC1", "SMPD2", "CEP63", "CEP72", "CEP89", "HAUS4", "HAUS5", "MZT1", "PCM1")

##GIST/MPNST cases
PID_file_GIST_MPNST <- PID_file[PID_file$sarcomatype %in% c("GIST", "MPNST"),]
##within cohort control
PID_file_control <- PID_file[PID_file$sarcomatype %nin% c("GIST", "MPNST"),]

##variant file(SKAT input)
fil_tab_noCH <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
fil_tab_isks <- fil_tab_noCH[grepl("ISKS_AR_AD", fil_tab_noCH$set),]

##Variant set hits in ISKS
var_sel <- function(var_type, gene_cpx){
  if(var_type == "inc_C3"){
fil_tab_isks_hits <- fil_tab_isks[fil_tab_isks$gene_symbol %in% gene_cpx,]
  }
  else if(var_type == "C4_C5"){
    fil_tab_isks_hits <- fil_tab_isks[fil_tab_isks$gene_symbol %in% gene_cpx &
                                        fil_tab_isks$auto_call %nin% "C3",]
  }
  return(fil_tab_isks_hits)
}

#isks_hits <- var_sel("inc_C3")
#isks_hits <- var_sel("C4_C5")
##Fisher's test
isks_cpx_fisher <- function(gene_cpx, var_type, cpx_name, pid_case, pid_control){
  inp_hit <- var_sel(var_type, gene_cpx)
  case_hits <- length(pid_case$pmn[pid_case$pmn %in% inp_hit$SAMPLE])
  case_tot <- length(pid_case$pmn)
  
  control_hits <- length(pid_control$pmn[pid_control$pmn %in% inp_hit$SAMPLE])
  control_tot <- length(pid_control$pmn)
  
  mpnst_inp <- c(case_hits, case_tot - case_hits , control_hits, control_tot - control_hits)
  sim_mat <- matrix(mpnst_inp ,nrow = 2, ncol = 2)
  colnames(sim_mat) <- c("case", "cont")
  rownames(sim_mat) <- c("hits", "no_hits")
  
  ft_mpnst <- fisher.test(sim_mat, conf.int = T, conf.level = 0.95)
  
  cpx_or_df <- cbind.data.frame("gene_cpx" =  cpx_name, "Cases" = case_hits,
                                "Controls" = control_hits,
                                "Fish_pval" = ft_mpnst$p.value,"CI_lower" = ft_mpnst$conf.int[1],
                                "CI_upper" = ft_mpnst$conf.int[2],
                                "OR_Fish" = ft_mpnst$estimate, "Coh" = "GIST_MPNSTvsOthers")
  return(cpx_or_df)
  
}

##
fish_df_list <- list()  
for(i in 1:length(cpx_list)){

 #fish_df_list[[i]] <- isks_cpx_fisher(cpx_list[[i]], "inc_C3", names(cpx_list)[i], 
 #                                     PID_file_GIST_MPNST, PID_file_control)
  fish_df_list[[i]] <- isks_cpx_fisher(cpx_list[[i]], "C4_C5", names(cpx_list)[i], 
                                       PID_file_GIST_MPNST, PID_file_control)

}

View(do.call("rbind.data.frame", fish_df_list))


##sarcoma genomicclass subtypes (complex versus others)
PID_file_complex <- PID_file[PID_file$genomicclass %in% c("Complex"),]
##within cohort control
PID_file_others <- PID_file[PID_file$genomicclass %nin% c("Complex"),]

fish_df_list <- list()  
for(i in 1:length(cpx_list)){
  
  #fish_df_list[[i]] <- isks_cpx_fisher(cpx_list[[i]], "inc_C3", names(cpx_list)[i], 
  #                                     PID_file_complex, PID_file_others)
  fish_df_list[[i]] <- isks_cpx_fisher(cpx_list[[i]], "C4_C5", names(cpx_list)[i], 
                                       PID_file_complex, PID_file_others)
  
}

View(do.call("rbind.data.frame", fish_df_list))
