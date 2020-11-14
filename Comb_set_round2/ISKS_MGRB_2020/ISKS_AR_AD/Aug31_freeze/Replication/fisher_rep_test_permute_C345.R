##Validation test
##use only modules enriched using PPI analysis
##generate 1000 permutations of case control splits for training(2/3) and test sets(2/3)
##estimate ORs for each module in training and test sets
##check frequency of significant ORs for all modules across training and test set
##permutation test for modules; sample data agnostic of label in the 2/3:1/3 split a 1000 times for null
##and derive p-values from the average frequency based test statistic derived from  case control splits.

.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)
fil_tab_noCH <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
Ex_samp_id <- unique(fil_tab_noCH$SAMPLE)
#binary phenotype vector 
p_vec_all <- ifelse(!is.na(as.numeric(as.character(Ex_samp_id))) | grepl("^CR|^LK",as.character(Ex_samp_id)), 1, 0)
#set.seed(3453466)


##Complexes updated on Sep.27.2020

Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")

Glutamate_NMDA <- c("DLG1", "GRIN2A", "GRIA1", "CAM2KB", "CAM2KD")

ZBTB16_complex <- c("ZBTB16", "ATG7", "ASB10", "ASB8", "FBXO7")

CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
CEP_HAUS_extn1 <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1", "CEP57", "PCM1")

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

Sarcoma_genes <- c("TP53", "SDHA", "SDHB", "SDHD")

TP53_control <- c("TP53") #positive control

MPNST_pos <- c("NF1", "LZTR1", "SDHA", "SDHB", "SDHD")

Breast_cancer_genes <- c("BRCA1", "BRCA2", "PALB2")
#Breast_cancer_genes <- c("BRCA1", "BRCA2", "PALB2", "RAD51C", "RAD51D")

library(readxl)
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

fil_tab_cpx_genes <- fil_tab_noCH[fil_tab_noCH$gene_symbol %in% cpx_genes,]


#####make genotype matrix function
make_geno_mat <- function(ftemp_file, p_vec, samp_id){
  maf_thresh <- 3.5/(2*length(p_vec))
  ftemp_file <- ftemp_file[ftemp_file$VAF >= 0.35 & ftemp_file$comb_score >= 5.6,]
  ftemp_tab_var_id <- unique(ftemp_file$VARIANT)
  samp_vec <- list()
  for(m in 1:length(ftemp_tab_var_id)){
    sam_gene_gt <- ftemp_file[ftemp_file$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,82,127:128)]
    sam_gene_gt <- unique(sam_gene_gt)
    # print(dim(sam_gene_gt)[1])
    if(dim(sam_gene_gt)[1] > 1 & length(unique(sam_gene_gt$vep_consequence)) > 1){
      # if(dim(sam_gene_gt)[1] > 1 & length(unique(sam_gene_gt$vep_consequence)) > 1 & length(unique(sam_gene_gt$comb_score)) > 1){
      vep_con <- unique(sam_gene_gt$vep_consequence)
      samp_vec_con <- list()
      for(k in 1:length(vep_con)){
        sam_gene_gt_con <- sam_gene_gt[sam_gene_gt$vep_consequence %in% vep_con[k],]
        sam_gene_gt_con <- unique(sam_gene_gt_con)
        cont_wt <- sum(sam_gene_gt_con[grepl("^[ABZ]",sam_gene_gt_con$SAMPLE),]$comb_score)
        case_wt <- sum(sam_gene_gt_con[!grepl("^[ABZ]",sam_gene_gt_con$SAMPLE),]$comb_score)
        maf_vec_cont <- sum(dim(sam_gene_gt_con[grepl("^[ABZ]",sam_gene_gt_con$SAMPLE),])[1])/(2*length(p_vec))
        maf_vec_case <- sum(dim(sam_gene_gt_con[!grepl("^[ABZ]",sam_gene_gt_con$SAMPLE),])[1])/(2*length(p_vec))
        
        ##genotype matrix  
        samp_id_set <- samp_id[p_vec]
        sam_gene_gt_con$add_mod <- as.numeric(sam_gene_gt_con$GT)
        sam10 <- ifelse(samp_id_set %in% sam_gene_gt_con$SAMPLE, 1, 0)
        sam10[which(sam10 != 0)] <- sam_gene_gt_con$add_mod ##additive model
        samp_vec_con[[k]] <- c(unique(sam_gene_gt_con$VARIANT),
                               unique(sam_gene_gt_con$gene_symbol), 
                               unique(sam_gene_gt_con$vep_consequence), 
                               unique(as.character(sam_gene_gt_con$auto_call)),
                               as.numeric(unique(sam_gene_gt_con$comb_score)), 
                               as.numeric(maf_vec_case),
                               as.numeric(maf_vec_cont), 
                               as.numeric(case_wt),
                               as.numeric(cont_wt),
                               sam10)
      }
      
    }
    else{
      ##compute cohort specific MAF
      cont_wt <- sum(sam_gene_gt[grepl("^[ABZ]",sam_gene_gt$SAMPLE),]$comb_score)
      case_wt <- sum(sam_gene_gt[!grepl("^[ABZ]",sam_gene_gt$SAMPLE),]$comb_score)
     # maf_vec_cont <- round(sum(ifelse(is.na(as.numeric(sam_gene_gt$SAMPLE)), 1, 0))/(2*length(p_vec)), 10)
      maf_vec_cont <- round(length(grep("^[ABZ]", sam_gene_gt$SAMPLE))/(2*length(p_vec)), 10)
      maf_vec_case <- round(sum(ifelse(!is.na(as.numeric(sam_gene_gt$SAMPLE)) | 
                                         grepl("^CR|^LK", as.character(sam_gene_gt$SAMPLE)), 1, 0))/(2*length(p_vec)),10)
      #maf_vec <- (maf_vec_cont + maf_vec_case)/(2*(1572 + 1110))
      ##genotype matrix  
      samp_id_set <- samp_id[p_vec]
      sam_gene_gt$add_mod <- as.numeric(sam_gene_gt$GT)
      sam10 <- ifelse(samp_id_set %in% sam_gene_gt$SAMPLE, 1, 0)
      sam10[which(sam10 != 0)] <- sam_gene_gt$add_mod ##additive model
      ##account for discrepancy in vep_consequence
      #vep_con <- names(sort(table(unlist(strsplit(sam_gene_gt$vep_consequence, split = "&"))),decreasing=TRUE)[1])
      # sam_gene_gt$vep_consequence <- vep_con
      
      
      samp_vec[[m]] <- c(ftemp_tab_var_id[m],
                         unique(sam_gene_gt$gene_symbol), 
                         unique(sam_gene_gt$vep_consequence), 
                         unique(as.character(sam_gene_gt$auto_call)),
                         as.numeric(unique(sam_gene_gt$comb_score)), 
                         as.numeric(maf_vec_case),
                         as.numeric(maf_vec_cont), 
                         as.numeric(case_wt),
                         as.numeric(cont_wt),
                         sam10)
    }
  }
  samp_vec_mat_uni <- do.call("rbind.data.frame", samp_vec)
  colnames(samp_vec_mat_uni) <- c("VARIANT", "gene_symbol", "vep_consequence", "auto_call", "comb_score", "coh_MAF_case",
                                  "coh_MAF_cont", "case_wt", "cont_wt", samp_id_set)
  if(exists("samp_vec_con")){
    samp_vec_mat_con <- do.call("rbind.data.frame", samp_vec_con)
    colnames(samp_vec_mat_con) <- c("VARIANT", "gene_symbol", "vep_consequence", "auto_call", "comb_score", "coh_MAF_case",
                                    "coh_MAF_cont", "case_wt", "cont_wt", samp_id_set)
    samp_vec_mat <- rbind.data.frame(samp_vec_mat_uni, samp_vec_mat_con)
  }else{
    samp_vec_mat <- samp_vec_mat_uni
  }
  
  ##Intracohort filter  : equivalent to ~ 5/1661*2 for case ; ~ 5/3209*2 for control
  ##changed to uniform intra_cohort_MAF filter
  ## change to MAF filter to 3/1661*2
  #     samp_vec_mat <- samp_vec_mat[as.numeric(as.character(samp_vec_mat$coh_MAF_case)) <= 0.0015 | 
  #                                    as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) <= 0.001,]
  samp_vec_mat <- samp_vec_mat[!(as.numeric(as.character(samp_vec_mat$coh_MAF_case)) > maf_thresh | 
                                   as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) > maf_thresh),]
  print(dim(samp_vec_mat))
  
  return(samp_vec_mat)
  
}

###Fisher's test function
##Included C3's and add null frequency
fish_test_fun <- function(df_mat_inp, gene_sym, coh, set, null_mat_inp){
  #df_mat_inp <- df_mat_inp[df_mat_inp$auto_call %nin% "C3",]
  # test1_inf <- test2[,c(1:9)]
  df_mat <- as.matrix(df_mat_inp[,-c(1:9)])
  null_mat <- as.matrix(null_mat_inp[,-c(1:9)])
  class(df_mat) <- "numeric"
  class(null_mat) <- "numeric"
  cont_test <- sum(colSums(df_mat[,grepl("^[ABZ]", colnames(df_mat))]))
  cont_test <- ifelse(cont_test == 0, 1, cont_test)
  case_test <- sum(colSums(df_mat[,!grepl("^[ABZ]", colnames(df_mat))]))
  cont_tot <- dim(df_mat[,grepl("^[ABZ]", colnames(df_mat))])[2]
  case_tot <- dim(df_mat[,!grepl("^[ABZ]", colnames(df_mat))])[2]
  null_count <- sum(colSums(null_mat))
  inp <- c(case_test, case_tot - case_test, 
           cont_test, cont_tot - cont_test)
  sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
  colnames(sim_mat) <- c("case", "cont")
  rownames(sim_mat) <- c("hits", "no_hits")
  #ft <- fisher.test(sim_mat, alternative = "greater")
  #ft <- fisher.test(sim_mat)
  ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.95)
  cbind.data.frame("gene" = gene_sym ,"Cases" = case_test, "Cases_comp" = inp[2],
                   "Controls" = cont_test, "Controls_comp" = inp[4], "null_freq" = null_count,
                   "Fish_pval" = ft$p.value,"CI_lower" = ft$conf.int[1],
                   "CI_upper" = ft$conf.int[2],
                   "OR_Fish" = ft$estimate, "Coh" = coh, "Set" = set)
  
}

###Function for parallel processing
library(plyr)
library(doParallel)
registerDoParallel(25)
##bootstrap function

val_test = function(data, phen_vec, samp_id, gene_sym, coh, seed = 46574689, B = 1000, ...)
{
  ft_df <- list()
  test_set_fish = llply(1:B, function(i)
  {
    perm_data = data
    set.seed(seed + i)
    
    p_vec_case_ind <- which(phen_vec == 1)
    p_vec_cont_ind <- which(phen_vec == 0)
    p_vec_ind_case_train <- sample(p_vec_case_ind, size = round(length(p_vec_case_ind)*2/3), replace = F)
    p_vec_ind_case_test <- p_vec_case_ind[p_vec_case_ind %nin% p_vec_ind_case_train] #1/3 split
    p_vec_ind_case_null_disc <- sample(1:length(phen_vec), size = round(length(p_vec_case_ind)*2/3), replace = F)
    #table(phen_vec[p_vec_ind_case_null_disc])
    
    p_vec_cont_ind_train <- sample(p_vec_cont_ind, size = round(length(p_vec_cont_ind)*2/3), replace = F)
    p_vec_cont_ind_test <- p_vec_cont_ind[p_vec_cont_ind %nin% p_vec_cont_ind_train] #1/3 split
    p_vec_ind_case_null_rep <- sample(1:length(phen_vec), size = round(length(p_vec_case_ind)*1/3), replace = F)
    
    train_ind_set <- c(p_vec_ind_case_train, p_vec_cont_ind_train)
    test_ind_set <- c(p_vec_ind_case_test, p_vec_cont_ind_test)
    test1 <- make_geno_mat(ftemp_file = data, p_vec = train_ind_set, samp_id = samp_id)
    test2 <- make_geno_mat(ftemp_file = data, p_vec = test_ind_set, samp_id = samp_id)
    null_case_disc <- make_geno_mat(ftemp_file = data, p_vec = p_vec_ind_case_null_disc, samp_id = samp_id)
    null_case_rep <- make_geno_mat(ftemp_file = data, p_vec = p_vec_ind_case_null_rep, samp_id = samp_id)
    
    list(fish_test_fun(test1, coh = coh, gene_sym = gene_sym, set = "disc", null_mat_inp = null_case_disc), 
         fish_test_fun(test2, coh = coh, gene_sym = gene_sym, set = "rep", null_mat_inp = null_case_rep))
    
  }, .parallel = TRUE)
}

# t1 <- val_test(data = ftemp_shel, phen_vec = p_vec, samp_id = Ex_samp_id)
# t2 <- val_test(data = ftemp_shel, phen_vec = p_vec, samp_id = Ex_samp_id)
for(k in 1:length(cpx_list)){
  ftemp_cpx <- fil_tab_cpx_genes[fil_tab_cpx_genes$gene_symbol %in% cpx_list[[k]],]
  t1000 <- val_test(data = ftemp_cpx, phen_vec = p_vec_all, samp_id = Ex_samp_id, gene_sym = names(cpx_list)[k],coh = "ISKSvsMGRB")
  disc_set <- do.call("rbind.data.frame", lapply(t1000, function(x)x[[1]]))
  rep_set <- do.call("rbind.data.frame", lapply(t1000, function(x)x[[2]]))
  
  write.table(disc_set, paste0("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Validation_new_cpx/",
                               names(cpx_list)[k],"_C345_fisher.disc_null.tsv"), sep = "\t", quote = F, row.names = F)
  write.table(rep_set, paste0("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Validation_new_cpx/",
                              names(cpx_list)[k],"_C345_fisher.rep_null.tsv"), sep = "\t", quote = F, row.names = F)
}


##modifications for metanalysis of ORs
##Mantel-Haenszel method (more Robust)
# CI_meta_MH_data <- function(cmpx_name){
#   library(meta)
#   library(metafor)
#   df1_disc <- read.delim(paste0("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Validation_new_cpx/",cmpx_name,"_C345_fisher.disc_null.tsv"),
#                          sep = "\t", header = T, stringsAsFactors = F )
#   df1_rep <- read.delim(paste0("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Validation_new_cpx/",cmpx_name,"_C345_fisher.rep_null.tsv"),
#                         sep = "\t", header = T, stringsAsFactors = F )
#   comb_Mant_fun <- function(df){
#     fixed_mod <- metabin(Cases, Cases_comp, Controls, Controls_comp, data = df,  
#                          studlab = gene,
#                          method.tau = "SJ",
#                          comb.fixed = FALSE,
#                          comb.random = TRUE,
#                          hakn = TRUE,
#                          prediction = TRUE,
#                          incr = 0.1,
#                          sm = "OR")
#     summ <- summary(fixed_mod)
#     
#     #   ret_df <- cbind.data.frame("cmpx_name" = cmpx_name, "Fish_OR_meta" = exp(summ$fixed$TE), 
#     #                              "lower.CI" = exp(summ$fixed$lower), "upper.CI" = exp(summ$fixed$upper))
#     ##use prediction interval for the forest plot
#     #https://www.erim.eur.nl/research-support/meta-essentials/interpret-results/the-forest-plot/prediction-interval/
#     #https://bmjopen.bmj.com/content/6/7/e010247
#     ret_df <- cbind.data.frame("cmpx_name" = cmpx_name, "Fish_OR_meta" = exp(summ$random$TE), 
#                                "lower.CI" = exp(summ$predict$lower), "upper.CI" = exp(summ$predict$upper))
#     
#   }
#   
#   
#   return(list(comb_Mant_fun(df1_disc), comb_Mant_fun(df1_rep)))
#   
# }
# 
# CI_MH_meta <- list()
# for(i in 1:length(cpx_list)){
#   CI_MH_meta[[i]] <- CI_meta_MH_data(names(cpx_list)[i])
# }
# 
# library(ggplot2)
# CI_disc_MH_meta <- lapply(CI_MH_meta, function(x)x[[1]])
# CI_disc_MH_meta_df <- do.call("rbind.data.frame", CI_disc_MH_meta)
# CI_disc_MH_meta_df$set <- "Disc"
# CI_rep_MH_meta <- lapply(CI_MH_meta, function(x)x[[2]])
# CI_rep_MH_meta_df <- do.call("rbind.data.frame", CI_rep_MH_meta)
# CI_rep_MH_meta_df$set <- "Rep"
# 
# CI_comb_MH_meta_df <- rbind.data.frame(CI_disc_MH_meta_df, CI_rep_MH_meta_df)
# CI_comb_MH_meta_df$Index <- rownames(CI_comb_MH_meta_df)
# 
# xname <- expression(paste(italic("Log"), "OR"))
# #setting up the basic plot
# p_meta <- ggplot(data=CI_comb_MH_meta_df, aes(y=cmpx_name, x=log(Fish_OR_meta), xmin=log(lower.CI), xmax=log(upper.CI))) + 
#   geom_point() + geom_point(data=subset(CI_comb_MH_meta_df, set=="Disc"), color="Black", size=2) + 
#   geom_errorbarh(height=.1) + scale_x_continuous(limits=c(-2,5), breaks = c(-2:5), name=xname) +
#   geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) + 
#   facet_grid(set~., scales= "free", space="free") + ggtitle("PPI cliques")+
#   theme_minimal() + theme(text=element_text(family="Times",size=18, color="black"))+
#   theme(panel.spacing = unit(1, "lines"))
# 
# p_meta
# 
# # ##For supp figure
# CI_comb_MH_meta_df_sub <- CI_comb_MH_meta_df[CI_comb_MH_meta_df$cmpx_name %in% c("SARC_genes", "BRCA_genes",
#                                                                                  "TP53_control", "Shelterin",
#                                                                                  "CEP_HAUS_core", "MPNST_pos"),]
# 
# p_sub <- ggplot(data=CI_comb_MH_meta_df_sub, aes(y=cmpx_name, x=log(Fish_OR_meta), xmin=log(lower.CI), xmax=log(upper.CI))) +
#   geom_point() + geom_point(data=subset(CI_comb_MH_meta_df_sub, set=="Disc"), color="Black", size=2) +
#   geom_errorbarh(height=.1) + scale_x_continuous(limits=c(-2,4), breaks = c(-2:4), name=xname) +
#   geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) +
#   facet_grid(set~., scales= "free", space="free") + ggtitle("PPI cliques")+
#   theme_minimal() + theme(text=element_text(family="Times",size=18, color="black"))+
#   theme(panel.spacing = unit(1, "lines"))
# 
# p_sub

