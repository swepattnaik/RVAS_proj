##Validation test
##use only modules enriched using PPI analysis
##generate 1000 permutations of case control splits for training(2/3) and test sets(2/3)
##estimate ORs for each module in training and test sets
##check frequency of significant ORs for all modules across training and test set
##permutation test for modules; sample data agnostic of label in the 2/3:1/3 split a 1000 times for null
##and derive p-values from the average frequency based test statistic derived from  case control splits.

.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
fil_tab_noCH <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT_disc/all_DISC_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
Ex_samp_id <- unique(fil_tab_noCH$SAMPLE)
#binary phenotype vector 
p_vec_all <- ifelse(!is.na(as.numeric(as.character(Ex_samp_id))) | grepl("^CR|^LK",as.character(Ex_samp_id)), 1, 0)
#set.seed(3453466)
`%nin%` = Negate(`%in%`)

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
      maf_vec_cont <- round(sum(ifelse(is.na(as.numeric(sam_gene_gt$SAMPLE)), 1, 0))/(2*length(p_vec)), 10)
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
##Excluded C3's and add null frequency
fish_test_fun <- function(df_mat_inp, gene_sym, coh, set, null_mat_inp){
  df_mat_inp <- df_mat_inp[df_mat_inp$auto_call %nin% "C3",]
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
  cbind.data.frame("gene" = gene_sym ,"Cases" = case_test,
                   "Controls" = cont_test, "null_freq" = null_count,
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
  
  write.table(disc_set, paste0("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31_disc/Validation/",
                               names(cpx_list)[k],"_fisher.disc_null.tsv"), sep = "\t", quote = F, row.names = F)
  write.table(rep_set, paste0("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31_disc/Validation/",
                              names(cpx_list)[k],"_fisher.rep_null.tsv"), sep = "\t", quote = F, row.names = F)
}

##plots
# library(ggplot2)
# rep_disc_dens <- list()
# rep_disc_dens_filt <- list()
# rep_disc_hist <- list()
# rep_disc_hist_filt <- list()
# for(i in 1:length(cpx_list)){
# df1_disc <- read.delim(paste0("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Validation/",names(cpx_list)[i],"_fisher.disc.tsv"),
#                   sep = "\t", header = T, stringsAsFactors = F )
# df1_rep <- read.delim(paste0("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Validation/",names(cpx_list)[i],"_fisher.rep.tsv"),
#                              sep = "\t", header = T, stringsAsFactors = F )
# #df1_disc$B_ind <- as.integer(as.character(rownames(df1_disc)))
# #df1_rep$B_ind <- as.integer(as.character(rownames(df1_rep)))
# df1_comb <- rbind.data.frame(df1_disc, df1_rep)
# df1_comb_filt <- df1_comb[df1_comb$Fish_pval < 0.05,]
# # cor.test(df1_disc$OR_Fish, df1_rep$OR_Fish, method = "spearman")
# # cor.test(df1_disc$OR_Fish, df1_rep$OR_Fish, method = "pearson")
# # cor.test(df1_disc$OR_Fish, df1_rep$OR_Fish, method = "kendall")
# 
# rep_disc_hist[[i]] <- ggplot(df1_comb, aes(x = OR_Fish, fill = Set)) + geom_histogram(alpha=0.5, position="identity") + ggtitle(names(cpx_list)[i])
# rep_disc_hist_filt[[i]] <- ggplot(df1_comb_filt, aes(x = OR_Fish, fill = Set)) + geom_histogram(alpha=0.5, position="identity") + ggtitle(names(cpx_list)[i])
# #rep_disc_dens[[i]] <- ggplot(df1_comb, aes(x = OR_Fish, fill = Set)) + geom_density(alpha = 0.5) + ggtitle(names(cpx_list)[i])
# #rep_disc_dens_filt[[i]] <- ggplot(df1_comb_filt, aes(x = OR_Fish, fill = Set)) + geom_density(alpha = 0.5) + ggtitle(names(cpx_list)[i])
# }
# 
# 
# 
# 
# 
# source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")
# multiplot(plotlist = rep_disc_hist, cols = 3)
# multiplot(plotlist = rep_disc_hist_filt, cols = 3)
# # multiplot(plotlist = rep_disc_dens, cols = 3)
# # multiplot(plotlist = rep_disc_dens_filt, cols = 3)
# # 
# # df1_disc <- read.delim(paste0("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Validation/",names(cpx_list)[i],"_fisher.disc.tsv"), 
# #                        sep = "\t", header = T, stringsAsFactors = F )
# # df1_rep <- read.delim(paste0("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/Validation/",names(cpx_list)[i],"_fisher.rep.tsv"), 
# #                       sep = "\t", header = T, stringsAsFactors = F )
# # df1_disc$B_ind <- as.integer(as.character(rownames(df1_disc)))
# # df1_rep$B_ind <- as.integer(as.character(rownames(df1_rep)))
# # df1_comb <- rbind.data.frame(df1_disc, df1_rep)
# # 
# # ggplot(df1_comb, aes(x=Fish_pval, fill=Set)) + geom_histogram(alpha=0.5, position="identity")
# # ggplot(df1_comb, aes(x=OR_Fish, fill=Set)) + geom_histogram(alpha=0.5, position="identity")
# # ggplot(df1_comb, aes(x=Fish_pval, fill=Set)) + geom_density(alpha=0.5, position="identity")
# # 
# # 
# # df1_disc_filt <- df1_disc[df1_disc$Fish_pval < 0.001,]
# # 
# # disc_set_filt <- disc_set[disc_set$Fish_pval < 0.001,]
# # ggplot(data.frame(OR = disc_set_filt$OR_Fish),aes(x=OR)) +
# #   geom_histogram(binwidth=0.25,aes(y=..density..)) +
# #   geom_density(color="red")
# # 
# # df1_rep_filt <- df1_rep[df1_rep$Fish_pval < 0.001,]
# # 
# # rep_set_filt <- rep_set[rep_set$Fish_pval < 0.005,]
# # ggplot(data.frame(OR = rep_set_filt$OR_Fish),aes(x=OR)) +
# #   geom_histogram(binwidth=0.25,aes(y=..density..)) +
# #   geom_density(color="green")
# # 
# # 
# # library(boot)
# # fc_or <- function(d, i){
# #   d2 <- d[i,]
# #   return(d2$OR_Fish)
# # }
# # 
# # set.seed(342113)
# # bootOR_disc <- boot(disc_set, fc_or, R=500)
# # plot(bootOR_disc)
# # boot.ci(boot.out = bootOR_disc, type = c("norm", "basic", "perc"))
# # bootOR_rep <- boot(rep_set, fc_or, R=500)
# # plot(bootOR_rep)
# # boot.ci(boot.out = bootOR_rep, type = c("norm", "basic", "perc"))
# 
# 

######Use this
# ##split data cumulative distribution for both disc and rep set
# library(ggpubr)
# library(reshape2)
# 
# plot_val_data <- function(cmpx_name,plot_type = c("box", "dens")){
#   dens_plot <- list()
#   box_plots <- list()
#   df1_disc <- read.delim(paste0("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Validation/",cmpx_name,"_fisher.disc_null.tsv"),
#                          sep = "\t", header = T, stringsAsFactors = F )
#   df1_rep <- read.delim(paste0("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Validation/",cmpx_name,"_fisher.rep_null.tsv"),
#                         sep = "\t", header = T, stringsAsFactors = F )
# 
#   disc_melt <- melt(df1_disc[,2:4])
#   rep_melt <- melt(df1_rep[,2:4])
# 
#   if(plot_type == "box"){
#     p <- ggboxplot(disc_melt, x = "variable", y = "value",
#                    color = "variable", palette = "jco",
#                    add = "jitter")
#     p <- p + stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less")) +
#       ggtitle(paste0(cmpx_name, "_disc"))
# 
#     p1 <- ggboxplot(rep_melt, x = "variable", y = "value",
#                     color = "variable", palette = "jco",
#                     add = "jitter")
#     p1 <- p1 + stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less")) +
#       ggtitle(paste0(cmpx_name, "_rep"))
# 
#     box_plots <- list(p, p1)
#     return(box_plots)
#   }
# 
#   else if(plot_type == "dens"){
#     p <- ggplot(disc_melt, aes(x=value, fill=variable)) + geom_density(alpha=0.5, position="identity") +
#       ggtitle(cmpx_name)
# 
#     p1 <- ggplot(rep_melt, aes(x=value, fill=variable)) + geom_density(alpha=0.5, position="identity") +
#       ggtitle(cmpx_name)
#     dens_plots <- list(p, p1)
#     return(dens_plots)
#   }
# 
# }
# 
# 
# ##multiplot
# source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")
# 
# box_plot_all <- list()
# for(i in 1:length(cpx_list)){
#   box_plot_all[[i]] <- plot_val_data(names(cpx_list)[i], "box")
# }
# 
# l1_disc <- lapply(box_plot_all, function(x)x[[1]])
# multiplot(plotlist = l1_disc, cols = 3)
# l1_rep <- lapply(box_plot_all, function(x)x[[2]])
# multiplot(plotlist = l1_rep, cols = 3)
# l1_all <- c(list(l1_disc), list(l1_rep))
# 
# multiplot(plotlist = l1_all, cols = 2)
