##gene level combined p-value from rank statistics
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)

ISKS_MGRB_2020_rnd3 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv",
                                  sep = "\t", header = T, stringsAsFactors = F)

ppi_res_fil_final <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/ppi_res_fil_final_SKATbin_comb_str_biog_Aug31_cyto.tsv", 
                                sep = "\t", header = T, stringsAsFactors = F)

##function for degree detection
library(igraph)
library(ggnet)
library(intergraph)
library(network)
library(org.Hs.eg.db)

can_net <- read.delim("~/VDLab_scripts/BioGrid/biogrid_db_all_subnet.sif", header = T, sep = " ", stringsAsFactor = F)
can_net_graph <- igraph::graph.data.frame(can_net, directed = F)
can_net_graph1 <- igraph::simplify(can_net_graph, remove.loops=T, remove.multiple = T)

##changed to cytoscape stringdb
  strindb_graph1 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/strindb_cyto_graph.rds")
strindb_biog_graph1 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/strindb_biog_graph_cyto.rds")

get_deg_rand <- function(gene_list, gene_name, graph){
  
  can_net1 <- igraph::as_data_frame(graph, what = "edges")
  prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% gene_list
                                 & as.character(can_net1$to) %in% gene_list, ])
  
  if(dim(prot_np_all[prot_np_all$from %in% gene_name | prot_np_all$to %in% gene_name,])[1] != 0 ){
    uniongraph <- igraph::graph.data.frame(prot_np_all, directed = F)
    return(igraph::degree(uniongraph, gene_name, mode="all"))
  }
  else{
    #    return(NULL)
    return(0)
  }
}

##Alternate Approach; scramble PvalSKAT and PvalPPI and generate null
##this version doesn't require permuting the networks; results look fine
#for skatbin*PPI_fisher
set.seed(473839)
##fix degree zero
ppi_res_fil_final[ppi_res_fil_final$degree == 0,]$PPI_p_val_wt <- sample(seq(0.6, 0.9, 0.0001), size = 572, replace = F)
ppi_res_fil_final[ppi_res_fil_final$degree_str == 0,]$PPI_p_val_wt_str <- sample(seq(0.6, 0.9, 0.0001), size = 672, replace = F)
ppi_res_fil_final[ppi_res_fil_final$degree_str_biog == 0,]$PPI_p_val_wt_str_biog <- sample(seq(0.6, 0.9, 0.0001), size = 419, replace = F)
##score for permutation test to generate individual gene level p-value
ppi_res_fil_final$new_comp_score_rand <- -log10(ppi_res_fil_final$pval_SKATbin/ppi_res_fil_final$wt_diff)*(-log10(ppi_res_fil_final$PPI_p_val_wt))
ppi_res_fil_final$new_comp_score_str_rand <- -log10(ppi_res_fil_final$pval_SKATbin/ppi_res_fil_final$wt_diff)*(-log10(ppi_res_fil_final$PPI_p_val_wt_str))
ppi_res_fil_final$new_comp_score_str_bio_rand <- -log10(ppi_res_fil_final$pval_SKATbin/ppi_res_fil_final$wt_diff)*(-log10(ppi_res_fil_final$PPI_p_val_wt_str_biog))
#ppi_res_fil_final$new_comp_score <- -log10(ppi_res_fil_final$pval_SKATbin)*(-log10(ppi_res_fil_final$PPI_p_val_wt_str_biog))

#scramble function
ddd <- seq(1:dim(ppi_res_fil_final)[1])
seed = 3453466
scrambled_skat <- list()
for (i in 1:25){
  set.seed(seed + i)
  scrambled_skat[[i]] <- sample(x = ddd, size = dim(ppi_res_fil_final)[1], replace = FALSE)
}

seed = 4534669
scrambled_ppi <- list()
for (i in 1:25){
  set.seed(seed + i)
  scrambled_ppi[[i]] <- sample(x = ddd, size = dim(ppi_res_fil_final)[1], replace = FALSE)
}

seed = 5347669
scrambled_wt <- list()
for (i in 1:25){
  set.seed(seed + i)
  scrambled_wt[[i]] <- sample(x = ddd, size = dim(ppi_res_fil_final)[1], replace = FALSE)
}

##null distribution for composite score
null_dist_biog <- list()
null_dist_str <- list()
null_dist_str_biog <- list()
for(i in 1:25){
  ##skatbin*PPI_fisher
  #  null_dist[[i]] <- -log10(ppi_res_fil_final$pval_SKATbin[scrambled_skat[[i]]]/(ppi_res_fil_final$wt_diff[scrambled_skat[[i]]]))*(-log10(ppi_res_fil_final$PPI_p_val_wt_str_biog[scrambled_ppi[[i]]]))
  null_dist_biog[[i]] <- -log10(ppi_res_fil_final$pval_SKATbin[scrambled_skat[[i]]]/(ppi_res_fil_final$wt_diff[scrambled_wt[[i]]]))*(-log10(ppi_res_fil_final$PPI_p_val_wt[scrambled_ppi[[i]]]))
  null_dist_str[[i]] <- -log10(ppi_res_fil_final$pval_SKATbin[scrambled_skat[[i]]]/(ppi_res_fil_final$wt_diff[scrambled_wt[[i]]]))*(-log10(ppi_res_fil_final$PPI_p_val_wt_str[scrambled_ppi[[i]]]))
  null_dist_str_biog[[i]] <- -log10(ppi_res_fil_final$pval_SKATbin[scrambled_skat[[i]]]/(ppi_res_fil_final$wt_diff[scrambled_wt[[i]]]))*(-log10(ppi_res_fil_final$PPI_p_val_wt_str_biog[scrambled_ppi[[i]]]))
}
##take mean and unlist null_dist for final null dist and calculation of sd
null_dist_biog <- lapply(null_dist_biog, function(x)mean(x))
null_dist_str <- lapply(null_dist_str, function(x)mean(x))
null_dist_str_biog <- lapply(null_dist_str_biog, function(x)mean(x))
hist(unlist(null_dist_biog))
hist(unlist(null_dist_str))
hist(unlist(null_dist_str_biog))

emp_pval_comp_scores <- function(gene_sym, null_dist, network, plot_hist = NULL){
  
  gene_add_rand_rank_PPI_scores <- unlist(null_dist)
  gene_add_rand_rank_PPI_scores[which(!is.finite(gene_add_rand_rank_PPI_scores))] <- 0
  
  
  #plot(hist(unlist(gene_add_rand_rank_PPI_scores)))
  if(plot_hist == 1){
    # par(mfrow=c(1,3))
    # hist(unlist(gene_add_degree))
    #hist(unlist(gene_add_rand_rank_score))
    # hist(unlist(gene_add_rand_rank_score))
    hist(gene_add_rand_rank_PPI_scores)
  }
  #top_500 SKAT
  if(network %in% "biog"){
    orig_deg <- get_deg_rand(ppi_res_fil_final$gene, gene_sym, can_net_graph1)
    orig_rank_PPI <- ppi_res_fil_final[ppi_res_fil_final$gene %in% gene_sym,]$new_comp_score_rand
  }
  else if(network %in% "str"){
    orig_deg <- get_deg_rand(ppi_res_fil_final$gene, gene_sym, strindb_graph1)
    orig_rank_PPI <- ppi_res_fil_final[ppi_res_fil_final$gene %in% gene_sym,]$new_comp_score_str_rand
  }
  else if(network %in% "str_bio") {
    orig_deg <- get_deg_rand(ppi_res_fil_final$gene, gene_sym, strindb_biog_graph1)
    orig_rank_PPI <- ppi_res_fil_final[ppi_res_fil_final$gene %in% gene_sym,]$new_comp_score_str_bio_rand
  }
  # orig_rank_PPI <- gender_PC123_cont_df_top[gender_PC123_cont_df_top$gene %in% gene_sym,]$rank_score * log(orig_deg + 2)
  a <- orig_rank_PPI
  s <- sd(gene_add_rand_rank_PPI_scores, na.rm = T)
  #n <- length(gene_rand[[1]])
  n <- length(unlist(null_dist))
  #n <- 100
  #xbar <- mean(gene_add_rank_PPI_scores)
  xbar <- mean(gene_add_rand_rank_PPI_scores)
  t <- (xbar-a)/(s/sqrt(n))
  ##one sided
  # t_out <- t.test(gene_add_rank_PPI_scores,mu=a,alternative="less")
  t_out <- t.test(gene_add_rand_rank_PPI_scores,mu=a,alternative="less")
  
  test_stat <- cbind.data.frame("orig_deg" = orig_deg, "orig_rank_PPI" = orig_rank_PPI,
                                "mean_rand_rank_PPI" = xbar,
                                "t_stat" = t, 
                                "pval_t-test" = t_out$p.value)
  
  return(test_stat)
}

###Function for parallel processing
library(doParallel)
library(doMC)
registerDoMC(30)

##parallel computation of gene based empirical p-values
genes <- unique(as.character(ppi_res_fil_final$gene))

res20_biog <- list()
res20_str <- list()
res20_str_biog <- list()
system.time(res20_biog <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
{emp_pval_comp_scores(genes[i], null_dist_biog, "biog", 0)})
system.time(res20_str <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
{emp_pval_comp_scores(genes[i], null_dist_str, "str", 0)})
system.time(res20_str_biog <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
{emp_pval_comp_scores(genes[i], null_dist_str_biog,"str_bio", 0)})
##process output
proc_out <- function(pval_list, network){
  res20_df_no_wt_diff <- do.call("rbind.data.frame", pval_list)
  res20_df_no_wt_diff$gene <- genes
  res20_df_no_wt_diff$t_stat <- -(res20_df_no_wt_diff$t_stat)
  colnames(res20_df_no_wt_diff)[4] <- paste0(colnames(res20_df_no_wt_diff)[4], "_", network)
  colnames(res20_df_no_wt_diff)[5] <- paste0(colnames(res20_df_no_wt_diff)[5], "_", network)
  return(res20_df_no_wt_diff)
}

res20_df_biog <- proc_out(res20_biog, "biog")
res20_df_str <- proc_out(res20_str, "str")
res20_df_str_biog <- proc_out(res20_str_biog, "str_bio")
res20_df_comb <- cbind.data.frame(res20_df_biog, res20_df_str, res20_df_str_biog)

g1 <- c("TP53", "NF1", "POT1", "TINF2", "TERF1",
        "DLG1", "GRIN2A", "GRIA1", "GRM4", "CAM2KB", "CAM2KD",
        "CEP63", "CEP72", "HAUS4", "HAUS5",
        "ZBTB16", "ATG7", "ASB10", "ASB8", "FBXO7", "TRIM4",
        "PCM1",
        "CENPC", "CENPF", "BUB1", "BUB1B", "AURKB", "RANGAP1", "RACGAP1", "MAD2L2",
        "TCP1", "CCT2", "CCT6A",
        "PRPF4", "DHX9", "DHX16", "SF3A1", "CPSF3", "PPIL3",
        "L3MBTL3")

t3 <- res20_df_biog[res20_df_biog$gene %in% g1,]
t3[t3$t_stat < 0,]$gene
t4 <- res20_df_str[res20_df_str$gene %in% g1,]
t4[t4$t_stat < 0,]$gene
t5 <- res20_df_str_biog[res20_df_str_biog$gene %in% g1,]
t5[t5$t_stat < 0,]$gene

ppi_res_fil_final[,50:55] <- res20_df_comb[match(ppi_res_fil_final$gene, res20_df_comb$gene),c(4:5,10:11,16:17)]

##Combine pvalues using Fisher method
library(metap)
#ppi_res_fil_final <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/scrambled_null_ppi_fisher_res_fil_final.tsv", 
#           sep = "\t", header = T, stringsAsFactors = F)
#inp_pval <- cbind.data.frame("pval_SKATbin" = ppi_res_fil_final$pval_SKATbin ,"PPI_p_val_wt_str_biog" = ppi_res_fil_final$PPI_p_val_wt_str_biog)
#ppi_res_fil_final$PPI_p_val_wt_str_biog <- ifelse(ppi_res_fil_final$degree_str_biog == 0, 1, ppi_res_fil_final$PPI_p_val_wt_str_biog)

ppi_res_fil_final$comb_pval_biog <- apply(ppi_res_fil_final[,c(14,22)], 1, function(x)sumlog(x)$p)
ppi_res_fil_final$comb_pval_str <- apply(ppi_res_fil_final[,c(14,28)], 1, function(x)sumlog(x)$p)
ppi_res_fil_final$comb_pval_str_biog <- apply(ppi_res_fil_final[,c(14,34)], 1, function(x)sumlog(x)$p)


##rank_score1
ppi_res_fil_final$rank_score1_biog <- apply(ppi_res_fil_final[,c(14,22)], 1, function(x)sum(-log10(x)))
ppi_res_fil_final$rank_score1_str <- apply(ppi_res_fil_final[,c(14,28)], 1, function(x)sum(-log10(x)))
ppi_res_fil_final$rank_score1_str_biog <- apply(ppi_res_fil_final[,c(14,34)], 1, function(x)sum(-log10(x)))

##rank_score2

ppi_res_fil_final$rank_score2_biog <- apply(ppi_res_fil_final[,c(6,14,22)], 1, function(x)(-log10(x[2]/x[1]) + (-1)*log10(x[3])))
ppi_res_fil_final$rank_score2_str <- apply(ppi_res_fil_final[,c(6,14,28)], 1, function(x)(-log10(x[2]/x[1]) + (-1)*log10(x[3])))
ppi_res_fil_final$rank_score2_str_biog <- apply(ppi_res_fil_final[,c(6,14,34)], 1, function(x)(-log10(x[2]/x[1]) + (-1)*log10(x[3])))

write.table(ppi_res_fil_final[,c(1:11,14,17:20,22,26:28,32:34,38:61)], "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/scrambled_cyto_null_ppi_fisher_res_fil_final_Aug31_09102020.tsv",
            sep = "\t", row.names = F, quote = F)



