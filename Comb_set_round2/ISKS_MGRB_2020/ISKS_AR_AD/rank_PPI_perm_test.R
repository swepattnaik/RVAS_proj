##gene level combined p-value from rank statistics
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )


ISKS_MGRB_2020_rnd3 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/comb_ISKS_MGRB_2020_C345_gnom_0002_all_genes_AD_ver4.tsv",
           sep = "\t", header = T, stringsAsFactors = F)

ppi_res_fil_final <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/ppi_res_fil_final_SKATbin.tsv", sep = "\t", header = T)
# hist(ppi_res_fil_final$rank_score)
# hist(ppi_res_fil_final$rank_PPI)
# hist(ppi_res_fil_final$degree)

#ppi_res_fil_final$new_rank_PPI <- -log10(ppi_res_fil_final$pval_SKATbin)*ppi_res_fil_final$degree
#ppi_res_fil_final$new_rank_PPI <- -log10(ppi_res_fil_final$pval_SKATbin)*log2(ppi_res_fil_final$degree + 2)

ppi_res_fil_final <- ppi_res_fil_final[ppi_res_fil_final$ISKS > 1,]

#ppi_res_fil_final <- ppi_res_fil_final[ppi_res_fil_final$ISKS > 1 & ppi_res_fil_final$degree > 0,]

##from composite_score_bootstrap.R script

##function for degree detection
library(igraph)
library(ggnet)
library(intergraph)
library(network)
library(org.Hs.eg.db)

can_net <- read.delim("~/VDLab_scripts/BioGrid/biogrid_db_all_subnet.sif", header = T, sep = " ", stringsAsFactor = F)


###Empirical p-value based on composite scores(rank_PPI)
##Degree randomisation must include all the genes while the rank_score sampling needs to be performed
##only from the list of genes with positive wt_diff

###focus on sampling conditions below to derive a robust random null distribution
##
##sample from list of genes with positive wt_diff (for rank_score sampling)
#gene_tot_rand <- as.character(gender_PC123_cont_df$gene)
gene_tot_rand <- as.character(ppi_res_fil_final$gene)
##degree randomisation (use entire gene list)
gene_tot_rand_deg <- as.character(ISKS_MGRB_2020_rnd3$gene)
##remove the topSKAT genes from the randomised list to showcase over-representation in the topset
#gene_tot_rand <- as.character(gender_PC123_cont_df$gene)[-c(1:1000)]
##build random samples for bootstrapping w/replacement
##perhaps use all the genes not in the topSKAT for estimating rank_PPI based z-score
##rank_PPI needs to be calculated for each random network
##top 1000 genes should provide a random network for the top 500/1000 SKAT genes
##for random network sample from genes that are not in the topSKAT or use the entire set of genes (think)
#1. for random_score randomisation: use the genes that have a positive rank_score 
gene_rand <- list()
for(i in 1:100){
  #  gene_rand_ind <- sample(length(gene_tot_rand), size = 85, replace = T) ##size based on clique filtering
  gene_rand_ind <- sample(length(gene_tot_rand), size = 100, replace = F)
  gene_rand[[i]] <- gene_tot_rand[gene_rand_ind]
}
#2. for degree randomization
deg_tol_rand <- list()
for(i in 1:100){
  #  gene_rand_ind <- sample(length(gene_tot_rand), size = 85, replace = T) ##size based on clique filtering
  gene_rand_deg_ind <- sample(length(gene_tot_rand_deg), size = 1209, replace = F)
  deg_tol_rand[[i]] <- gene_tot_rand_deg[gene_rand_deg_ind]
}

##generate degree from random network
can_net_graph <- igraph::graph.data.frame(can_net, directed = F)
can_net_graph1 <- igraph::simplify(can_net_graph, remove.loops=T, remove.multiple = T)
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

###Function for parallel processing
library(doParallel)
library(doMC)
registerDoMC(30)

##precomputing permuted degree with log transform (compute once and save)
##speeds up the computation enormously

perm_deg <- function(gene_sym){
  gene_add <- lapply(deg_tol_rand, function(x)c(gene_sym, x)) ##precompute
  ##degree randomization
  gene_add_degree <- lapply(gene_add, function(x)log(get_deg_rand(x, gene_sym, can_net_graph1) + 2)) ##precompute
  return(gene_add_degree)
}

genes <- unique(as.character(ppi_res_fil_final$gene))
#g1 <- genes[1:3]
system.time(precomp_perm_net <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
{perm_deg(genes[i])})

precomp_perm_net1 <- lapply(precomp_perm_net, function(x)unlist(x))
names(precomp_perm_net1) <- genes
#saveRDS(precomp_perm_net1, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/precomp_perm_net_no_wtdiff.rds",
#        compress = T)
saveRDS(precomp_perm_net1, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/precomp_perm_net_no_wtdiff_ver2.rds",
        compress = T)
##
#load precomputed degree from permuted networks
#precomp_perm_net1 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/precomp_perm_net_no_wtdiff.rds")
precomp_perm_net1 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/precomp_perm_net_no_wtdiff_ver2.rds")


##Alternate Approach; scramble PvalSKAT and PvalPPI and generate null
##this version doesn't require permuting the networks; results look fine
#for skatbin*PPI_fisher

ppi_res_fil_final[ppi_res_fil_final$degree == 0,]$PPI_p_val_wt <- sample(seq(0.35, 0.7, 0.0001), size = 565, replace = F)
ppi_res_fil_final$new_comp_score <- -log10(ppi_res_fil_final$pval_SKATbin)*(-log10(ppi_res_fil_final$PPI_p_val_wt))

#ppi_res_fil_final$new_comp_score <- -log10(ppi_res_fil_final$pval_SKATbin)+ abs(-log10(ppi_res_fil_final$PPI_p_val_wt))
#ppi_res_fil_final$degree <- ifelse(ppi_res_fil_final$degree == 0, 2, ppi_res_fil_final$PPI_p_val_wt)

##skatbin*degree log transform
#ppi_res_fil_final$new_comp_score <- ifelse(ppi_res_fil_final$degree == 0, -log10(ppi_res_fil_final$pval_SKATbin),
#                                           -log10(ppi_res_fil_final$pval_SKATbin)*(log2(ppi_res_fil_final$degree + 2)))

##skatbin*degree no transform
#ppi_res_fil_final$degree <- ppi_res_fil_final$degree + 1
#ppi_res_fil_final$new_comp_score <- -log10(ppi_res_fil_final$pval_SKATbin)*(ppi_res_fil_final$degree)

#scramble function
ddd <- seq(1:dim(ppi_res_fil_final)[1])
seed = 3453466
scrambled_skat <- list()
for (i in 1:25){
  set.seed(seed + i)
  #scrambled_skat[[i]] <- sample(x = ddd, size = 100, replace = FALSE)
  scrambled_skat[[i]] <- sample(x = ddd, size = dim(ppi_res_fil_final)[1], replace = FALSE)
}
  
seed = 4534669
scrambled_ppi <- list()
for (i in 1:25){
  set.seed(seed + i)
  #scrambled_ppi[[i]] <- sample(x = ddd, size = 100, replace = FALSE)
  scrambled_ppi[[i]] <- sample(x = ddd, size = dim(ppi_res_fil_final)[1], replace = FALSE)
}

##null distribution for composite score
null_dist <- list()
for(i in 1:25){
  ##skatbin*PPI_fisher
  null_dist[[i]] <- -log10(ppi_res_fil_final$pval_SKATbin[scrambled_skat[[i]]])*(-log10(ppi_res_fil_final$PPI_p_val_wt[scrambled_ppi[[i]]]))

 
  ##skatbin*degree log transform
  # null_dist[[i]] <- -log10(ppi_res_fil_final$pval_SKATbin[scrambled_skat[[i]]])*
 #   (log2(ppi_res_fil_final$degree[scrambled_ppi[[i]]] + 2))

  ##skatbin*degree no transform  
#  null_dist[[i]] <- -log10(ppi_res_fil_final$pval_SKATbin[scrambled_skat[[i]]])*(ppi_res_fil_final$degree[scrambled_ppi[[i]]])
  
}
##take mean and unlist null_dist for final null dist and calculation of sd
null_dist <- lapply(null_dist, function(x)mean(x))
hist(unlist(null_dist))
##New version
##add options for rank scores: with/without wt_diff
emp_pval_comp_scores <- function(gene_sym, opt = NULL, plot_hist = NULL){
  ##add gene of interest to compute degree in randomly sampled genes
#  gene_add <- lapply(deg_tol_rand, function(x)c(gene_sym, x)) ##precompute
  ##degree randomization
 # gene_add_degree <- lapply(gene_add, function(x)log(get_deg_rand(x, gene_sym) + 2)) ##precompute
  
  gene_add_degree <- precomp_perm_net1[gene_sym] 
#  gene_add_degree <- lapply(gene_add, function(x)get_deg_rand(x, gene_sym))
  ##sampling distribution from rank score (derived from rankscores with positive wt_diff)
  #gene_add_rand_rank_score <- lapply(gene_rand, function(x)mean(gender_PC123_cont_df[gender_PC123_cont_df$gene %in% x,]$rank_score))
  ##opt == 1 corresponds to score with wt_diff
  if(opt == 1){
  gene_add_rand_rank_score <- lapply(gene_rand, function(x)mean(ppi_res_fil_final[ppi_res_fil_final$gene %in% x,]$rank_score))
  gene_add_rand_rank_PPI_scores <- unlist(gene_add_degree)*unlist(gene_add_rand_rank_score)
  gene_add_rand_rank_PPI_scores[which(!is.finite(gene_add_rand_rank_PPI_scores))] <- 0
  orig_rank_PPI <- ppi_res_fil_final[ppi_res_fil_final$gene %in% gene_sym,]$rank_PPI
  
  }
  ##opt == 2 corresponds to score without wt_diff
  else if (opt == 2){
    gene_add_rand_rank_score <- lapply(gene_rand, function(x)mean(-log10(ppi_res_fil_final[ppi_res_fil_final$gene %in% x,]$pval_SKATbin)))
    gene_add_rand_rank_PPI_scores <- unlist(gene_add_degree)*unlist(gene_add_rand_rank_score)
    gene_add_rand_rank_PPI_scores[which(!is.finite(gene_add_rand_rank_PPI_scores))] <- 0
    orig_rank_PPI <- ppi_res_fil_final[ppi_res_fil_final$gene %in% gene_sym,]$new_rank_PPI
    }
    #gene_add_rand_rank_score <- lapply(gene_add, function(x)mean(gender_PC123_cont_df[gender_PC123_cont_df$gene %in% x,]$rank_score))
  else if (opt == 3){
    gene_add_rand_rank_PPI_scores <- unlist(null_dist)
    gene_add_rand_rank_PPI_scores[which(!is.finite(gene_add_rand_rank_PPI_scores))] <- 0
    orig_rank_PPI <- ppi_res_fil_final[ppi_res_fil_final$gene %in% gene_sym,]$new_comp_score
  }
 
  #plot(hist(unlist(gene_add_rand_rank_PPI_scores)))
  if(plot_hist == 1){
 # par(mfrow=c(1,3))
 # hist(unlist(gene_add_degree))
  #hist(unlist(gene_add_rand_rank_score))
 # hist(unlist(gene_add_rand_rank_score))
  hist(gene_add_rand_rank_PPI_scores)
  }
  #top_500 SKAT
  orig_deg <- get_deg_rand(ppi_res_fil_final$gene, gene_sym, can_net_graph1)
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
  if (opt == 1 | opt == 2){
  test_stat <- cbind.data.frame("orig_deg" = orig_deg, "orig_rank_PPI" = orig_rank_PPI, 
                                "mean_rand_rank_score" = mean(unlist(gene_add_rand_rank_score)),
                                "mean_rand_rank_PPI" = xbar,
                                "t_stat" = abs(t), 
                                "pval_t-test" = t_out$p.value)
  }
  else{
    test_stat <- cbind.data.frame("orig_deg" = orig_deg, "orig_rank_PPI" = orig_rank_PPI,
                                  "mean_rand_rank_PPI" = xbar,
                                  "t_stat" = t, 
                                  "pval_t-test" = t_out$p.value)
  }
  return(test_stat)
}

##parallel computation of gene based empirical p-values
genes <- unique(as.character(ppi_res_fil_final$gene))
# system.time(res20 <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
# {emp_pval_comp_scores(genes[i], 1, 0)})
# 
# res20_df_wt_diff <- do.call("rbind.data.frame", res20)
# res20_df_wt_diff$gene <- rownames(res20_df_wt_diff)
# ppi_res_fil_final[,27:29] <- res20_df_wt_diff[match(res20_df_wt_diff$gene, ppi_res_fil_final$gene),c(4:6)]

res20 <- list()
system.time(res20 <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
{emp_pval_comp_scores(genes[i], 3, 0)})
res20_df_no_wt_diff <- do.call("rbind.data.frame", res20)
res20_df_no_wt_diff$gene <- genes
res20_df_no_wt_diff$t_stat <- -(res20_df_no_wt_diff$t_stat)
#res20_df_no_wt_diff$gene <- rownames(res20_df_no_wt_diff)
##QC
# t2 <- res20_df_no_wt_diff[res20_df_no_wt_diff$gene %in% c("TP53", "NF1", "POT1", "TINF2", "TERF1", 
#                                                     "DLG1", "GRIN2A", "GRIA1", "CAM2KB", "CAM2KD",
#                                                     "CEP63", "CEP72", "HAUS4", "HAUS5",
#                                                     "ZBTB16", "ATG7", "ASB10", "ASB8", "FBXO7", "TRIM4",
#                                                     "PCM1",
#                                                     "CENPC", "CENPF", "BUB1", "BUB1B", "AURKB", "RANGAP1", "RACGAP1", "MAD2L2",
#                                                     "TCP1", "CCT2", "CCT6A",
#                                                     "PRPF4", "DHX9", "DHX16", "SF3A1", "CPSF3", "PPIL3",
#                                                     "L3MBTL3"),]
t3 <- res20_df_no_wt_diff[res20_df_no_wt_diff$gene %in% c("TP53", "NF1", "POT1", "TINF2", "TERF1",
                                                          "DLG1", "GRIN2A", "GRIA1", "CAM2KB", "CAM2KD",
                                                          "CEP63", "CEP72", "HAUS4", "HAUS5",
                                                          "ZBTB16", "ATG7", "ASB10", "ASB8", "FBXO7", "TRIM4",
                                                          "PCM1",
                                                          "CENPC", "CENPF", "BUB1", "BUB1B", "AURKB", "RANGAP1", "RACGAP1", "MAD2L2",
                                                          "TCP1", "CCT2", "CCT6A",
                                                          "PRPF4", "DHX9", "DHX16", "SF3A1", "CPSF3", "PPIL3",
                                                          "L3MBTL3"),]
# lost_degree_no_log <- t2[t2$t_stat < 0,]$gene
lost_degree_PPI_fisher <- t3[t3$t_stat < 0,]$gene
#ppi_res_fil_final[,30:32] <- res20_df_no_wt_diff[match(res20_df_no_wt_diff$gene, ppi_res_fil_final$gene),c(4:6)]
#colnames(ppi_res_fil_final)[27:29] <- paste0("wt_", colnames(ppi_res_fil_final)[27:29])
#write.table(ppi_res_fil_final, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/permute_pval_ppi_res_fil_final.tsv",
#            sep = "\t", row.names = F, quote = F)
ppi_res_fil_final[,27:29] <- res20_df_no_wt_diff[match(res20_df_no_wt_diff$gene, ppi_res_fil_final$gene),c(4:6)]

##Combine pvalues using Fisher method
library(metap)
#ppi_res_fil_final <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/scrambled_null_ppi_fisher_res_fil_final.tsv", 
#           sep = "\t", header = T, stringsAsFactors = F)
inp_pval <- cbind.data.frame("pval_SKATbin" = ppi_res_fil_final$pval_SKATbin ,"PPI_p_val_wt" = ppi_res_fil_final$PPI_p_val_wt)
ppi_res_fil_final$PPI_p_val_wt <- ifelse(ppi_res_fil_final$degree == 0, 1, ppi_res_fil_final$PPI_p_val_wt)
comb_pval <- apply(inp_pval, 1, function(x)sumlog(x)$p)
ppi_res_fil_final$comb_fisher <- comb_pval

#write.table(ppi_res_fil_final[,c(1:12,14,17:20,22,26:28)], "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/scrambled_null_ppi_res_fil_final.tsv",
#            sep = "\t", row.names = F, quote = F)
write.table(ppi_res_fil_final[,c(1:12,14,17:20,22,26:28,30)], "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/scrambled_null_ppi_fisher_res_fil_final.tsv",
            sep = "\t", row.names = F, quote = F)
