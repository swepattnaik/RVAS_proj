##use pnorm
##Induce graph including random networks and enriched genes
##compute centrality/degree
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Exome_skat_wsing_load123_noCH_C5eqC4_nonmds_gt_isksrisc_sept05_rect_ASP_cmaf_new_score.rds")

df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})
Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Exome_pc123_srt_SKAT_case_enr_nCH_iskrisc_05Sept_rect_ASP_cmaf_new_score.rds")
##Add p-value
Exome_pc123_srt_SKAT_case_enr_nCH_pval <- Map(cbind.data.frame, Exome_pc123_srt_SKAT_case_enr_nCH, Exome_pc123_srt_SKAT)
##extract genes to sample from
#gene_tot <- names(V(can_net_graph1))
##enriched genes
Exome_pc123_srt_SKAT_case_enr_nCH_pval01 <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_pval, function(x)
  x[x$pval_SKATbin < 0.1 & x$wt_diff > 0,])
Exome_pc123_srt_SKAT_case_enr_nCH_pvalgt01 <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_pval, function(x)
  x[x$pval_SKATbin > 0.1 & x$wt_diff > 0,])
Exome_pc123_srt_SKAT_case_neg_enr_nCH_pval <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_pval, function(x)
  x[x$wt_diff < 0,])
#Exome_pc123_srt_SKAT_case_pos_enr_nCH_pval <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_pval, function(x)
#  x[x$wt_diff > 0,])
##for sane log2 wt_diff values
Exome_pc123_srt_SKAT_case_pos_enr_nCH_pval <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_pval, function(x)
  x[x$wt_diff > 1,])
##function for gene specific empirical p-value
# Exome_pc123_srt_SKAT_case_enr_nCH_ranked <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_pval01, 
#                                                    function(x) cbind.data.frame(x, "rank_score" = -log10(x[,14]) * x[,6])) ##added intra cohort MAF filtered variant positions
# 
# Exome_pc123_srt_SKAT_case_enr_nCH_ranked <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_ranked, function(x)
#   x[order(x[,17], decreasing = T),])
# gender_PC123_cont_df <- Exome_pc123_srt_SKAT_case_enr_nCH_ranked[[4]]
# 
# Exome_pc123_srt_SKAT_case_enr_nCH_ranked <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_pval01, 
#                                                    function(x) cbind.data.frame(x, "rank_score" = -log10(x[,14]) * x[,6])) ##added intra cohort MAF filtered variant positions

##don't log transform the degree
##positive wt_diff
Exome_pc123_srt_SKAT_case_pos_enr_nCH_ranked <- lapply(Exome_pc123_srt_SKAT_case_pos_enr_nCH_pval, 
                                                   function(x) cbind.data.frame(x, "rank_score" = -log10(x[,14]) * x[,6]))
Exome_pc123_srt_SKAT_case_pos_enr_nCH_ranked <- lapply(Exome_pc123_srt_SKAT_case_pos_enr_nCH_ranked, function(x)
  x[order(x[,17], decreasing = T),])

##all rank_score to entire list
Exome_pc123_srt_SKAT_case_all_enr_nCH_ranked <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_pval, 
                                                       function(x) cbind.data.frame(x, "rank_score" = -log10(x[,14]) * x[,6]))
Exome_pc123_srt_SKAT_case_all_enr_nCH_ranked <- lapply(Exome_pc123_srt_SKAT_case_all_enr_nCH_ranked, function(x)
  x[order(x[,17], decreasing = T),])

##Add rank-score to all genes
#Exome_pc123_srt_SKAT_case_enr_nCH_ranked_all <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_pval, 
#                                                       function(x) cbind.data.frame(x, "rank_score" = -log10(x[,14]) * log(x[,6])))
                                                       
gender_PC123_cont_df <- Exome_pc123_srt_SKAT_case_pos_enr_nCH_ranked[[4]]
gender_PC123_cont_df_all <- Exome_pc123_srt_SKAT_case_all_enr_nCH_ranked[[4]]

##function for degree detection
library(igraph)
library(ggnet)
library(intergraph)
library(network)
library(org.Hs.eg.db)

can_net <- read.delim("~/VDLab_scripts/BioGrid/biogrid_db_all_subnet.sif", header = T, sep = " ", stringsAsFactor = F)

get_topSKAT_degree <- function(enr_df){
  
  skat_genes_top100 <- as.character(enr_df[,3])
  
  can_net <- read.delim("~/VDLab_scripts/BioGrid/biogrid_db_all_subnet.sif", header = T, sep = " ", stringsAsFactor = F)
  #can_net <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Biogrid_latest/Biog_net_hs.sif", header = T, sep = "\t", stringsAsFactor = F)
  can_net_graph <- igraph::graph.data.frame(can_net, directed = F)
  can_net_graph1 <- igraph::simplify(can_net_graph, remove.loops=T, remove.multiple = T)
  can_net1 <- igraph::as_data_frame(can_net_graph1, what = "edges")
  prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% skat_genes_top100 
                                 & as.character(can_net1$to) %in% skat_genes_top100, ])
  uniongraph <- igraph::graph.data.frame(prot_np_all, directed = F)
  
  
  
  ret_df <- as.data.frame(igraph::degree(uniongraph))
  ret_df$gene <- rownames(ret_df)
  colnames(ret_df)[1] <- c("degree")
  #Biogrid_df <- as.data.frame(igraph::degree(can_net_graph1)) 
  #Biogrid_df$gene <- rownames(Biogrid_df)
  #colnames(Biogrid_df)[1] <- c("degree")
  #ret_df$bioG_deg <- Biogrid_df[match(ret_df$gene, Biogrid_df$gene), 1]
  #ret_df$PPI_ratio <- (ret_df$degree)/(ret_df$bioG_deg)
  
  l1 <- lapply(ret_df[,2], function(x)as_ids(adjacent_vertices(uniongraph, v = x)[[1]]))
  l2 <- lapply(l1, function(x)paste(x, sep=",", collapse=","))
  ret_df$interactors <- unlist(l2)
  
  return(ret_df)
}


#for gender_PC123_cont_df
df_degree <- get_topSKAT_degree(gender_PC123_cont_df)
gender_PC123_cont_df$degree <- df_degree[match(gender_PC123_cont_df$gene, df_degree$gene), 1]
gender_PC123_cont_df$int <- df_degree[match(gender_PC123_cont_df$gene, df_degree$gene), 3]
gender_PC123_cont_df$degree <- ifelse(is.na(gender_PC123_cont_df$degree), 0, gender_PC123_cont_df$degree)
gender_PC123_cont_df$rank_PPI <- gender_PC123_cont_df$rank_score * log(gender_PC123_cont_df$degree + 2)
#gender_PC123_cont_df$PPI_ratio <- df_degree[match(gender_PC123_cont_df$gene, df_degree$gene), 4]
#gender_PC123_cont_df$rank_PPI <- gender_PC123_cont_df$rank_score * (gender_PC123_cont_df$PPI_ratio))
gender_PC123_cont_df <- gender_PC123_cont_df[order(gender_PC123_cont_df$rank_PPI, decreasing = T),]

#for gender_PC123_cont_df_all (not needed)
# df_degree <- get_topSKAT_degree(gender_PC123_cont_df_all)
# gender_PC123_cont_df_all$degree <- df_degree[match(gender_PC123_cont_df_all$gene, df_degree$gene), 1]
# gender_PC123_cont_df_all$int <- df_degree[match(gender_PC123_cont_df_all$gene, df_degree$gene), 3]
# gender_PC123_cont_df_all$degree <- ifelse(is.na(gender_PC123_cont_df_all$degree), 0, gender_PC123_cont_df_all$degree)

#############
##remove from gene_tot for bootstrapping: substitute for Fisher's test of PPI enrichment
##This is more robust
#gene_tot <- c(as.character(Exome_pc123_srt_SKAT_case_neg_enr_nCH_pval[[4]]$gene), 
#              as.character(Exome_pc123_srt_SKAT_case_enr_nCH_pvalgt01[[4]]$gene))
#gene_tot_rand <- gene_tot[!(gene_tot %in% as.character(Exome_pc123_srt_SKAT_case_enr_nCH_pval01[[4]]$gene))]
gene_tot_rand <- c(as.character(Exome_pc123_srt_SKAT_case_neg_enr_nCH_pval[[4]]$gene), 
                   as.character(Exome_pc123_srt_SKAT_case_enr_nCH_pvalgt01[[4]]$gene))
##build random samples for bootstrapping w/replacement
#rand_samp <- sapply(gene_tot, function(x) {sample(1:1000, size = 905, prob = 1/905)})
gene_rand <- list()
for(i in 1:50){
  #  gene_rand_ind <- sample(length(gene_tot_rand), size = 85, replace = T) ##size based on clique filtering
  gene_rand_ind <- sample(length(gene_tot_rand), size = 500, replace = F)
  gene_rand[[i]] <- gene_tot[gene_rand_ind]
}

##function for report degree of specific vertex of a subgraph
can_net <- read.delim("~/VDLab_scripts/BioGrid/biogrid_db_all_subnet.sif", header = T, sep = " ", stringsAsFactor = F)

##generate degree from random network
get_deg_rand <- function(gene_list, gene_name){
  can_net_graph <- igraph::graph.data.frame(can_net, directed = F)
  can_net_graph1 <- igraph::simplify(can_net_graph, remove.loops=T, remove.multiple = T)
  can_net1 <- igraph::as_data_frame(can_net_graph1, what = "edges")
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

##degree based enrichment of gene (does not sample from rank_score distribution)
emp_pval_gene <- function(gene_sym){

  gene_add <- lapply(gene_rand, function(x)c(gene_sym, x))
  gene_add_degree <- lapply(gene_add, function(x)get_deg_rand(x, gene_sym))
  plot(hist(unlist(gene_add_degree)))
  a <- gender_PC123_cont_df[gender_PC123_cont_df$gene %in% gene_sym,]$degree
  s <- sd(unlist(gene_add_degree), na.rm = T)
  n <- length(gene_rand[[1]])
  xbar <- mean(unlist(gene_add_degree))
  t <- (xbar-a)/(s/sqrt(n))
  ##one sided
  t_out <- t.test(unlist(gene_add_degree),mu=a,alternative="less")
  return(t_out$p.value)
}

#################Start here############
###Empirical p-value based on composite scores(rank_PPI)
gene_tot_rand <- as.character(gender_PC123_cont_df$gene)
##Degree randomisation must include all the genes while the rank_score sampling needs to be performed
##only from the list of genes with positive wt_diff

###focus on sampling conditions below to derive a robust random null distribution
##
##sample from list of genes with positive wt_diff (for rank_score sampling)
#gene_tot_rand <- as.character(gender_PC123_cont_df$gene)
gene_tot_rand <- as.character(gender_PC123_cont_df$gene)
##degree randomisation (use entire gene list)
gene_tot_rand_deg <- as.character(gender_PC123_cont_df_all$gene)
##remove the topSKAT genes from the randomised list to showcase over-representation in the topset
#gene_tot_rand <- as.character(gender_PC123_cont_df$gene)[-c(1:1000)]
##build random samples for bootstrapping w/replacement
##perhaps use all the genes not in the topSKAT for estimating rank_PPI based z-score
##rank_PPI needs to be calculated for each random network
##top 1000 genes should provide a random network for the top 500/1000 SKAT genes
##for random network sample from genes that are not in the topSKAT or use the entire set of genes (think)
#1. for random_score randomisation: use the genes that have a positive rank_score 
gene_rand <- list()
for(i in 1:50){
  #  gene_rand_ind <- sample(length(gene_tot_rand), size = 85, replace = T) ##size based on clique filtering
  gene_rand_ind <- sample(length(gene_tot_rand), size = 1000, replace = F)
  gene_rand[[i]] <- gene_tot_rand[gene_rand_ind]
}
#2. for degree randomization
deg_tol_rand <- list()
for(i in 1:50){
  #  gene_rand_ind <- sample(length(gene_tot_rand), size = 85, replace = T) ##size based on clique filtering
  gene_rand_deg_ind <- sample(length(gene_tot_rand_deg), size = 1000, replace = F)
  deg_tol_rand[[i]] <- gene_tot_rand_deg[gene_rand_deg_ind]
}

##generate degree from random network
get_deg_rand <- function(gene_list, gene_name){
  can_net_graph <- igraph::graph.data.frame(can_net, directed = F)
  can_net_graph1 <- igraph::simplify(can_net_graph, remove.loops=T, remove.multiple = T)
  can_net1 <- igraph::as_data_frame(can_net_graph1, what = "edges")
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
##composite score generation (rankPPI scores) ##not necessary

# emp_pval_comp_scores <- function(gene_sym){
#   
#  # gene_add <- lapply(gene_rand, function(x)c(gene_sym, x))
#   gene_add_degree <- lapply(gene_rand, function(x)mean(gender_PC123_cont_df[gender_PC123_cont_df$gene %in% x,]$rank_PPI))
#   plot(hist(unlist(gene_add_degree)))
#   a <- gender_PC123_cont_df[gender_PC123_cont_df$gene %in% gene_sym,]$rank_PPI
#   s <- sd(unlist(gene_add_degree), na.rm = T)
#   n <- length(gene_rand[[1]])
#   xbar <- mean(unlist(gene_add_degree))
#   t <- (xbar-a)/(s/sqrt(n))
#   ##one sided
#   t_out <- t.test(unlist(gene_add_degree),mu=a,alternative="less")
#   return(t_out$p.value)
# }

##compute rank_PPI on the fly

##New version
emp_pval_comp_scores <- function(gene_sym){
  ##add gene of interest to compute degree in randomly sampled genes
  gene_add <- lapply(deg_tol_rand, function(x)c(gene_sym, x))
  ##degree randomization
  gene_add_degree <- lapply(gene_add, function(x)get_deg_rand(x, gene_sym))
  #gene_add_degree_log <- log(unlist(gene_add_degree))
  ##sampling distribution from rank score (derived from rankscores with positive wt_diff)
  gene_add_rand_rank_score <- lapply(gene_rand, function(x)mean(gender_PC123_cont_df[gender_PC123_cont_df$gene %in% x,]$rank_score))
  #gene_add_rand_rank_score <- lapply(gene_add, function(x)mean(gender_PC123_cont_df[gender_PC123_cont_df$gene %in% x,]$rank_score))
  
  #plot(hist(unlist(gene_add_rand_rank_score)))
  #gene_add_rank_PPI <- lapply(gene_add_degree, function(x) log(x) * gender_PC123_cont_df[gender_PC123_cont_df$gene %in% gene_sym,]$rank_score)
  #gene_add_rank_PPI_scores <- unlist(gene_add_rank_PPI)
  #gene_add_rank_PPI_scores[which(!is.finite(gene_add_rank_PPI_scores))] <- 0
  #plot(hist(gene_add_rank_PPI_scores))
  gene_add_rand_rank_PPI <- Map('*',gene_add_degree,gene_add_rand_rank_score)
  gene_add_rand_rank_PPI_scores <- unlist(gene_add_rand_rank_PPI)
  gene_add_rand_rank_PPI_scores[which(!is.finite(gene_add_rand_rank_PPI_scores))] <- 0
  #plot(hist(unlist(gene_add_rand_rank_PPI_scores)))
  par(mfrow=c(1,3))
  hist(unlist(gene_add_degree))
  hist(unlist(gene_add_rand_rank_score))
  hist(gene_add_rand_rank_PPI_scores)
  #top_500 SKAT
  gender_PC123_cont_df_top <- gender_PC123_cont_df[1:1000,]
  orig_deg <- get_deg_rand(gender_PC123_cont_df_top$gene, gene_sym)
  orig_rank_PPI <- gender_PC123_cont_df_top[gender_PC123_cont_df_top$gene %in% gene_sym,]$rank_score * log(orig_deg)
  a <- orig_rank_PPI
  s <- sd(gene_add_rand_rank_PPI_scores, na.rm = T)
  n <- length(gene_rand[[1]])
  #xbar <- mean(gene_add_rank_PPI_scores)
  xbar <- mean(gene_add_rand_rank_PPI_scores)
  t <- (xbar-a)/(s/sqrt(n))
  ##one sided
  # t_out <- t.test(gene_add_rank_PPI_scores,mu=a,alternative="less")
  t_out <- t.test(gene_add_rand_rank_PPI_scores,mu=a,alternative="less")
  return(t_out$p.value)
}




##Old version
###Empirical p-value based on composite scores(rank_PPI)
gender_PC123_cont_df <- gender_PC123_cont_df[gender_PC123_cont_df$rank_score > 1,] #to avoid negative logarithm values
gene_tot_rand <- as.character(gender_PC123_cont_df$gene)
##Degree randomisation must include all the genes while the rank_score sampling needs to be performed
##only from the list of genes with positive wt_diff
gene_rand <- list()
for(i in 1:50){
  #  gene_rand_ind <- sample(length(gene_tot_rand), size = 85, replace = T) ##size based on clique filtering
  gene_rand_ind <- sample(length(gene_tot_rand), size = 1000, replace = F)
  gene_rand[[i]] <- gene_tot_rand[gene_rand_ind]
}

emp_pval_comp_scores <- function(gene_sym){
  
  gene_add <- lapply(gene_rand, function(x)c(gene_sym, x))
  ##degree randomization
  gene_add_degree <- lapply(gene_add, function(x)get_deg_rand(x, gene_sym))
  ##sampling distribution from rank score
  gene_add_rand_rank_score <- lapply(gene_add, function(x)mean(log(gender_PC123_cont_df[gender_PC123_cont_df$gene %in% x,]$rank_score)))
  #plot(hist(unlist(gene_add_rand_rank_score)))
  #gene_add_rank_PPI <- lapply(gene_add_degree, function(x) log(x) * gender_PC123_cont_df[gender_PC123_cont_df$gene %in% gene_sym,]$rank_score)
  #gene_add_rank_PPI_scores <- unlist(gene_add_rank_PPI)
  #gene_add_rank_PPI_scores[which(!is.finite(gene_add_rank_PPI_scores))] <- 0
  #plot(hist(gene_add_rank_PPI_scores))
  gene_add_rand_rank_PPI <- Map('*',gene_add_degree,gene_add_rand_rank_score)
  gene_add_rand_rank_PPI_scores <- unlist(gene_add_rand_rank_PPI)
  gene_add_rand_rank_PPI_scores[which(!is.finite(gene_add_rand_rank_PPI_scores))] <- 0
  
  #plot(hist(unlist(gene_add_rand_rank_PPI_scores)))
  par(mfrow=c(1,3))
  hist(unlist(gene_add_degree))
  hist(unlist(gene_add_rand_rank_score))
  hist(gene_add_rand_rank_PPI_scores)
  #top_500 SKAT
  gender_PC123_cont_df_top <- gender_PC123_cont_df[1:1000,]
  orig_deg <- get_deg_rand(gender_PC123_cont_df_top$gene, gene_sym)
  orig_rank_PPI <- log(gender_PC123_cont_df_top[gender_PC123_cont_df_top$gene %in% gene_sym,]$rank_score) * orig_deg
#  orig_rank_PPI <- gender_PC123_cont_df_top[gender_PC123_cont_df_top$gene %in% gene_sym,]$rank_score * log(orig_deg)
  a <- orig_rank_PPI
  s <- sd(gene_add_rand_rank_PPI_scores, na.rm = T)
  n <- length(gene_rand[[1]])
  #xbar <- mean(gene_add_rank_PPI_scores)
  xbar <- mean(gene_add_rand_rank_PPI_scores)
  t <- (xbar-a)/(s/sqrt(n))
  ##one sided
 # t_out <- t.test(gene_add_rank_PPI_scores,mu=a,alternative="less")
  t_out <- t.test(gene_add_rand_rank_PPI_scores,mu=a,alternative="less")
  return(t_out$p.value)
}


emp_pval_comp_scores("TP53")
emp_pval_comp_scores("TERF1")
emp_pval_comp_scores("POT1")
emp_pval_comp_scores("TINF2")

##use fisher's p-value to generate composite score
#use function in summary document
#top_500 SKAT
gender_PC123_cont_df_pval01 <- gender_PC123_cont_df[gender_PC123_cont_df$pval_SKATbin < 0.1,]
top_genes_enr <- as.character(gender_PC123_cont_df_pval01$gene)
para_fisher_fun_new <- function(gene_sym){
  gene_sym <- gene_sym[3]
  # gene_sym <- as.character(gene_sym[,3]) #for serial loop
  print(gene_sym)
  can_net_graph <- igraph::graph.data.frame(can_net, directed = F)
  can_net_graph1 <- igraph::simplify(can_net_graph, remove.loops=T, remove.multiple = T)
  can_net1 <- igraph::as_data_frame(can_net_graph1, what = "edges")
  prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% top_genes_enr 
                                 & as.character(can_net1$to) %in% top_genes_enr, ])
  uniongraph <- igraph::graph.data.frame(prot_np_all, directed = F)
  
  ##degree biogrid
  deg_biog <- igraph::degree(can_net_graph1)[which(names(igraph::degree(can_net_graph1)) %in% gene_sym)]
  deg_union <- igraph::degree(uniongraph)[which(names(igraph::degree(uniongraph)) %in% gene_sym)]
  #tot_biog <- igraph::vcount(can_net_graph1) ##should be edge count not vertex count
  #tot_union <- igraph::vcount(uniongraph)
  tot_biog <- igraph::ecount(can_net_graph1) 
  tot_union <- igraph::ecount(uniongraph)
  
  ##test
  if(length(deg_union) > 0 & length(deg_biog) > 0){
    inp <- c(deg_union, deg_biog, tot_union - deg_union, tot_biog - deg_biog)
    mgrb_tsg <- matrix(inp ,nrow = 2, ncol = 2)
    colnames(mgrb_tsg) <- c("deg_enr", "deg_bio")
    rownames(mgrb_tsg) <- c("Enriched", "Biog")
    #ft <- fisher.test(mgrb_tsg, alternative = "greater")
    ft <- fisher.test(mgrb_tsg)
    ft_df <- cbind.data.frame("gene" =  gene_sym, "PPI_p_val_wt" = ft$p.value,
                              "CI_lower" = ft$conf.int[1],
                              "CI_upper" = ft$conf.int[2],
                              "OR" = ft$estimate)
    
    return(ft_df)
  }else if(length(deg_biog) == 0 | length(deg_union) == 0){
    ft_df <- cbind.data.frame("gene" =  gene_sym, "PPI_p_val_wt" = 1,
                              "CI_lower" = NA,
                              "CI_upper" = NA,
                              "OR" = NA)
    return(ft_df)
  }
}



library(parallel)
cl <- makeCluster(25)
#para_pheno <- parLapply(cl, 1:nrow(IBD_dist1), function(i) pick_pinf(IBD_dist1[i,]))
clusterExport(cl, c("can_net", "top_genes_enr", "para_fisher_fun_new"))
system.time(para_fish <- parApply(cl, gender_PC123_cont_df_pval01, 1, para_fisher_fun_new))
stopCluster(cl)
fisher_res <- do.call("rbind.data.frame",para_fish)
ppi_res_fil_final <- cbind.data.frame(gender_PC123_cont_df_pval01, fisher_res)

ppi_res_fil_final$comp_score <- -log10(ppi_res_fil_final$pval_SKATbin * ppi_res_fil_final$PPI_p_val_wt)
ppi_res_fil_final <- ppi_res_fil_final[order(ppi_res_fil_final$comp_score, decreasing = T),]
#t <- (xbar-a)/(s/sqrt(n))
ppi_res_fil_final$comp_Z_score <- (mean(ppi_res_fil_final$comp_score) - ppi_res_fil_final$comp_score)/(sd(ppi_res_fil_final$comp_score)/sqrt(length(ppi_res_fil_final$comp_score)))
ppi_res_fil_final$comp_score_pval <- pnorm(ppi_res_fil_final$comp_Z_score, lower.tail = T)
##diagnostics: qqplot
library("car")
qqPlot(ppi_res_fil_final$comp_score)
hist(ppi_res_fil_final$comp_score)
table(is.na(ppi_res_fil_final$OR))
table(ppi_res_fil_final$degree == 0)

# sample from H0 separately, no assumption about equal variance
p_val_vec <- c()
for(k in 1:length(ppi_res_fil_final$comp_score)){
#  newscore <- ppi_res_fil_final$comp_score - mean(ppi_res_fil_final$comp_score) + ppi_res_fil_final$comp_score[k]
bstrap <- c()
for (i in 1:1000){
newsample <- sample(ppi_res_fil_final$comp_score, 40, replace=T)
bstrap <- c(bstrap, mean(newsample))}
#hist(bstrap)
#qqnorm(bstrap)
#qqline(bstrap)
p_val_vec <- c(p_val_vec, t.test(bstrap,alternative="less",mu=ppi_res_fil_final$comp_score[k])$p.value)
}

qqnorm(-log10(p_val_vec))

##try network propogation approaches for gene prioritization