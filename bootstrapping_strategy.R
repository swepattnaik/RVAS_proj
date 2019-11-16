##bootstrapping PPI enrichment
##pool of 1000 random gene samples (with replacement)##build this once
#loop through random genes and induce subgraph for each of the top 906 genes
#derive empirical p-value for test statistic (degree/centrality here)
##ensure that the number of genes sampled is equivalent to the list that is bootstrapped
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
##remove from gene_tot for bootstrapping
gene_tot <- c(as.character(Exome_pc123_srt_SKAT_case_neg_enr_nCH_pval[[4]]$gene), 
              as.character(Exome_pc123_srt_SKAT_case_enr_nCH_pvalgt01[[4]]$gene))
gene_tot_rand <- gene_tot[!(gene_tot %in% as.character(Exome_pc123_srt_SKAT_case_enr_nCH_pval01[[4]]$gene))]
##build random samples for bootstrapping w/replacement
#rand_samp <- sapply(gene_tot, function(x) {sample(1:1000, size = 905, prob = 1/905)})
gene_rand <- list()
for(i in 1:100){
#  gene_rand_ind <- sample(length(gene_tot_rand), size = 85, replace = T) ##size based on clique filtering
  gene_rand_ind <- sample(length(gene_tot_rand), size = 100, replace = T)
  gene_rand[[i]] <- gene_tot[gene_rand_ind]
}

##function for report degree of specific vertex of a subgraph
can_net <- read.delim("~/VDLab_scripts/BioGrid/biogrid_db_all_subnet.sif", header = T, sep = " ", stringsAsFactor = F)
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

##function for gene specific empirical p-value
Exome_pc123_srt_SKAT_case_enr_nCH_ranked <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_pval01, 
                                                   function(x) cbind.data.frame(x, "rank_score" = -log10(x[,14]) * x[,6])) ##added intra cohort MAF filtered variant positions

Exome_pc123_srt_SKAT_case_enr_nCH_ranked <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_ranked, function(x)
  x[order(x[,17], decreasing = T),])
gender_PC123_cont_df <- Exome_pc123_srt_SKAT_case_enr_nCH_ranked[[4]]

##function for degree detection
library(igraph)
library(ggnet)
library(intergraph)
library(network)
library(org.Hs.eg.db)


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
  
  l1 <- lapply(ret_df[,2], function(x)as_ids(adjacent_vertices(uniongraph, v = x)[[1]]))
  l2 <- lapply(l1, function(x)paste(x, sep=",", collapse=","))
  ret_df$interactors <- unlist(l2)
  
  return(ret_df)
}

df_degree <- get_topSKAT_degree(gender_PC123_cont_df)
gender_PC123_cont_df$degree <- df_degree[match(gender_PC123_cont_df$gene, df_degree$gene), 1]
gender_PC123_cont_df$int <- df_degree[match(gender_PC123_cont_df$gene, df_degree$gene), 3]
gender_PC123_cont_df$degree <- ifelse(is.na(gender_PC123_cont_df$degree), 0, gender_PC123_cont_df$degree)
gender_PC123_cont_df$rank_PPI <- gender_PC123_cont_df$rank_score * log(gender_PC123_cont_df$degree + 2)
gender_PC123_cont_df <- gender_PC123_cont_df[order(gender_PC123_cont_df$rank_PPI, decreasing = T),]

##clique based filtering

top_genes <- as.character(gender_PC123_cont_df$gene)
can_net_graph <- igraph::graph.data.frame(can_net, directed = F)
can_net_graph1 <- igraph::simplify(can_net_graph, remove.loops=T, remove.multiple = T)
can_net1 <- igraph::as_data_frame(can_net_graph1, what = "edges")
prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% top_genes 
                               & as.character(can_net1$to) %in% top_genes, ])

uniongraph <- igraph::graph.data.frame(prot_np_all, directed = F)
te1 <- cliques(uniongraph, min=3) ##used for selection of top 84 genes; this number also used for bootstrapping
#te1 <- cliques(uniongraph, min=4) ##use for final list
##use te1 for network enrichment
cpx <- names(unlist(te1))
cpx_genes <- unique(cpx)
cpx_np_all <- unique(can_net1[as.character(can_net1$from) %in% cpx 
                              & as.character(can_net1$to) %in% cpx, ])
cpx_graph <- igraph::graph.data.frame(cpx_np_all, directed = F)

##draw graph
library(ggnet)
net_mat_t_list <- get.adjacency(cpx_graph, type=c("both"), attr=NULL, names=TRUE, sparse = FALSE)
net_new_t_list <- network(net_mat_t_list, directed = FALSE)
network.vertex.names(net_new_t_list) <- V(cpx_graph)$name
gnet <- ggnet2(net_new_t_list, alpha = 0.75, edge.alpha = 0.5, label = TRUE, label.size = 3,  mode = "kamadakawai") 

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
s1 <- emp_pval_gene("TP53") #works
# plot(hist(unlist(s1)))
# plot(ecdf(unlist(s1)))
# 1 - ecdf(unlist(s1))(6)

##use t-distribution for skewed distribution
# 
a <- gender_PC123_cont_df[gender_PC123_cont_df$gene %in% "POT1",]$degree
 s <- sd(unlist(s1))
 n <- 85
xbar <- mean(unlist(s1))
t <- (xbar-a)/(s/sqrt(n))
# 2*pt(-abs(t),df=n-1)
# ##one sided
# t.test(unlist(s1),mu=5,alternative="less")

##use pnorm
##Induce graph including random networks and enriched genes
##compute centrality/degree

