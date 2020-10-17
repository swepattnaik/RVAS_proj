##Visualisation
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(knitr)
library(igraph)
library(ggnet)
library(intergraph)
library(network)
library(org.Hs.eg.db)
library(topGO)

`%nin%` = Negate(`%in%`)

cpx_list <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Cliques/Cliques_max101.rds")

cpx_genes <- unique(unlist(lapply(cpx_list, function(x)names(x))))
##filter list based on burden test (adjusted Pval < 0.01)
comb_burd_SKAT <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Cliques/Clique_OR_SKAT.tsv",
                             sep = "\t", header = T, stringsAsFactors = F)
sig_cliq <- comb_burd_SKAT[comb_burd_SKAT$adj_pval < 0.01,]$clique

cpx_list_filt <- cpx_list[sig_cliq]
cpx_genes_filt <- unique(unlist(lapply(cpx_list_filt, function(x)names(x))))

##graph

strindb_biog_graph1 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/strindb_biog_graph_cyto.rds")

can_net1 <- igraph::as_data_frame(strindb_biog_graph1, what = "edges")
prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% cpx_genes_filt 
                               & as.character(can_net1$to) %in% cpx_genes_filt, ])
save_sif <- prot_np_all
colnames(save_sif) <- c("node1", "node2")
write.table(save_sif, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Cliques/top_cliques.sif",
            row.names = F, sep = " ", quote = F)

cpx_graph <- igraph::graph.data.frame(prot_np_all, directed = F)

net_mat_t_list <- get.adjacency(cpx_graph, type=c("both"), attr=NULL, names=TRUE, sparse = FALSE)
net_new_t_list <- network(net_mat_t_list, directed = FALSE)
network.vertex.names(net_new_t_list) <- V(cpx_graph)$name

gnet <- ggnet2(net_new_t_list, alpha = 0.75, edge.alpha = 0.5, label = TRUE, label.size = 3,  mode = "kamadakawai")
gnet

net_new_t_list %v% "mod" <- as.character(comb_genes_df_new[match(network.vertex.names(net_new_t_list), as.character(comb_genes_df_new$genes)),2])
set.seed(312016)

ggnet2(net_new_t_list, color = "mod",
       size = 6, palette = "Set1", edge.alpha = 0.25, label = TRUE, label.size = 3,
       mode = "fruchtermanreingold",
       edge.color = c("color", "grey50"),
       color.legend = "modules") +
  theme(legend.position = "bottom")


plot_func_cliq <- function(cliq_list){
  cpx_genes <- names(cliq_list)
  can_net1 <- igraph::as_data_frame(strindb_biog_graph1, what = "edges")
  prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% cpx_genes 
                                 & as.character(can_net1$to) %in% cpx_genes, ])
  
  cpx_graph <- igraph::graph.data.frame(prot_np_all, directed = F)
  
  net_mat_t_list <- get.adjacency(cpx_graph, type=c("both"), attr=NULL, names=TRUE, sparse = FALSE)
  net_new_t_list <- network(net_mat_t_list, directed = FALSE)
  network.vertex.names(net_new_t_list) <- V(cpx_graph)$name
  
  cpx_plt <- ggnet2(net_new_t_list, alpha = 0.75, edge.alpha = 0.5, label = TRUE, label.size = 3,  mode = "kamadakawai")
return(cpx_plt)
}

all_cliq_plot <- lapply(cpx_list_filt, function(x)plot_func_cliq(x))
##multiplot
source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")
multiplot(plotlist = all_cliq_plot, cols = 3)

cpx_list_filt_gt4 <- cpx_list_filt[names(cpx_list_filt[which(sapply(cpx_list_filt, function(x)length(x) > 4))])]
all_cliq_plot_gt4 <- lapply(cpx_list_filt_gt4, function(x)plot_func_cliq(x))
multiplot(plotlist = all_cliq_plot_gt4, cols = 4)
