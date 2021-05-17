library(ggplot2)
library(dplyr)
library(ggrepel)
library(knitr)
library(igraph)
library(ggnet)
library(intergraph)
library(network)
library(diffuStats)

ppi_res_fil_final_comb <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/isks_cyto_ppifisher_ASRB_comb_genes_top_new_rand_comp_score_Aug31.tsv", sep = "\t", header = T,
                                     stringsAsFactors = F)

##filter1
filt1_df <- ppi_res_fil_final_comb[ppi_res_fil_final_comb$ISKS > 1 & ppi_res_fil_final_comb$MGRB <= 14,]

##filter2: based on pvalue of ASRB genes
filt2_df <- filt1_df[filt1_df$pval_SKATbin.1 > 0.1 | is.na(filt1_df$pval_SKATbin.1),]
top_genes <- as.character(filt2_df$gene)
##filter3 : based on hypergeometric score of PPI (can be used for further fine-tuning and ranking)
##not needed for gene selection


strindb_biog_graph1 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/strindb_biog_graph_cyto.rds")

can_net1 <- igraph::as_data_frame(strindb_biog_graph1, what = "edges")
prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% top_genes 
                               & as.character(can_net1$to) %in% top_genes, ])
# prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% cep_genes 
#                                & as.character(can_net1$to) %in% cep_genes, ])

uniongraph <- igraph::graph.data.frame(prot_np_all, directed = F)

#score_val <- filt2_df[match(V(uniongraph)$name, filt2_df$gene), 6]
#score_val <- filt2_df[match(V(uniongraph)$name, filt2_df$gene), 44]
score_val <- filt2_df[match(V(uniongraph)$name, filt2_df$gene), 21]
#set_vertex_attr(uniongraph, input_vec, index = V(uniongraph), score_val)

#uniongraph1 <- uniongraph %>%
#  set_vertex_attr("input_vec", value = score_val)
uniongraph1 <- uniongraph %>%
  set_graph_attr("i_vec", value = score_val)

score_vec <- uniongraph1$i_vec
names(score_vec) <- V(uniongraph1)$name
#data("graph_toy")

op_vec <- diffuStats::diffuse(
  graph = uniongraph1, 
  method = "raw", scores = score_vec)

# K_diff <- diffusionKernel(uniongraph1)
# 
# op_vec <- diffuStats::diffuse(
#   graph = uniongraph1, 
#   method = "raw", scores = score_vec)

nodes_top <- head(names(op_vec[order(op_vec, decreasing = T)]),90)
nodes_top_scores <- op_vec[names(op_vec) %in% nodes_top]
#nodes_top_scores[order(nodes_top_scores, decreasing = T)]
plot(subgraph(uniongraph1, nodes_top))

igraph::plot.igraph(
  subgraph(uniongraph1, nodes_top), 
  vertex.color = diffuStats::scores2colours(nodes_top_scores),
  main = "Diffusion scores in our lattice"
)

op_vec <- diffuStats::diffuse(
  graph = uniongraph1, 
  method = "raw", scores = K_diff)

igraph::plot.igraph(
  uniongraph1, 
  vertex.color = diffuStats::scores2colours(op_vec),
  vertex.shape = diffuStats::scores2shapes(K_diff),
  main = "Diffusion scores in our lattice"
)
head(output_vec, 15)
