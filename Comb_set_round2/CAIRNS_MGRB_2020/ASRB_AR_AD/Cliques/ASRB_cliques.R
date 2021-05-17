---
title: "ASRB_cliques"
author: "Swetansu Pattnaik"
date: "15/04/2021"
output: html_document
---

```{r setup, include=FALSE}
`%nin%` = Negate(`%in%`)
knitr::opts_chunk$set(echo = TRUE)
```

##load libraries
```{r echo = FALSE, warning = FALSE, message = FALSE, figure.align = "center", figure.height = 10, figure.width = 12}
library(readxl)
library(ggplot2)
library(ggdendro)
library(dendextend)
library(dplyr)
library(ggrepel)
library(knitr)
library(igraph)
library(ggnet)
library(intergraph)
library(network)
library(org.Hs.eg.db)
library(topGO)
```


##Clique detection (cliq_OR_SKAT.R and SKATOPINT_package.Rmd; line:701-711)
##Select ASRB genes without controlling for age
```{r}
Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CAIRNS/round2/ASRB_AR_AD/SKAT/Exome_skat_para_result_CAIRNS_combset2020_uni_MAF_PC1234_ver4.rds")

df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})
Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CAIRNS/round2/ASRB_AR_AD/SKAT/Exome_para_pc123_SKAT_Enriched_CAIRNS_2020.rds")
##Add p-value
Exome_pc123_srt_SKAT_case_enr_nCH_pval <- Map(cbind.data.frame, Exome_pc123_srt_SKAT_case_enr_nCH, Exome_pc123_srt_SKAT)

cairns_mgrb_genes <- Exome_pc123_srt_SKAT_case_enr_nCH_pval[[4]]
cairns_mgrb_genes_top <- cairns_mgrb_genes[cairns_mgrb_genes$wt_diff > 0 & cairns_mgrb_genes$pval_SKATbin < 0.1,]

filt1_df <- cairns_mgrb_genes_top[cairns_mgrb_genes_top$CAIRNS > 1 & cairns_mgrb_genes_top$MGRB <= 14,]
top_genes <- as.character(filt1_df$gene)
##filter3 : based on hypergeometric score of PPI (can be used for further fine-tuning and ranking)
##not needed for gene selection


strindb_biog_graph1 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/strindb_biog_graph_cyto.rds")

can_net1 <- igraph::as_data_frame(strindb_biog_graph1, what = "edges")
prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% top_genes 
                               & as.character(can_net1$to) %in% top_genes, ])
# prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% cep_genes 
#                                & as.character(can_net1$to) %in% cep_genes, ])

uniongraph <- igraph::graph.data.frame(prot_np_all, directed = F)

te3 <- max_cliques(uniongraph, min=3)
names(te3) <- paste("c", 1:length(te3), sep = "_")

##OR use a more generic connected component based filter
groups <- components(uniongraph)
group_tab <- as.data.frame(groups$membership)
group_tab$genes <- names(groups$membership)
#saveRDS(te3, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Cliques/Cliques_max101.rds", compress = T)
saveRDS(te3, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/CAIRNS/round2/ASRB_AR_AD/SKAT/Cliques/Cliques_min4.rds", compress = T)

```

