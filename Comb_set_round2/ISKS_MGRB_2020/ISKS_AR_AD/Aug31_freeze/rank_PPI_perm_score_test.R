##gene level combined p-value from rank statistics
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)

ISKS_MGRB_2020_rnd3 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv",
                                  sep = "\t", header = T, stringsAsFactors = F)

ppi_res_fil_final <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/ppi_res_fil_final_SKATbin_comb_str_biog_Aug31.tsv", 
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

strindb_graph1 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/strindb_graph.rds")
strindb_biog_graph1 <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/strindb_biog_graph1.rds")

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
ppi_res_fil_final[ppi_res_fil_final$degree_str == 0,]$PPI_p_val_wt_str <- sample(seq(0.6, 0.9, 0.0001), size = 268, replace = F)
ppi_res_fil_final[ppi_res_fil_final$degree_str_biog == 0,]$PPI_p_val_wt_str_biog <- sample(seq(0.6, 0.9, 0.0001), size = 170, replace = F)
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

ppi_res_fil_final[,44:49] <- res20_df_comb[match(ppi_res_fil_final$gene, res20_df_comb$gene),c(4:5,10:11,16:17)]

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

# inp_pval1 <- inp_pval
# inp_pval1$pval_SKATbin_wt <- inp_pval1$pval_SKATbin/ppi_res_fil_final$wt_diff
# inp_pval1 <- inp_pval1[,-1]
# ppi_res_fil_final$fisher_rank_score2 <- apply(inp_pval1, 1, function(x)sum(-log10(x)))

#write.table(ppi_res_fil_final[,c(1:12,14,17:20,22,26:28)], "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/scrambled_null_ppi_res_fil_final.tsv",
#            sep = "\t", row.names = F, quote = F)
write.table(ppi_res_fil_final[,c(1:12,14,17:20,22,26:28,32:34,38:58)], "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/scrambled_null_ppi_fisher_res_fil_final_Aug31.tsv",
            sep = "\t", row.names = F, quote = F)


#################################

###Analysis and QC : make a new script with the following section(s)

####cliques  in g1
top_genes <- g1
can_net1 <- igraph::as_data_frame(strindb_biog_graph1, what = "edges")
prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% top_genes 
                               & as.character(can_net1$to) %in% top_genes, ])

uniongraph <- igraph::graph.data.frame(prot_np_all, directed = F)

te1 <- cliques(uniongraph, min=3)
cpx <- names(unlist(te1))
cpx_np_all <- unique(can_net1[as.character(can_net1$from) %in% cpx 
                              & as.character(can_net1$to) %in% cpx, ])
cpx_graph <- igraph::graph.data.frame(cpx_np_all, directed = F)
cpx_graph <- igraph::simplify(cpx_graph, remove.loops=T, remove.multiple = T)
net_mat_t_list <- get.adjacency(cpx_graph, type=c("both"), attr=NULL, names=TRUE, sparse = FALSE)
net_new_t_list <- network(net_mat_t_list, directed = FALSE)
network.vertex.names(net_new_t_list) <- V(cpx_graph)$name
#gnet <- ggnet2(net_new_t_list, alpha = 0.75, edge.alpha = 0.5, label = TRUE, label.size = 3,  mode = #"kamadakawai") + theme(legend.position = "bottom") + theme(legend.title = element_blank())

gnet <- ggnet2(net_new_t_list, alpha = 0.75, edge.alpha = 0.5, color = "green",
               label = TRUE, label.size = 3, size = 7, mode = "fruchtermanreingold")
gnet

######rank sorted plots
ppi_res_fil_final <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/scrambled_null_ppi_fisher_res_fil_final.tsv",
                                sep = "\t", header = T, stringsAsFactors = F)


library(ggplot2)
library(ggrepel)
library(ggpubr)
scores_ppi_srt <- ppi_res_fil_final[order(ppi_res_fil_final$new_comp_score, decreasing = T),]
scores_ppi_srt$gene <- factor(scores_ppi_srt$gene, levels = scores_ppi_srt$gene)


ggscatter(scores_ppi_srt, x= "gene", y= "new_comp_score", label=ifelse(scores_ppi_srt$gene %in% g1,as.character(scores_ppi_srt$gene),""),
         repel = TRUE, size = 1, color = ifelse(scores_ppi_srt$gene %in% g1, "red", "grey50")) + theme(axis.title.x=element_blank(),
                                         axis.text.x=element_blank(),
                                         axis.ticks.x=element_blank()) + geom_hline(yintercept=17.26, linetype="dashed", 
                                                                                    color = "blue")

scores_ppi_rank <- ppi_res_fil_final[order(ppi_res_fil_final$rank_PPI, decreasing = T),]
scores_ppi_rank$gene <- factor(scores_ppi_rank$gene, levels = scores_ppi_rank$gene)

ggscatter(scores_ppi_rank, x= "gene", y= "new_comp_score", label=ifelse(scores_ppi_rank$gene %in% g1,as.character(scores_ppi_rank$gene),""),
          repel = TRUE, size = 1, color = ifelse(scores_ppi_rank$gene %in% g1, "red", "grey50")) + 
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + geom_hline(yintercept=10, linetype="dashed", color = "blue")


scores_ppi_rank <- ppi_res_fil_final[order(ppi_res_fil_final$fisher_rank_score1, decreasing = T),]
scores_ppi_rank$gene <- factor(scores_ppi_rank$gene, levels = scores_ppi_rank$gene)

ggscatter(scores_ppi_rank, x= "gene", y= "rank_PPI", label=ifelse(scores_ppi_rank$gene %in% g1,as.character(scores_ppi_rank$gene),""),
          repel = TRUE, size = 1, color = ifelse(scores_ppi_rank$gene %in% g1, "red", "grey50")) + 
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + geom_hline(yintercept=10, linetype="dashed", color = "blue")

scores_fisher_rank1 <- ppi_res_fil_final[order(ppi_res_fil_final$fisher_rank_score1, decreasing = T),]
scores_fisher_rank1$gene <- factor(scores_fisher_rank1$gene, levels = scores_fisher_rank1$gene)

ggscatter(scores_fisher_rank1, x= "gene", y= "fisher_rank_score1", label=ifelse(scores_fisher_rank1$gene %in% g1,as.character(scores_fisher_rank1$gene),""),
          repel = TRUE, size = 1, color = ifelse(scores_fisher_rank1$gene %in% g1, "red", "grey50")) + 
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + geom_hline(yintercept=7, linetype="dashed", color = "blue")

scores_fisher_rank2 <- ppi_res_fil_final[order(ppi_res_fil_final$fisher_rank_score2, decreasing = T),]
scores_fisher_rank2$gene <- factor(scores_fisher_rank2$gene, levels = scores_fisher_rank2$gene)

ggscatter(scores_fisher_rank2, x= "gene", y= "fisher_rank_score2", label=ifelse(scores_fisher_rank2$gene %in% g1,as.character(scores_fisher_rank2$gene),""),
          repel = TRUE, size = 1, color = ifelse(scores_fisher_rank2$gene %in% g1, "red", "grey50")) + 
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + geom_hline(yintercept=8.53, linetype="dashed", color = "blue")


###scatter plots for genes using only Biogrid network

ppi_res_fil_final_comb <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/isks_ppifisher_ASRB_comb_genes_top_new_comp_score_RankScores1_RankScore2.tsv",
                                sep = "\t", header = T, stringsAsFactors = F)

ppi_res_fil_final_comb <- ppi_res_fil_final_comb[ppi_res_fil_final_comb$gene %in% ppi_res_fil_final$gene,]

ppi_res_fil_final_comb$RankScore1 <- as.numeric(ppi_res_fil_final_comb$RankScore1)
scores_rank1 <- ppi_res_fil_final_comb[order(ppi_res_fil_final_comb$RankScore1, decreasing = T),]
scores_rank1$gene <- factor(scores_rank1$gene, levels = scores_rank1$gene)


ggscatter(scores_rank1, x= "gene", y= "RankScore1", label=ifelse(scores_rank1$gene %in% g1,as.character(scores_rank1$gene),""),
          repel = TRUE, size = 1, color = ifelse(scores_rank1$gene %in% g1, "red", "grey50")) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  geom_hline(yintercept=mean(scores_rank1$RankScore1), linetype="dashed", color = "blue")

ppi_res_fil_final_comb$RankScore2 <- as.numeric(ppi_res_fil_final_comb$RankScore2)
scores_rank2 <- ppi_res_fil_final_comb[order(ppi_res_fil_final_comb$RankScore2, decreasing = T),]
scores_rank2$gene <- factor(scores_rank2$gene, levels = scores_rank2$gene)


ggscatter(scores_rank2, x= "gene", y= "RankScore2", label=ifelse(scores_rank2$gene %in% g1,as.character(scores_rank2$gene),""),
          repel = TRUE, size = 1, color = ifelse(scores_rank1$gene %in% g1, "red", "grey50")) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  geom_hline(yintercept=mean(scores_rank2$RankScore2), linetype="dashed", color = "blue")

########correlation between pval_SKATbin and PPI_p_val_wt_str_biog

df_score_merge <- ppi_res_fil_final[,c(13,24,23)]
df_score_merge[df_score_merge$degree_str_biog == 0,]$PPI_p_val_wt_str_biog <- sample(seq(0.6, 0.9, 0.0001), size = 169, replace = F)
df_1 <- apply(df_score_merge[,c(1:2)], 2, function(x)(-1)*log10(x))
df_1 <- as.data.frame(df_1)
#df_2 <- df_1[df_1$PPI_p_val_wt_str_biog <40,]
sp <- ggscatter(df_1, x = "PPI_p_val_wt_str_biog", y = "pval_SKATbin",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)

# Add correlation coefficient
p1 <- sp + stat_cor(method = "spearman", label.x = 40, label.y = 6)

########correlation between pval_SKATbin and wt_diff

df_score_merge_wt_diff <- ppi_res_fil_final[,c(6,13,24,23)]
df_2 <- df_score_merge_wt_diff[,c(1:2)]
df_2$pval_SKATbin = -log10(df_1$pval_SKATbin)
#df_2 <- df_1[df_1$PPI_p_val_wt_str_biog <40,]
sp2 <- ggscatter(df_2, x = "wt_diff", y = "pval_SKATbin",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
 p2 <- sp2 + stat_cor(method = "spearman", label.x = 250, label.y = 6)

 source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")

 multiplot(p1, p2, cols = 2) 
 
 ##perform PCA and extract loading using PPI_p_val_wt_str_biog, pval_SKATbin and wt_diff


 library(factoextra)
 df_12 <- cbind.data.frame(df_1, "wt_diff" = df_2[,1])
 res.pca <- prcomp(df_12, scale = TRUE)
 fviz_eig(res.pca)

 
 # Eigenvalues
 eig.val <- get_eigenvalue(res.pca)
 eig.val
 
 # Results for Variables
 res.var <- get_pca_var(res.pca)
 res.var$coord          # Coordinates
 res.var$contrib        # Contributions to the PCs
 res.var$cos2           # Quality of representation 
 # Results for individual genes
 res.ind <- get_pca_ind(res.pca)
 res.ind$coord          # Coordinates
 res.ind$contrib        # Contributions to the PCs
 res.ind$cos2           # Quality of representation 
 
 ##score using res.ind$contrib
 score_mat <- res.ind$contrib
 rownames(score_mat) <- ppi_res_fil_final$gene
 pca_scores_gene <- rowSums(score_mat[,1:2]) ##add contributions from dim1 and dim2
 
 ppi_res_fil_final$pca_scores_gene <- pca_scores_gene
 
 
 pca_scores <- ppi_res_fil_final[order(ppi_res_fil_final$pca_scores_gene, decreasing = T),]
 pca_scores$gene <- factor(pca_scores$gene, levels = pca_scores$gene)
 
 ggscatter(pca_scores, x= "gene", y= "pca_scores_gene", label=ifelse(pca_scores$gene %in% g1,as.character(pca_scores$gene),""),
           repel = TRUE, size = 1, color = ifelse(pca_scores$gene %in% g1, "red", "grey50")) + 
   theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + geom_hline(yintercept=mean(ppi_res_fil_final$pca_scores_gene), linetype="dashed", color = "blue")
 

 
 
 # ggplot(scores_ppi_srt, aes(x= gene, y= new_comp_score, label=gene))+
#   geom_point(color = ifelse(scores_ppi_srt$gene %in% g1, "red", "grey50"), size = 1.5) +
#   geom_text(aes(label=ifelse(gene %in% g1,as.character(gene),'')),hjust=-0.5,vjust=1, angle = 45, 
#             size = 2.5, position=position_jitter(width=2,height=3)) + 
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# 
# ggplot(scores_ppi_srt, aes(x= gene, y= new_comp_score, label=gene))+
#   geom_point()+
#   geom_text_repel(data = scores_ppi_srt[scores_ppi_srt$gene %in% g1,], 
#                   aes(x= gene, y= new_comp_score, label=gene))
# 
# library(dplyr)
# ggplot(scores_ppi_srt, aes(x= gene, y= new_comp_score, label=gene)) +
#   geom_point() +
#   geom_text_repel(data = . %>% 
#                     mutate(label = ifelse(gene %in% g1,as.character(gene),"")),
#                   aes(label = label, color = "red"), 
#                   box.padding = 1, show.legend = FALSE) + theme(axis.title.x=element_blank(),
#                                                                 axis.text.x=element_blank(),
#                                                                 axis.ticks.x=element_blank()) + scale_size_manual(values = c(3.5))
# #  theme_bw()
# 
# ggplot(scores_ppi_srt, aes(x= gene, y= new_comp_score, label=gene)) +
#   geom_point(color = ifelse(scores_ppi_srt$gene %in% g1, "red", "grey50")) +
#   geom_text_repel(data = . %>% 
#                     mutate(label = ifelse(gene %in% g1,as.character(gene),"")),
#                   aes(label = label, color = "red"), 
#                   box.padding = 1, show.legend = FALSE, segment.size = 0.25, xlim = c(200, 1200)) + theme(axis.title.x=element_blank(),
#                                                                                                           axis.text.x=element_blank(),
#                                                                                                           axis.ticks.x=element_blank())


################(Not needed)
##Exhaustive network test
# ISKS_MGRB_2020_rnd3 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/comb_ISKS_MGRB_2020_C345_gnom_0002_all_genes_AD_ver4.tsv",
#                                   sep = "\t", header = T, stringsAsFactors = F)
# 
# gene_tot_rand <- as.character(ppi_res_fil_final$gene)
# ##degree randomisation (use entire gene list)
# gene_tot_rand_deg <- as.character(ISKS_MGRB_2020_rnd3$gene)
# 
# `%nin%` = Negate(`%in%`)
# 
# gene_tot_rand_deg <- gene_tot_rand_deg[gene_tot_rand_deg %nin% gene_tot_rand]
# 
# #1. for random_score randomisation: use the genes that have a positive rank_score 
# gene_rand <- list()
# for(i in 1:100){
#   #  gene_rand_ind <- sample(length(gene_tot_rand), size = 85, replace = T) ##size based on clique filtering
#   gene_rand_ind <- sample(length(gene_tot_rand), size = 100, replace = F)
#   gene_rand[[i]] <- gene_tot_rand[gene_rand_ind]
# }
# #2. for degree randomization
# deg_tol_rand <- list()
# for(i in 1:100){
#   #  gene_rand_ind <- sample(length(gene_tot_rand), size = 85, replace = T) ##size based on clique filtering
#   gene_rand_deg_ind <- sample(length(gene_tot_rand_deg), size = 1177, replace = F) ##size = 1177 number of genes in ppi_res_fil_final input
#   deg_tol_rand[[i]] <- gene_tot_rand_deg[gene_rand_deg_ind]
# }
# 
# ##generate degree from random network
# #can_net_graph <- igraph::graph.data.frame(can_net, directed = F)
# #can_net_graph1 <- igraph::simplify(can_net_graph, remove.loops=T, remove.multiple = T)
# get_deg_rand <- function(gene_list, gene_name, graph){
#  # gene_list <- c(gene_name, gene_list)
#   can_net1 <- igraph::as_data_frame(graph, what = "edges")
#   prot_np_all <- unique(can_net1[as.character(can_net1$from) %in% gene_list
#                                  & as.character(can_net1$to) %in% gene_list, ])
#   
#   if(dim(prot_np_all[prot_np_all$from %in% gene_name | prot_np_all$to %in% gene_name,])[1] != 0 ){
#     uniongraph <- igraph::graph.data.frame(prot_np_all, directed = F)
#     return(igraph::degree(uniongraph, gene_name, mode="all"))
#   }
#   else{
#     #    return(NULL)
#     return(0)
#   }
# }
# 
# ###Function for parallel processing
# library(doParallel)
# library(doMC)
# registerDoMC(30)
# 
# ##precomputing permuted degree with log transform (compute once and save)
# ##speeds up the computation enormously
# 
# perm_deg <- function(gene_sym, graph){
#   gene_add <- lapply(deg_tol_rand, function(x)c(gene_sym, x)) ##precompute
#   ##degree randomization
#   gene_add_degree <- lapply(gene_add, function(x)get_deg_rand(x, gene_sym, graph)) ##precompute
#   return(gene_add_degree)
# }
# 
# genes <- unique(as.character(ppi_res_fil_final$gene))
# #g1 <- genes[1:3]
# system.time(precomp_perm_net <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
# {perm_deg(genes[i], strindb_biog_graph1)})
# 
# precomp_perm_net1 <- lapply(precomp_perm_net, function(x)unlist(x))
# names(precomp_perm_net1) <- genes
# 
# ##test degree enrichment
# ##add options for rank scores: with/without wt_diff
# emp_pval_degree <- function(gene_sym,plot_hist = NULL){
#   ##add gene of interest to compute degree in randomly sampled genes
#   #  gene_add <- lapply(deg_tol_rand, function(x)c(gene_sym, x)) ##precompute
#   ##degree randomization
#   # gene_add_degree <- lapply(gene_add, function(x)log(get_deg_rand(x, gene_sym) + 2)) ##precompute
#   
#   gene_add_degree <- unlist(precomp_perm_net1[gene_sym])
#   
#   #plot(hist(unlist(gene_add_rand_rank_PPI_scores)))
#   if(plot_hist == 1){
#     # par(mfrow=c(1,3))
#     # hist(unlist(gene_add_degree))
#     #hist(unlist(gene_add_rand_rank_score))
#     # hist(unlist(gene_add_rand_rank_score))
#     hist(gene_add_degree)
#   }
#   #top_500 SKAT
#   orig_deg <- get_deg_rand(ppi_res_fil_final$gene, gene_sym, strindb_biog_graph1)
#   # orig_rank_PPI <- gender_PC123_cont_df_top[gender_PC123_cont_df_top$gene %in% gene_sym,]$rank_score * log(orig_deg + 2)
#   
#   a <- orig_deg
#   s <- sd(gene_add_degree, na.rm = T)
#   #n <- length(gene_rand[[1]])
#   #n <- length(unlist(null_dist))
#   n <- 100
#   #xbar <- mean(gene_add_rank_PPI_scores)
#   xbar <- mean(gene_add_degree)
#   t <- (xbar-a)/(s/sqrt(n))
#   ##one sided
#   # t_out <- t.test(gene_add_rank_PPI_scores,mu=a,alternative="less")
#   t_out <- t.test(gene_add_degree,mu=a,alternative="less")
#  
#     test_stat <- cbind.data.frame("orig_deg" = orig_deg, "mean_rand_rank_PPI" = xbar,
#                                   "t_stat" = t, 
#                                   "pval_t-test" = t_out$p.value)
#   return(test_stat)
# }

