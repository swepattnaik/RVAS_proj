
`%nin%` = Negate(`%in%`)
Exome_skat <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_skat_para_result_isks_combset2020_uni_MAF_PC1234_ver4_clinrect_Aug31.rds")
df_skat <- lapply(Exome_skat, function(x)do.call("rbind.data.frame", x))
Exome_pc123_srt_SKAT <- lapply(df_skat, function(x) x[order(x$pval_SKATbin, decreasing = F),])
Exome_pc123_srt_SKAT <- lapply(Exome_pc123_srt_SKAT, function(x){colnames(x)[1] <- c("symbol"); return(x)})
Exome_pc123_srt_SKAT_case_enr_nCH <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/Exome_para_pc1234_SKAT_Enriched_ISKS_2020_Aug31.rds")
##Add p-value
Exome_pc123_srt_SKAT_case_enr_nCH_pval <- Map(cbind.data.frame, Exome_pc123_srt_SKAT_case_enr_nCH, Exome_pc123_srt_SKAT)

# Exome_pc123_srt_SKAT_case_enr_nCH_pval <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_pval, function(x) x[x$pval_SKATbin < 0.1 & x$wt_diff > 0,])
# 
# ##Better option for Enriched list: sort by -log(p_val) * weight difference ; appropriately factors influence of variant on ISKS phenotype
# ##use SKATO p-value; go back to using SKATbin (x[,14])
Exome_pc123_srt_SKAT_case_enr_nCH_ranked <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_pval,
                                                   function(x) cbind.data.frame(x, "rank_score" = -log10(x[,14])* x[,6]))
Exome_pc123_srt_SKAT_case_enr_nCH_ranked <- lapply(Exome_pc123_srt_SKAT_case_enr_nCH_ranked, function(x)
  x[order(x[,17], decreasing = F),])

Exome_pc123_srt_SKAT_case_enr_nCH_ranked_all = Exome_pc123_srt_SKAT_case_enr_nCH_ranked[[4]] 

ranked_list = Exome_pc123_srt_SKAT_case_enr_nCH_ranked_all$rank_score
names(ranked_list) = Exome_pc123_srt_SKAT_case_enr_nCH_ranked_all$gene

##Cliques/Genesets
library(igraph)
library(fgsea)
#data(examplePathways)
#data(exampleRanks)
cpx_list <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Cliques/Cliques_max64.rds")
comb_burd_SKAT <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/Aug31/Cliques/Clique_OR_SKAT_968.tsv",
                             sep = "\t", header = T, stringsAsFactors = F)
sig_cliq <- comb_burd_SKAT[comb_burd_SKAT$adj_pval < 0.01,]$clique

cpx_list <- cpx_list[sig_cliq]
cpx_list <- lapply(cpx_list, function(x)names(x))
#cpx_genes <- unique(unlist(lapply(cpx_list, function(x)names(x))))

fgseaRes <- fgsea(pathways = cpx_list, 
                  stats    = ranked_list,
                  minSize  = 3,
                  maxSize  = 500, nperm = 1000, nproc = 8)
View(fgseaRes)


##merge cliques
##obtained from script: cliques_for_paper.Rmd (lines 119-178)

## Jaccard similarity function
library(arrangements)
library(ggdendro)
library(ggplot2)
make_matrix <- function(cliq_list) {
  
  # comb <- combinations(names(cliq_list), 2, replace = F)
  comb <- combinations(names(cliq_list), 2, replace = T) #include diagonals
  prct <- list()
  a_len <- list()
  b_len <- list()
  a_len_ov <- list()
  b_len_ov <- list()
  ov_tot <- list()
  for (i in 1:nrow(comb)) {
    #  a <- unlist(strsplit(as.character(exemp_df_L1[exemp_df_L1$L1 == comb[i, 1],]$mods), ","))
    #  b <- unlist(strsplit(as.character(exemp_df_L1[exemp_df_L1$L1 == comb[i, 2],]$mods), ","))
    a <- cliq_list[[comb[i,1]]]
    b <- cliq_list[[comb[i,2]]]
    prct[[i]] <- 2 * length(intersect(a, b)) / (length(a) + length(b))
    a_len[[i]] <- length(a) 
    b_len[[i]] <- length(b)
    a_len_ov[[i]] <- length(a[a %in% b])/length(a)
    b_len_ov[[i]] <- length(b[b %in% a])/length(b)
    ov_tot[[i]] <- length(intersect(a, b))
    # cat("\nMatching between", comb[i, 1], "and", comb[i, 2], "is", prct)
  }
  
  comb1 <- cbind.data.frame(comb, "perc" = unlist(prct), "c1" = unlist(a_len), "c2" = unlist(b_len),
                            "c1_perc_ov" = unlist(a_len_ov), "c2_perc_ov" = unlist(b_len_ov), 
                            "ov_tot" = unlist(ov_tot), stringsAsFactors=FALSE)
  colnames(comb1)[1] <- "m1"
  colnames(comb1)[2] <- "m2"
  
  return(comb1)
  
}

cliq_jacc_sim <- make_matrix(cpx_list)
#library(reshape2)
#raw_mat <- acast(cliq_jacc_sim, m1 ~ m2, value.var="ov_tot")

g <- graph.data.frame(cliq_jacc_sim[,c(1:3)], directed=FALSE)
raw_mat <- get.adjacency(g, attr="perc", sparse=FALSE)
##function to extract cluster after specifying a cutpoint
cluster_mem <- function(mat_inp){
  hc <- hclust(dist(mat_inp), "ave")
  #plot(hc, hang = -1, cex = 0.8)
  #abline(h = 1.5, lty = 2)
  groups <- cutree(hc, h=1.5)
  
  #cliq_clusters <- cbind("cliq_name"= rownames(raw_mat),groups)
  
  ###start here
  groups <- cutree(hc, h=1.5)
  plot(hc, hang = -1, cex = 0.8)
  gg <- rect.hclust(hc,k=length(unique(groups)))
  return(gg)
}
gg1 <- cluster_mem(raw_mat)
collapsed_cliques <- lapply(gg1, function(x)unique(unlist(lapply(cpx_list[x], function(y)y))))
names(collapsed_cliques) <- paste("Clust", 1:6, sep = "_")
#collapsed_cliques[[7]] = c("BRCA1", "BRCA2", "PALB2", "ATM", "FANCD2", "FANCI") #negative control
#collapsed_cliques[[7]] = c("BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2", "CDH1")
#collapsed_cliques[[7]] = c("CRYGA", "PTP4A3", "PAICS", "DNAH1", "RBP3", "NAA25")
set.seed(94032)
#collapsed_cliques[[7]] = names(sample(ranked_list[which(ranked_list > 0)], size = 6, replace = F))
tt = Exome_pc123_srt_SKAT_case_enr_nCH_ranked_all
tt$zscore = (tt$rank_score - mean(tt$rank_score))/sd(tt$rank_score)
collapsed_cliques[[7]] = sample(as.character(tt[tt$wt_diff > 0 & tt$pval_SKATbin < 0.1,]$gene), size = 6, replace = F)
names(collapsed_cliques)[7] = "Clust_7"

##test rank list (wt_diff > 0)
tt1_rank = tt$zscore
names(tt1_rank) = as.character(tt$gene)
#https://www.biostars.org/p/429563/ ##GSEA usage
set.seed(4532)
fgseaRes_col <- fgsea(pathways = collapsed_cliques, 
                  stats    = tt1_rank,
                  minSize  = 4,
                  maxSize  = 500, nperm = 10000, nproc = 8)
View(fgseaRes_col)
