
## Jaccard similarity function
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
library(ggdendro)
library(arrangements)
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

cliq_jacc_sim <- make_matrix(cpx_list_filt)
#library(reshape2)
#raw_mat <- acast(cliq_jacc_sim, m1 ~ m2, value.var="ov_tot")

g <- graph.data.frame(cliq_jacc_sim[,c(1:3)], directed=FALSE)
raw_mat <- get.adjacency(g, attr="perc", sparse=FALSE)
hc <- hclust(dist(raw_mat), "ave")
#plot(hc, hang = -1, cex = 0.8)
#groups <- cutree(hc, h=1.5)

#cliq_clusters <- cbind("cliq_name"= rownames(raw_mat),groups)
groups <- cutree(hc, h=1.5)
plot(hc, hang = -1, cex = 0.8)
gg <- rect.hclust(hc,k=length(unique(groups)))
dev.off()
clust.gr<-data.frame(num=unlist(gg), clust=rep(c(paste0("Clust_",1:length(unique(groups)))),times=sapply(gg,length)))
dendr <- dendro_data(hc, type="rectangle") 
text.df<-merge(label(dendr),clust.gr,by.x="label",by.y="row.names")
#print("plot clusters")
merge_clust <- ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=text.df, aes(x=x, y=y, label=label, hjust=0,color=clust), size=3) +
  geom_hline(yintercept=1.5, linetype="dashed") + 
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

merge_clust