mgrb_nc_var_filt4 = read.delim("~/RVAS/ISKS_MGRB_NC/results/mgrb_nc_var_filt_rm_dup_NC_scored.tsv", 
                               sep = "\t", header = T, stringsAsFactors = F)
#mgrb_nc_var_filt4$QC_id = apply(mgrb_nc_var_filt4[,3:5], 1 , paste, collapse = ":")
mgrb_nc_var_filt5 = mgrb_nc_var_filt4[mgrb_nc_var_filt4$NC_scores >= 95,] ##Use 95; similar to 95% CI
mgrb_nc_var_filt5 = mgrb_nc_var_filt5[!is.na(mgrb_nc_var_filt5$NC_scores),]
isks_nc_var_filt4 = read.delim("~/RVAS/ISKS_MGRB_NC/results/isks_nc_var_filt_rm_dup_NC_scored.tsv", 
                               sep = "\t", header = T, stringsAsFactors = F)
#isks_nc_var_filt4$QC_id = apply(isks_nc_var_filt4[,3:5], 1 , paste, collapse = ":")
isks_nc_var_filt5 = isks_nc_var_filt4[isks_nc_var_filt4$NC_scores >= 95,]
isks_nc_var_filt5 = isks_nc_var_filt5[!is.na(isks_nc_var_filt5$NC_scores),]
toMatch = c("synonymous_variant", "stop_gained", "missense_variant", "frameshift_variant",
            "splice_acceptor_variant", "splice_donor_variant")

isks_nc_var_filt6 = isks_nc_var_filt5[!grepl(paste(toMatch, collapse = "|"), 
                                             isks_nc_var_filt5$Consequence),]
mgrb_nc_var_filt6 = mgrb_nc_var_filt5[!grepl(paste(toMatch, collapse = "|"), 
                                             mgrb_nc_var_filt5$Consequence),]
##Gene check
##SKAT input
fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)

Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")

CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
CEP_HAUS_core_mod <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1")

MPNST_pos <- c("NF1", "LZTR1", "SDHA", "SDHB", "SDHC", "SDHD")


c1 = fil_tab[fil_tab$gene_symbol %in% "TP53",]$VARIANT
var_all1 = unlist(lapply(strsplit(c1, split=":"),function(x)paste(x[1],x[2],sep=":")))
s1 = fil_tab[fil_tab$gene_symbol %in% "TP53",]$SAMPLE
isks_nc_var_filt4[isks_nc_var_filt4$QC_id %in% c1,]$sid
mgrb_nc_var_filt4[mgrb_nc_var_filt4$QC_id %in% c1,]$sid

#########
inp = c(length(unique(isks_nc_var_filt6$sid)), 
            1644 - length(unique(isks_nc_var_filt6$sid)),
            length(unique(mgrb_nc_var_filt6$sid)),
            3205 - length(unique(mgrb_nc_var_filt6$sid)))
sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
colnames(sim_mat) <- c("case", "cont")
rownames(sim_mat) <- c("hits", "no_hits")
ft = fisher.test(sim_mat,conf.int = T, conf.level = 0.95)
ft_df <- cbind.data.frame("gene" = "Sarc_genes" ,"Cases" = sim_mat[1,1],
                          "Controls" = sim_mat[1,2],
                          "Fish_pval" = ft$p.value,"CI_lower" = ft$conf.int[1],
                          "CI_upper" = ft$conf.int[2],
                          "OR_Fish" = ft$estimate, "case_coh_size" = sum(sim_mat[,1]),
                          "control_coh_size" = sum(sim_mat[,2]),
                          "Coh" = "NC_ISKSvsMGRB")


##complex related enrichment
Shelterin <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "SMARCAL1", "STAG3", "TIMELESS")

CEP_HAUS_core <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1", "SSNA1")
CEP_HAUS_core_mod <- c("CEP63", "CEP72", "HAUS4", "HAUS5", "MZT1")

MPNST_pos <- c("NF1", "LZTR1", "SDHA", "SDHB", "SDHC", "SDHD")
HBOC_genes = c("BRCA1", "BRCA2", "PALB2")


##Function for fisher test in loop
fish_test <- function(isks_df, mgrb_df, coh_size_isks, coh_size_mgrb, genes){
  isks_df_sub = isks_df[isks_df$SYMBOL %in% genes, ]
  mgrb_df_sub = mgrb_df[mgrb_df$SYMBOL %in% genes, ]
  inp = c(length(unique(isks_df_sub$sid)), 
          coh_size_isks - length(unique(isks_df_sub$sid)),
          length(unique(mgrb_df_sub$sid)),
          coh_size_mgrb - length(unique(mgrb_df_sub$sid)))
  sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
  colnames(sim_mat) <- c("case", "cont")
  rownames(sim_mat) <- c("hits", "no_hits")
  ft = fisher.test(sim_mat,conf.int = T, conf.level = 0.95)
 # genes = deparse(substitute(genes))
  ft_df <- cbind.data.frame("gene" = genes ,"Cases" = sim_mat[1,1],
                            "Controls" = sim_mat[1,2],
                            "Fish_pval" = ft$p.value,"CI_lower" = ft$conf.int[1],
                            "CI_upper" = ft$conf.int[2],
                            "OR_Fish" = ft$estimate, "case_coh_size" = sum(sim_mat[,1]),
                            "control_coh_size" = sum(sim_mat[,2]),
                            "Coh" = "NC_ISKSvsMGRB")
  return(ft_df)
}



all_genes = unique(c(isks_nc_var_filt6$SYMBOL, mgrb_nc_var_filt6$SYMBOL))

genes_ft  = lapply(all_genes, function(x)fish_test(isks_nc_var_filt6, mgrb_nc_var_filt6,
                                                   1644, 3205, x))
genes_ft_df = do.call("rbind.data.frame", genes_ft)
genes_ft_df$adj.pval = p.adjust(genes_ft_df$Fish_pval, method = "bonferroni", n = length(genes_ft_df$Fish_pval))
write.table(genes_ft_df, "~/RVAS/ISKS_MGRB_NC/results/Fish_test_all_genes.tsv",
            sep = "\t", quote = F, row.names = F)
##Fisher test complex
fish_test_cpx <- function(isks_df, mgrb_df, coh_size_isks, coh_size_mgrb, cpx_name){
  isks_df_sub = isks_df[isks_df$SYMBOL %in% cpx_name, ]
  mgrb_df_sub = mgrb_df[mgrb_df$SYMBOL %in% cpx_name, ]
  inp = c(length(unique(isks_df_sub$sid)), 
          coh_size_isks - length(unique(isks_df_sub$sid)),
          length(unique(mgrb_df_sub$sid)),
          coh_size_mgrb - length(unique(mgrb_df_sub$sid)))
  sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
  colnames(sim_mat) <- c("case", "cont")
  rownames(sim_mat) <- c("hits", "no_hits")
  ft = fisher.test(sim_mat,conf.int = T, conf.level = 0.95)
  genes = deparse(substitute(cpx_name))
  ft_df <- cbind.data.frame("gene" = genes ,"Cases" = sim_mat[1,1],
                            "Controls" = sim_mat[1,2],
                            "Fish_pval" = ft$p.value,"CI_lower" = ft$conf.int[1],
                            "CI_upper" = ft$conf.int[2],
                            "OR_Fish" = ft$estimate, "case_coh_size" = sum(sim_mat[,1]),
                            "control_coh_size" = sum(sim_mat[,2]),
                            "Coh" = "NC_ISKSvsMGRB")
  return(ft_df)
}

genes_ft_Shel  = fish_test_cpx(isks_nc_var_filt6, mgrb_nc_var_filt6, 1644, 3205, Shelterin)
genes_ft_CEP  = fish_test_cpx(isks_nc_var_filt6, mgrb_nc_var_filt6, 1644, 3205, CEP_HAUS_core)
genes_ft_MPNST  = fish_test_cpx(isks_nc_var_filt6, mgrb_nc_var_filt6, 1644, 3205, MPNST_pos)
genes_ft_HBOC  = fish_test_cpx(isks_nc_var_filt6, mgrb_nc_var_filt6, 1644, 3205, HBOC_genes)

genes_ft_df_cpx = rbind.data.frame(genes_ft_Shel, genes_ft_CEP, genes_ft_MPNST, genes_ft_HBOC)
genes_ft_df_cpx$adj.pval = p.adjust(genes_ft_df_cpx$Fish_pval, method = "bonferroni", n = length(genes_ft_df_cpx$Fish_pval))
write.table(genes_ft_df_cpx, "~/RVAS/ISKS_MGRB_NC/results/Fish_test_complexes.tsv",
            sep = "\t", quote = F, row.names = F)
#write.table(genes_ft_df_cpx, "~/RVAS/ISKS_MGRB_NC/results/Fish_test_complexes_NC99.tsv",
#            sep = "\t", quote = F, row.names = F)

##Forest plot
library(ggplot2)
fp <- ggplot(data=genes_ft_df_cpx, aes(x=gene, y=OR_Fish, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("genesets") + ylab("Mean (95% CI)") + ggtitle("ISKS vs MGRB") + 
  theme_bw()
  
####
samp_mpnst = read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/MPNST_GIST/mpnst_samples.txt",
                        sep = "", stringsAsFactors = F)
mpnst_isks_nc_var_filt6 = isks_nc_var_filt6[isks_nc_var_filt6$rect_sam %in% samp_mpnst$x,]

mpnst_genes_ft_Shel  = fish_test_cpx(mpnst_isks_nc_var_filt6, mgrb_nc_var_filt6, 157, 3205, Shelterin)
mpnst_genes_ft_CEP  = fish_test_cpx(mpnst_isks_nc_var_filt6, mgrb_nc_var_filt6, 157, 3205, CEP_HAUS_core)
mpnst_genes_ft_MPNST  = fish_test_cpx(mpnst_isks_nc_var_filt6, mgrb_nc_var_filt6, 157, 3205, MPNST_pos)
mpnst_genes_ft_HBOC  = fish_test_cpx(mpnst_isks_nc_var_filt6, mgrb_nc_var_filt6, 157, 3205, HBOC_genes)

mpnst_genes_ft_df_cpx = rbind.data.frame(mpnst_genes_ft_Shel, mpnst_genes_ft_CEP, mpnst_genes_ft_MPNST, mpnst_genes_ft_HBOC)
mpnst_genes_ft_df_cpx$adj.pval = p.adjust(mpnst_genes_ft_df_cpx$Fish_pval, method = "bonferroni", n = length(mpnst_genes_ft_df_cpx$Fish_pval))
write.table(mpnst_genes_ft_df_cpx, "~/RVAS/ISKS_MGRB_NC/results/mpnst_Fish_test_complexes.tsv",
            sep = "\t", quote = F, row.names = F)
#write.table(mpnst_genes_ft_df_cpx, "~/RVAS/ISKS_MGRB_NC/results/mpnst_Fish_test_complexes_NC99.tsv",
#            sep = "\t", quote = F, row.names = F)

fp_mpnst <- ggplot(data=mpnst_genes_ft_df_cpx, aes(x=gene, y=OR_Fish, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("genesets") + ylab("Mean (95% CI)") + ggtitle("MPNST vs MGRB") + 
  theme_bw()

source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")

multiplot(fp, fp_mpnst, cols = 2)
####################check for all variants

c1 = fil_tab[fil_tab$gene_symbol %in% "TP53",]$VARIANT
c2 = unlist(lapply(strsplit(c1, split = ":"), function(x)paste(x[1], x[2], sep = ":")))

isks_nc_var <- readRDS("~/RVAS/ISKS_MGRB_NC/vep_annot_comb.rds")
#isks_nc_var$QC_id = apply(isks_nc_var[,3:5], 1 , paste, collapse = ":")
isks_nc_var$gnomAD_NFE_AF <- as.numeric(isks_nc_var$gnomAD_NFE_AF)
isks_nc_var$gnomAD_NFE_AF <- ifelse(is.na(isks_nc_var$gnomAD_NFE_AF), 0, isks_nc_var$gnomAD_NFE_AF)
colnames(isks_nc_var)[1:10] <- c("CHROM", "POS", "ID", "REF", "ALT", "GT", "GQ", "PL", "DP", "AD")
isks_nc_var$AD <- as.numeric(unlist(lapply(strsplit(isks_nc_var$AD, split = ","), function(x)x[2])))
isks_nc_var <- isks_nc_var[isks_nc_var$AD > 0, ]

isks_nc_var[isks_nc_var$QC_id %in% c1,]$QC_id
mgrb_nc_var_filt4[mgrb_nc_var_filt4$QC_id %in% c1,]$sid

dim(isks_nc_var[isks_nc_var$Location %in% c2,])
