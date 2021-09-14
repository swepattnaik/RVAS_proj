##R2Q8
##demonstrate by gridsearch that the cutoff of 0.0002 reduces the C4's in controls vs cases
##in cognate phenotypes (MPNIST)
##This cutoff mainly applies to C4/5 variants, gnomad_NFE == 0 is used to filter C3s
#1. extract variants from sarc genes (Extended list if needed) from raw file
#2. extract all C4's 

`%nin%` = Negate(`%in%`)

Ex_tab_mgrb <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/MGRB/all_mgrb_combset2020_variants_AR_AD_all_fields_clingene_rnd3.tsv", sep = "\t", header = T, stringsAsFactors = F)
Ex_tab_isks <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/all_isks_combset2020_variants_AR_AD_all_fields_clingene_rnd3_freeze.tsv", sep = "\t", header = T, stringsAsFactors = F)
Ex_tab_mgrb$set <- gsub("ISKS", "MGRB", Ex_tab_mgrb$set)

Ex_tab <- rbind.data.frame(Ex_tab_mgrb, Ex_tab_isks)
Ex_tab$gnomad_AF <- gsub("\\[", "", Ex_tab$gnomad_AF)
Ex_tab$gnomad_AF <- gsub("\\]", "", Ex_tab$gnomad_AF)
class(Ex_tab$gnomad_AF) <- "numeric"
Ex_tab$gnomad_AF_NFE <- gsub("\\[", "", Ex_tab$gnomad_AF_NFE)
Ex_tab$gnomad_AF_NFE <- gsub("\\]", "", Ex_tab$gnomad_AF_NFE)
class(Ex_tab$gnomad_AF_NFE) <- "numeric"
Ex_tab$gnomad_AF <- ifelse(is.na(Ex_tab$gnomad_AF), 0, Ex_tab$gnomad_AF)
Ex_tab$swegen_AF <- ifelse(is.na(Ex_tab$swegen_AF), 0, Ex_tab$swegen_AF)
Ex_tab$gnomad_AF_NFE <- ifelse(is.na(Ex_tab$gnomad_AF_NFE), 0, Ex_tab$gnomad_AF_NFE)
##GQ filter

Ex_tab <- Ex_tab[Ex_tab$GQ >= 80,]
Ex_tab <- Ex_tab[!is.na(Ex_tab$SAMPLE),]

##sarc genes
sarc_cont <- c("TP53", "NF1", "BRCA2", "ERCC2", "EXT1", "EXT2","SDHA", "SDHB", "SDHD")
Ex_tab_sarc <- Ex_tab[Ex_tab$gene_symbol %in% sarc_cont,]
write.table(Ex_tab_sarc, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/R2Q8/Ex_tab_sarc_genes_rnd3_freeze.tsv",
            sep = "\t", quote = F, row.names = F)

##start gridsearch
Ex_tab_sarc <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/review_doc/R2Q8/Ex_tab_sarc_genes_rnd3_freeze.tsv",
                          sep = "\t", header = T, stringsAsFactors = F)
#grid_scenario <- c(0, 0.0001, 0.0002, 0.0004, 0.0016, 0.0128)
grid_scenario <- numeric()
  for (i in 0:7){grid_scenario[i+1] = (0.0001 * 2^i)} #covers the entire range of rare variants
grid_scenario = c(0, grid_scenario)
Ex_tab_sarc_C45 <- Ex_tab_sarc[Ex_tab_sarc$auto_call %in% c("C4", "C5"),]
table(Ex_tab_sarc[Ex_tab_sarc$auto_call %in% c("C4", "C5"),]$set)
summary(Ex_tab_sarc_C45$gnomad_AF_NFE)
mgrb_isks_list <- list()
for(i in 1:length(grid_scenario)){
  df_mgrb_isks <- as.data.frame(table(Ex_tab_sarc_C45[Ex_tab_sarc_C45$gnomad_AF_NFE <= grid_scenario[i],]$set))
  df_mgrb_isks$gnomad_AF_NFE = grid_scenario[i]
  mgrb_isks_list[[i]] <- df_mgrb_isks
}

mgrb_isks_list_df <- do.call("rbind.data.frame", mgrb_isks_list)
#mgrb_isks_list_df$log_AF = -log2(mgrb_isks_list_df$gnomad_AF_NFE)
library(ggplot2)
gg1 <- ggplot(data=mgrb_isks_list_df,
       aes(x=gnomad_AF_NFE, y=Freq, colour=Var1)) + scale_x_continuous(trans='log2') + 
  geom_line() +
  theme(legend.position = "bottom")
##show difference in counts of ISKS - MGRB
diff_freq <- sapply(mgrb_isks_list, function(x)x[1,2] - x[2,2])
df_diff_freq <- cbind.data.frame("gnomad_AF_NFE" = grid_scenario, "diff"= diff_freq)
gg2 <- ggplot(data=df_diff_freq,
       aes(x=gnomad_AF_NFE, y=diff)) + scale_x_continuous(trans='log2') + geom_line()

source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")
multiplot(gg1, gg2, cols = 1)

##Fisher test
source("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/review_add/R2Q8_NFEcut/pos_con_diagnostic.R")
cpx_OR_fisher_sarc(40, 15, 1644, 3205, "Sarc_genes")
mgrb_isks_list_ft <- lapply(mgrb_isks_list, function(x)cpx_OR_fisher_sarc(x[1,2], x[2,2], 1644, 3205, "Sarc_genes"))
mgrb_isks_list_ft_df <- do.call("rbind.data.frame", mgrb_isks_list_ft)
mgrb_isks_list_ft_df$NFE <- grid_scenario

##final plot
library(ggplot2)

scaleFUN <- function(x) sprintf("%.5f", x)
  fp <- ggplot(data=mgrb_isks_list_ft_df, aes(x=NFE, y=log2(OR_Fish), ymin=log2(CI_lower), 
                                              ymax=log2(CI_upper))) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  scale_x_continuous(trans='log2', labels = scaleFUN) + 
  xlab("gnomAD_AF_NFE") + ylab("log(OR) Mean (95% CI)") +
  theme_bw() + # use a white background
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank())
  fp1 = fp + xlim(-0.00001, 0.013) + geom_vline(xintercept = 0.0002, lty = 2, color = "red")
  fp2 = fp + geom_vline(xintercept = 0.0002, lty = 2, color = "red") 
  #fp + theme(axis.text.x = element_text(vjust = 0.5)) + geom_hline(yintercept = 0.0002, lty = 2, color = "red")
  #fp + annotation_logticks()
  multiplot(fp1, fp2, cols = 2)
  
  
  #   scale_x_continuous(breaks = log2(c(0, 0.0001, 0.0002, 0.0004, 0.0016, 0.0128),
  #labels = c(1, 0.0001, 0.0002, 0.0004, 0.0016, 0.0128))) + 
