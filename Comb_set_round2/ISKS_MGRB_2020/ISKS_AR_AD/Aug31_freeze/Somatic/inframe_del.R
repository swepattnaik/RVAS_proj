##inframe indels greater than size 15nt

`%nin%` = Negate(`%in%`)

Ex_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/all_shards_ISKS_combset2020_gnomadAF05.tsv", sep = "\t",
                     header = T, stringsAsFactors = F)
Ex_tab <- Ex_tab[!is.na(Ex_tab$SAMPLE),]
#remove variants ending with Asterisk; these donot have vep annotation
Ex_tab <- Ex_tab[!is.na(Ex_tab$gene_symbol),]

Ex_tab <- Ex_tab[!is.na(Ex_tab$DP),] ##There are records with NA in DP attribute ; file: sept05_rect

Ex_tab <- Ex_tab[!(Ex_tab$VAF %in% "NaN"),]
class(Ex_tab$VAF) <- "numeric" ##remove headers from shards in the next run
Ex_tab <- Ex_tab[!is.na(Ex_tab$VAF),] ##some VAFs are NAs in sept05_rect file
#Ex_tab1 <- Ex_tab[as.numeric(Ex_tab$VAF) >= 0.1 & as.numeric(Ex_tab$VAF) <= 0.35 & as.numeric(Ex_tab$DP) >= 10,]
#Ex_tab1 <- Ex_tab[as.numeric(Ex_tab$VAF) >= 0.35 & as.numeric(Ex_tab$DP) >= 10,]
Ex_tab1 <- Ex_tab[as.numeric(Ex_tab$DP) >= 10,] ##include both germline and somatic
##remove duplicates; CR57(1961)
rect_sam_dat <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/ISKS_RISC_LIONS_final_freeze.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
Ex_tab1$rect_SAMPLE <- rect_sam_dat[match(Ex_tab1$SAMPLE, rect_sam_dat$JCInputID), 9]
Ex_tab1 <- Ex_tab1[!is.na(Ex_tab1$rect_SAMPLE),]
Ex_tab1 <- Ex_tab1[,-1]
Ex_tab1$SAMPLE <- Ex_tab1$rect_SAMPLE
Ex_tab1 <- Ex_tab1[,c(120,1:118)]
Ex_tab1_inframe_del <- Ex_tab1[grepl("^inframe",Ex_tab1$vep_consequence), ]
Ex_tab1_inframe_del_sp <-Ex_tab1_inframe_del[grepl("^inframe_deletion", Ex_tab1_inframe_del$vep_consequence),]

Ex_tab1_inframe_del_sp$del_len <- unlist(lapply(strsplit(Ex_tab1_inframe_del_sp$VARIANT, split = ":"),
                                      function(x)nchar(x[3])))
Ex_tab1_inframe_del_sel <- Ex_tab1_inframe_del_sp[Ex_tab1_inframe_del_sp$del_len >= 15,]
Ex_tab1_inframe_del_sel$gnomad_AF <- gsub("\\[", "", Ex_tab1_inframe_del_sel$gnomad_AF)
Ex_tab1_inframe_del_sel$gnomad_AF <- gsub("\\]", "", Ex_tab1_inframe_del_sel$gnomad_AF)
class(Ex_tab1_inframe_del_sel$gnomad_AF) <- "numeric"
Ex_tab1_inframe_del_sel$gnomad_AF_NFE <- gsub("\\[", "", Ex_tab1_inframe_del_sel$gnomad_AF_NFE)
Ex_tab1_inframe_del_sel$gnomad_AF_NFE <- gsub("\\]", "", Ex_tab1_inframe_del_sel$gnomad_AF_NFE)
class(Ex_tab1_inframe_del_sel$gnomad_AF_NFE) <- "numeric"
Ex_tab1_inframe_del_sel$gnomad_AF <- ifelse(is.na(Ex_tab1_inframe_del_sel$gnomad_AF), 0, Ex_tab1_inframe_del_sel$gnomad_AF)
Ex_tab1_inframe_del_sel$swegen_AF <- ifelse(is.na(Ex_tab1_inframe_del_sel$swegen_AF), 0, Ex_tab1_inframe_del_sel$swegen_AF)
Ex_tab1_inframe_del_sel$gnomad_AF_NFE <- ifelse(is.na(Ex_tab1_inframe_del_sel$gnomad_AF_NFE), 0, Ex_tab1_inframe_del_sel$gnomad_AF_NFE)

#Ex_tab1_inframe_del_sel[Ex_tab1_inframe_del_sel$gene_symbol %in% c("CEP57", "CEP192"),]

write.table(Ex_tab1_inframe_del_sel, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/ISKS_GL_som_inframe_del_gte15nt.tsv",
            sep = "\t", row.names = F, quote = F)

###MGRB

Ex_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/MGRB/all_shards_MGRB_combset2020_gnomadAF05.tsv", sep = "\t",
                     header = T, stringsAsFactors = F)
Ex_tab <- Ex_tab[!is.na(Ex_tab$SAMPLE),]
#remove variants ending with Asterisk; these donot have vep annotation
Ex_tab <- Ex_tab[!is.na(Ex_tab$gene_symbol),]

Ex_tab <- Ex_tab[!is.na(Ex_tab$DP),] ##There are records with NA in DP attribute ; file: sept05_rect

#remove duplicates
dup_samp <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/fin_samp_dup_drop.txt", sep = "", header = T,
                       stringsAsFactors = F)
`%nin%` = Negate(`%in%`)
Ex_tab <- Ex_tab[Ex_tab$SAMPLE %nin% dup_samp$x,]
Ex_tab <- Ex_tab[!(Ex_tab$VAF %in% "NaN"),]
class(Ex_tab$VAF) <- "numeric" ##remove headers from shards in the next run
Ex_tab <- Ex_tab[!is.na(Ex_tab$VAF),] ##some VAFs are NAs in sept05_rect file
#Ex_tab1 <- Ex_tab[as.numeric(Ex_tab$VAF) >= 0.35 & as.numeric(Ex_tab$DP) >= 10,]
Ex_tab1 <- Ex_tab[as.numeric(Ex_tab$DP) >= 10,] ##include both germline and somatic
Ex_tab1_inframe_del <- Ex_tab1[grepl("^inframe",Ex_tab1$vep_consequence), ]
Ex_tab1_inframe_del_sp <-Ex_tab1_inframe_del[grepl("^inframe_deletion", Ex_tab1_inframe_del$vep_consequence),]

Ex_tab1_inframe_del_sp$del_len <- unlist(lapply(strsplit(Ex_tab1_inframe_del_sp$VARIANT, split = ":"),
                                                function(x)nchar(x[3])))
Ex_tab1_inframe_del_sel <- Ex_tab1_inframe_del_sp[Ex_tab1_inframe_del_sp$del_len >= 15,]
Ex_tab1_inframe_del_sel$gnomad_AF <- gsub("\\[", "", Ex_tab1_inframe_del_sel$gnomad_AF)
Ex_tab1_inframe_del_sel$gnomad_AF <- gsub("\\]", "", Ex_tab1_inframe_del_sel$gnomad_AF)
class(Ex_tab1_inframe_del_sel$gnomad_AF) <- "numeric"
Ex_tab1_inframe_del_sel$gnomad_AF_NFE <- gsub("\\[", "", Ex_tab1_inframe_del_sel$gnomad_AF_NFE)
Ex_tab1_inframe_del_sel$gnomad_AF_NFE <- gsub("\\]", "", Ex_tab1_inframe_del_sel$gnomad_AF_NFE)
class(Ex_tab1_inframe_del_sel$gnomad_AF_NFE) <- "numeric"
Ex_tab1_inframe_del_sel$gnomad_AF <- ifelse(is.na(Ex_tab1_inframe_del_sel$gnomad_AF), 0, Ex_tab1_inframe_del_sel$gnomad_AF)
Ex_tab1_inframe_del_sel$swegen_AF <- ifelse(is.na(Ex_tab1_inframe_del_sel$swegen_AF), 0, Ex_tab1_inframe_del_sel$swegen_AF)
Ex_tab1_inframe_del_sel$gnomad_AF_NFE <- ifelse(is.na(Ex_tab1_inframe_del_sel$gnomad_AF_NFE), 0, Ex_tab1_inframe_del_sel$gnomad_AF_NFE)

#Ex_tab1_inframe_del_sel[Ex_tab1_inframe_del_sel$gene_symbol %in% c("CEP57", "CEP192"),]

write.table(Ex_tab1_inframe_del_sel, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/MGRB/MGRB_GL_som_inframe_del_gte15nt.tsv",
            sep = "\t", row.names = F, quote = F)

##Enrichment analysis
`%nin%` = Negate(`%in%`)
library(dplyr)
ISKS_inf_del <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/ISKS_GL_som_inframe_del_gte15nt.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
ISKS_inf_del$set <- "ISKS"
ISKS_inf_del$var_ind <- paste(ISKS_inf_del$SAMPLE, ISKS_inf_del$VARIANT, ISKS_inf_del$del_len, sep = "_")
MGRB_inf_del <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/MGRB/MGRB_GL_som_inframe_del_gte15nt.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
MGRB_inf_del$set <- "MGRB"
MGRB_inf_del$var_ind <- paste(MGRB_inf_del$SAMPLE, MGRB_inf_del$VARIANT, MGRB_inf_del$del_len, sep = "_")
##remove multiple(greater than 3) inframe deletions within sample within gene
ISKS_inf_del_mult <- ISKS_inf_del %>% group_by(SAMPLE, gene_symbol) %>% summarise(n = n())

MGRB_inf_del_mult <- MGRB_inf_del %>% group_by(SAMPLE, gene_symbol) %>% summarise(n = n())

ISKS_inf_del_mult <- as.data.frame(ISKS_inf_del_mult)

MGRB_inf_del_mult <- as.data.frame(MGRB_inf_del_mult)

ISKS_inf_del_mult <- ISKS_inf_del_mult[ISKS_inf_del_mult$n >= 2,]
MGRB_inf_del_mult <- MGRB_inf_del_mult[MGRB_inf_del_mult$n >= 2,]

##function to remove redundant inframe deletions

inf_trim <- function(mult_infdel_inp, inf_del_varfile){
  sub_inf_sam_gen_ind <- list()
  for(i in 1:dim(mult_infdel_inp)[1]){
    sub_inf_sam_gen <- inf_del_varfile[inf_del_varfile$gene_symbol %in% mult_infdel_inp$gene_symbol[i] &
                                         inf_del_varfile$SAMPLE %in% mult_infdel_inp$SAMPLE[i], ]
    sub_inf_sam_gen_ind[[i]] <- sub_inf_sam_gen[order(sub_inf_sam_gen$del_len, decreasing = T),]$var_ind[-1]
    
  }
  names(sub_inf_sam_gen_ind) <- paste(mult_infdel_inp$SAMPLE, mult_infdel_inp$gene_symbol, sep = "_")
  return(sub_inf_sam_gen_ind)
}

rm_ISKS_ind <- inf_trim(ISKS_inf_del_mult, ISKS_inf_del)
rm_ISKS_ind_vec <- unlist(rm_ISKS_ind)
ISKS_inf_del <- ISKS_inf_del[ISKS_inf_del$var_ind %nin% rm_ISKS_ind_vec,]

rm_MGRB_ind <- inf_trim(MGRB_inf_del_mult, MGRB_inf_del)
rm_MGRB_ind_vec <- unlist(rm_MGRB_ind)
MGRB_inf_del <- MGRB_inf_del[MGRB_inf_del$var_ind %nin% rm_MGRB_ind_vec,]

##EDA count MGRB and ISKS 
##frequency versus length of Inf_del (All som+GL)
library(ggplot2)
hist(ISKS_inf_del$del_len)
comb_inf_del_df <- rbind.data.frame(ISKS_inf_del, MGRB_inf_del)
p<- ggplot(comb_inf_del_df, aes(x=del_len, color=set)) + geom_density() +
  geom_vline(aes(xintercept=50), linetype="dashed")
p

##check germline
GL_comb_inf_del_df <-  comb_inf_del_df %>% filter(VAF >= 0.35)
p_GL<- ggplot(GL_comb_inf_del_df, aes(x=del_len, color=set)) + geom_density() +
  geom_vline(aes(xintercept=50), linetype="dashed") + ggtitle("Germline")

#p_GL

##check somatic
Som_comb_inf_del_df <-  comb_inf_del_df %>% filter(VAF >= 0.1 & VAF <= 0.35)
p_Som<- ggplot(Som_comb_inf_del_df, aes(x=del_len, color=set)) + geom_density() +
  geom_vline(aes(xintercept=50), linetype="dashed") + ggtitle("Somatic")

#p_Som

source("~/APCluster_graphs/case_TCGA_new_PRAD_3431/multi_plot.R")
multiplot(p_GL, p_Som, cols = 2)

##Allele count
##Germline analysis
library(reshape2)
GL_fun <- function(GL_var_df){
GL_comb_inf_del_VAR <- unique(GL_var_df$VARIANT)
GL_comb_inf_del_SAM <- unique(GL_var_df$SAMPLE)

match_var_GL <- function(var_inp){
  samp_id <- unique(GL_var_df[GL_var_df$VARIANT %in% var_inp, ]$SAMPLE)
  vec_sam <- ifelse(GL_comb_inf_del_SAM %in% samp_id, 1, 0)
  return(vec_sam)
}
#1:1577121:CCTTCCTCCTCCTCCT:C
#GL_comb_inf_del_SAM[which(match_var("1:1577121:CCTTCCTCCTCCTCCT:C") == 1)]

vec_list <- list()
vec_list <- lapply(GL_comb_inf_del_VAR, function(x)match_var_GL(x))

vec_list_mat <- do.call("rbind", vec_list)
colnames(vec_list_mat) <- GL_comb_inf_del_SAM
rownames(vec_list_mat) <- GL_comb_inf_del_VAR

##MAF 
Var_MAF_ISKS <- rowSums(vec_list_mat[,grep("^[ABZ]", colnames(vec_list_mat), invert = T)])
# Var_MAF_ISKS_df <- as.data.frame(Var_MAF_ISKS)
# Var_MAF_ISKS_df$set <- "ISKS"
# colnames(Var_MAF_ISKS_df)[1] <- "Allele_Count"
Var_MAF_MGRB <- rowSums(vec_list_mat[,grep("^[ABZ]", colnames(vec_list_mat))])
# Var_MAF_MGRB_df <- as.data.frame(Var_MAF_MGRB)
# Var_MAF_MGRB_df$set <- "MGRB"
# colnames(Var_MAF_MGRB_df)[1] <- "Allele_Count"
# Var_MAF_ISKS_MGRB_df <- rbind.data.frame(Var_MAF_ISKS_df, Var_MAF_MGRB_df)
Var_MAF_Tot <- rowSums(vec_list_mat)

Var_MAF_ISKS_MGRB_df <- as.data.frame(cbind(Var_MAF_ISKS, Var_MAF_MGRB, Var_MAF_Tot))
return(Var_MAF_ISKS_MGRB_df)
}

GL_var_MAF_ISKS_MGRB_df <- GL_fun(GL_comb_inf_del_df)
GL_var_MAF_ISKS_MGRB_df$VAR <- rownames(GL_var_MAF_ISKS_MGRB_df)
#GL_var_MAF_ISKS_MGRB_df$Tot_AC <- GL_var_MAF_ISKS_MGRB_df$Var_MAF_ISKS + GL_var_MAF_ISKS_MGRB_df$Var_MAF_MGRB
# library(reshape2) ##melting is not equivalent to SV analysis as it ignores the Allele Count across both MGRB and ISKS
# GL_melt <- melt(GL_var_MAF_ISKS_MGRB_df)
# colnames(GL_melt) <- c("VARIANT", "Allele_Count", "set", "Coh_Allele_Count")
# GL_melt$set <- gsub("Var_MAF_", "", GL_melt$set)
DF_add_AC <- function(AC_df, var_df)
{
  df_inp1 <- list()
  for(k in 1:length(unique(AC_df$Var_MAF_Tot))){
    AC_vec <- unique(AC_df$Var_MAF_Tot)
  AC_df_sub <- AC_df[AC_df$Var_MAF_Tot == AC_vec[k], ]
  df_sub <- var_df[var_df$VARIANT %in% AC_df_sub$VAR,]
  df_sub$Allele_Count <- AC_vec[k]
  df_inp1[[k]] <- df_sub
  }
  df_AC <- do.call("rbind.data.frame", df_inp1)
  return(df_AC)
}

GL_comb_inf_del_df_AC <- DF_add_AC(GL_var_MAF_ISKS_MGRB_df, GL_comb_inf_del_df)
##CDF
cdf_all_GL <- ggplot(GL_comb_inf_del_df_AC %>% filter(Allele_Count > 0), aes(x = Allele_Count)) +
  stat_ecdf(aes(color = set,linetype = set), 
            geom = "step", size = 1.5) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  labs(y = "cdf(Allele Count)") + ggtitle("GL CDF All events")

cdf_all_GL_trunc <- ggplot(GL_comb_inf_del_df_AC %>% filter(Allele_Count > 0 & Allele_Count < 5), aes(x = Allele_Count)) +
  stat_ecdf(aes(color = set,linetype = set), 
            geom = "step", size = 1.5) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  labs(y = "cdf(Allele Count)") + ggtitle("GL CDF Rare events")

multiplot(cdf_all_GL, cdf_all_GL_trunc, cols = 2)
##relationship to size
# GL_VAR_del_len_df <- unique(GL_comb_inf_del_df[,c(2,120)])
# GL_melt$del_len <- GL_VAR_del_len_df[match(GL_melt$VARIANT,GL_VAR_del_len_df$VARIANT),2]
t1 <- GL_comb_inf_del_df_AC %>% filter(Allele_Count > 0 & Allele_Count < 5)
GL_size_AC <- ggplot(t1 , aes(x=as.factor(Allele_Count), y=del_len, fill = set)) + 
  geom_boxplot() + 
  ggtitle("GL Size vs AC") + labs(x = "Allele Count", y = "Deletion size")

##Fisher test
#GL_dist <-  GL_melt %>% filter(Allele_Count > 0 & Allele_Count < 5)
GL_dist <- GL_comb_inf_del_df_AC %>% filter(Allele_Count > 0 & Allele_Count < 5)

fish_test_fun <- function(var_df, AC){
# nISKS <- length(unique(var_df[var_df$VARIANT %in% melt_df$VARIANT & var_df$set %in% "ISKS",]$SAMPLE))
# nMGRB <- length(unique(var_df[var_df$VARIANT %in% melt_df$VARIANT & var_df$set %in% "MGRB",]$SAMPLE))

nISKS <- length(unique(var_df[var_df$Allele_Count == AC & var_df$set %in% "ISKS",]$SAMPLE))
nMGRB <- length(unique(var_df[var_df$Allele_Count == AC & var_df$set %in% "MGRB",]$SAMPLE))

mat <- matrix(c(nISKS, 1644 - nISKS, nMGRB, 3205 - nMGRB), nrow = 2)
colnames(mat) <- c("ISKS", "MGRB")
rownames(mat) <- c("Yes", "No")
#print(mat)
#fisher.test(mat)
ft <- fisher.test(mat, conf.int = T, conf.level = 0.95)
ft_df <- cbind.data.frame("Allele_Count" =  AC, "P_Value" = ft$p.value,
                          "CI_lower" = ft$conf.int[1],
                          "CI_upper" = ft$conf.int[2],
                          "OR" = ft$estimate)

return(ft_df)
}
#GL_fish <- fish_test_fun(GL_dist, GL_comb_inf_del_df, "GL")
GL_fish <- fish_test_fun(GL_comb_inf_del_df_AC, 1)
##Map GL to CGC genes
# ISKS_genes <- unique(GL_comb_inf_del_df[GL_comb_inf_del_df$VARIANT %in% GL_dist$VARIANT & GL_comb_inf_del_df$set %in% "ISKS", ]$gene_symbol)
# MGRB_genes <- unique(GL_comb_inf_del_df[GL_comb_inf_del_df$VARIANT %in% GL_dist$VARIANT & GL_comb_inf_del_df$set %in% "MGRB", ]$gene_symbol)

ISKS_genes <- unique(GL_dist[GL_dist$set %in% "ISKS", ]$gene_symbol)
MGRB_genes <- unique(GL_dist[GL_dist$set %in% "MGRB", ]$gene_symbol)

mirabello_genes <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/mirabello_genes.txt",
                              sep = "", header = F, stringsAsFactors = F)
##remove trailing and leading white spaces
library(stringr)
mira_genes <- mirabello_genes %>% 
  mutate(V1 = str_trim(mirabello_genes$V1, side = "both"))
##mitotic check point genes
mito_chk_point = read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/mitotic_checkpoint/mitocheck_point.txt",
                            sep = "", header = F, skip = 1, stringsAsFactors = F)
##cep_haus C3C4C5
cep_haus_all <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/CEP_HAUS_C345.txt",
                           sep = "", header = F, stringsAsFactors = F)

top_SKAT_cgc <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/DT_skatbin_cgc.txt",
                           sep = "", header = F, stringsAsFactors = F)
cgc_genes <- read.delim("~/RVAS/cancer_gene_census_hg37.csv", sep = ",", header = T, stringsAsFactors = F)
#remove fusion genes from cgc list
cgc_genes$Gene.Symbol[(grep("fusion", cgc_genes$Role.in.Cancer))]
sheltrin_complex <- c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "ACD")
Telo_extension <- c("TIMELESS", "TIPIN", "FANCM", "BRCA1", "BLM") ##replication stress is a source of telomere recombination
Sheltrin_comp_extn = c("ACD", "POT1", "TERF1", "TERF2", "TERF2IP", "TINF2", "ATM", 
                       "BAG3", "BLM", "BRCA1", "CALD1", "CLK3", "DCLRE1B", "FANCD2", 
                       "FBXO4", "HSPA4", "KIAA1191", "MRE11A", "NBN", "PINX1", "PRKDC", 
                       "RAD50", "SLX4", "STUB1", "TNKS", "TNKS2", "U2AF2", "UCHL1", 
                       "WRN", "XRCC5", "XRCC6")

##Add mito_chk_point
##removed centrosome maturation genes
gene_sub <- unique(c(top_SKAT_cgc$V1, cgc_genes$Gene.Symbol, sheltrin_complex, Telo_extension, Sheltrin_comp_extn, 
                     mito_chk_point$V1, mira_genes$V1, cep_haus_all$V1))

table(ISKS_genes %in% gene_sub)
table(MGRB_genes %in% gene_sub)
gene_sub[gene_sub %in% ISKS_genes]
gene_sub[gene_sub %in% MGRB_genes]
ISKS_cgc_df <- GL_dist[GL_dist$set %in% "ISKS", ]
ISKS_cgc_df <-ISKS_cgc_df[ISKS_cgc_df$gene_symbol %in% gene_sub,]
write.table(ISKS_cgc_df, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/ISKS_cgc_inframe_del.tsv",
            sep = "\t", row.names = F, quote = F)
MGRB_cgc_df <- GL_dist[GL_dist$set %in% "MGRB", ]
MGRB_cgc_df <-MGRB_cgc_df[MGRB_cgc_df$gene_symbol %in% gene_sub,]
write.table(MGRB_cgc_df, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/MGRB_cgc_inframe_del.tsv",
            sep = "\t", row.names = F, quote = F)
##MGRB+ISKS Combined
GL_dist_cgc <-  GL_dist[GL_dist$gene_symbol %in% gene_sub,]
write.table(GL_dist_cgc, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/ISKS_MGRB_cgc_rare_inframe_del.tsv",
            sep = "\t", row.names = F, quote = F)
##CDF of CGC variants F(AC) vs AC
#GL_dist$gene <-  GL_comb_inf_del_df[match(GL_dist$VARIANT, GL_comb_inf_del_df$VARIANT), 9]

cdf_cgc_GL_trunc <- ggplot(GL_dist_cgc, aes(x = Allele_Count)) +
  stat_ecdf(aes(color = set,linetype = set), 
            geom = "step", size = 1.5) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  labs(y = "cdf(Allele Count)") + ggtitle("GL CGC: Rare Inframe Deletions")

##Fisher loop
fish_loop <- function(var_df, range){
  fish_out_df <- list()
  for(i in 1:range){
    print(i)
    var_df1 <- var_df[var_df$Allele_Count == i,]
    print(dim(var_df1))
    print(table(var_df1$set))
    fish_out_df[[i]] <- fish_test_fun(var_df1, i)
  }
  fish_out <- do.call("rbind.data.frame", fish_out_df)
  return(fish_out)
}

GL_cgc_fish <- fish_loop(GL_dist_cgc,4)
##Forest plot
p_forest <- ggplot(data=GL_cgc_fish, aes(y=Allele_Count, x=OR, xmin=CI_lower, xmax=CI_upper))+ 
  geom_point()+ geom_errorbarh(height=.1)+
  
  #sets the scales
  #note that I reverse the y axis to correctly order the effect #sizes based on my index variable
  scale_x_continuous(limits=c(-0.5,6), breaks = c(-0.5:6))+
  scale_y_continuous(name = "", breaks=1:4, labels = GL_cgc_fish$Allele_Count)+
  
  #adding a vertical line at the effect = 0 mark
  geom_vline(xintercept=1, color="black", linetype="dashed", alpha=.5)+
  ggtitle("Inframe Deletion Pathogenicity: Allele_Count vs OR")+
  theme_minimal()+
  theme(text=element_text(family="Arial",size=12, color="black"))+ labs(x = "OR", y = "Allele_Count")



##Somatic analysis
#Som_comb_inf_del_df

som_fun <- function(som_var_df){
  Som_comb_inf_del_VAR <- unique(som_var_df$VARIANT)
  Som_comb_inf_del_SAM <- unique(som_var_df$SAMPLE)
  
  match_var_Som <- function(var_inp){
    samp_id <- unique(som_var_df[som_var_df$VARIANT %in% var_inp, ]$SAMPLE)
    vec_sam <- ifelse(Som_comb_inf_del_SAM %in% samp_id, 1, 0)
    return(vec_sam)
  }
                    vec_list <- lapply(Som_comb_inf_del_VAR, function(x)match_var_Som(x))
                    
                    vec_list_mat <- do.call("rbind", vec_list)
                    colnames(vec_list_mat) <- Som_comb_inf_del_SAM
                    rownames(vec_list_mat) <- Som_comb_inf_del_VAR
                    print(dim(vec_list_mat))
                    ##MAF 
                    Var_MAF_ISKS <- rowSums(vec_list_mat[,grep("^[ABZ]", colnames(vec_list_mat), invert = T)])
                    # Var_MAF_ISKS_df <- as.data.frame(Var_MAF_ISKS)
                    # Var_MAF_ISKS_df$set <- "ISKS"
                    # colnames(Var_MAF_ISKS_df)[1] <- "Allele_Count"
                    Var_MAF_MGRB <- rowSums(vec_list_mat[,grep("^[ABZ]", colnames(vec_list_mat))])
                    # Var_MAF_MGRB_df <- as.data.frame(Var_MAF_MGRB)
                    # Var_MAF_MGRB_df$set <- "MGRB"
                    # colnames(Var_MAF_MGRB_df)[1] <- "Allele_Count"
                    # Var_MAF_ISKS_MGRB_df <- rbind.data.frame(Var_MAF_ISKS_df, Var_MAF_MGRB_df)
                    Var_MAF_ISKS_MGRB_df <- as.data.frame(cbind(Var_MAF_ISKS, Var_MAF_MGRB))
                    return(Var_MAF_ISKS_MGRB_df)
}

Som_var_MAF_ISKS_MGRB_df <- som_fun(Som_comb_inf_del_df)
Som_var_MAF_ISKS_MGRB_df$VAR <- rownames(Som_var_MAF_ISKS_MGRB_df)
library(reshape2)
som_melt <- melt(Som_var_MAF_ISKS_MGRB_df)
colnames(som_melt) <- c("VARIANT", "set", "Allele_Count")
som_melt$set <- gsub("Var_MAF_", "", som_melt$set)

##CDF
Som_all_GL <- ggplot(som_melt %>% filter(Allele_Count > 0), aes(x = Allele_Count)) +
  stat_ecdf(aes(color = set,linetype = set), 
            geom = "step", size = 1.5) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  labs(y = "cdf(Allele Count)") + ggtitle("Som CDF All events")

Som_all_GL_trunc <- ggplot(som_melt %>% filter(Allele_Count > 0 & Allele_Count < 5), aes(x = Allele_Count)) +
  stat_ecdf(aes(color = set,linetype = set), 
            geom = "step", size = 1.5) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  labs(y = "cdf(Allele Count)") + ggtitle("Som CDF Rare events")

multiplot(Som_all_GL, Som_all_GL_trunc, cols = 2)

##relationship to size
VAR_del_len_df <- unique(Som_comb_inf_del_df[,c(2,120)])
som_melt$del_len <- VAR_del_len_df[match(som_melt$VARIANT,VAR_del_len_df$VARIANT),2]
t1 <- som_melt %>% filter(Allele_Count > 0 & Allele_Count < 5)
Som_size_AC <- ggplot(t1 , aes(x=as.factor(Allele_Count), y=del_len, fill = set)) + 
  geom_boxplot() + 
  ggtitle("Somatic Size vs AC") + labs(x = "Allele Count", y = "Deletion size")


##fisher test
Som_dist <-  som_melt %>% filter(Allele_Count > 0 & Allele_Count < 5)
Som_fish <- fish_test_fun(Som_dist, Som_comb_inf_del_df, "Som")



##Summarise GL
library(ggpubr)
fisher_op <- rbind.data.frame(GL_fish, Som_fish)
tab <- ggtexttable(fisher_op, rows = NULL, theme = ttheme("classic",base_size = 12))
p_tab <- ggarrange(tab,
                labels = c("G"),ncol = 1)
gg_sum <- ggarrange(p_GL, p_Som, cdf_all_GL, cdf_all_GL_trunc, Som_all_GL, Som_all_GL_trunc,p_tab,
                labels = c("A", "B", "C", "D", "E", "F", "G"),ncol = 2, nrow = 4)

##GL plot
GLcgc_tab <- ggtexttable(GL_cgc_fish, rows = NULL, theme = ttheme("classic",base_size = 12))
gg_GL_cgc_tab <- ggarrange(cdf_cgc_GL_trunc, GLcgc_tab, labels = c("A", "B"),ncol = 1, nrow = 2)
multiplot(gg_GL_cgc, cols = 1)
gg_GL_cgc_trunc <- ggarrange(cdf_cgc_GL_trunc, p_forest, labels = c("A", "B"),ncol = 1, nrow = 2)
multiplot(gg_GL_cgc_trunc, cols = 1)

##check overlap with shortlisted+CGC genes

