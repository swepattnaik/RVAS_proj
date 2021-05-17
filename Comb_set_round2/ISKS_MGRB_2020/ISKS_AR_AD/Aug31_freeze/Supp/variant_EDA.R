##vep_input demo data is provided; use sarcoma internal control genes
var_file_isks <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/all_isks_combset2020_variants_AR_AD_all_fields_clingene_rnd3.tsv",
                            sep = "\t", header = T, stringsAsFactors = F)
var_file_isks <- var_file_isks[!is.na(var_file_isks$auto_call),]
var_file_isks$gnomad_AF_NFE <- as.numeric(var_file_isks$gnomad_AF_NFE)
var_file_isks$gnomad_AF_NFE <- ifelse(is.na(var_file_isks$gnomad_AF_NFE), 0, var_file_isks$gnomad_AF_NFE)
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(ggridges)

p2 <- ggplot(data=var_file_isks, aes(x=gnomad_AF_NFE, group=auto_call, fill=auto_call)) +
  geom_density(adjust=1.5, alpha=.4) 
p <- ggplot(data=var_file_isks, aes(x=gnomad_AF_NFE, group=auto_call, fill=auto_call)) +
  geom_density(adjust=1.5, position="fill")


rg_gph <- ggplot(data=var_file_isks, aes(x=gnomad_AF_NFE, y=as.factor(auto_call), fill=auto_call)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")
rg_gph
