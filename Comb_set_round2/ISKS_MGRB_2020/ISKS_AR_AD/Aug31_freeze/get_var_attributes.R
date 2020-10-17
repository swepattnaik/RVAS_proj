##get genomic and protein change and add to PID file specific to sheltrin table file
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)

library(readxl)
##read PID variant file
fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/comb_clin_isks_mgrb_2020_C3C4C5_NFE0002_AD_AR_all_fields_rnd3_Aug31.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)
fil_tab$is_case <- ifelse(grepl("^ISKS", fil_tab$set), 1, 0)
fil_tab_isks_set <- fil_tab[fil_tab$is_case == 1,]
fil_tab_isks_set$auto_call <- ifelse(fil_tab_isks_set$is_AR == 1 & fil_tab_isks_set$GT == 2, paste0(fil_tab_isks_set$auto_call, "_ar"), fil_tab_isks_set$auto_call)

##read Phenofile
##Mandy phenotype
#comb_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/Copy of ISKS_RisC_PID file 170120.xlsx", sheet = 2)
comb_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/PID_Master_file_290420AgeBloodTaken_Aug5.xlsx", sheet = 1, col_types = c("list"))
comb_pheno <- as.data.frame(comb_pheno)
comb_pheno1 <- sapply(comb_pheno, unlist)
colnames(comb_pheno1) <- colnames(comb_pheno)
comb_pheno <- comb_pheno1
comb_pheno <- as.data.frame(comb_pheno, stringsAsFactors = F)
comb_pheno <- unique(comb_pheno)
comb_pheno <- comb_pheno[!is.na(comb_pheno$pid),]
comb_pheno$`age at dateExtracted` <- as.numeric(comb_pheno$`age at dateExtracted`)
comb_pheno$AgeatSarcoma <- as.numeric(comb_pheno$AgeatSarcoma)
comb_pheno$SubjectAgeCancer <- as.numeric(comb_pheno$SubjectAgeCancer)

##Additional phenotypes added on Sept.8
add_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/extraPIDSdataset_100920.xlsx", sheet = 1, col_types = c("list"))
add_pheno <- as.data.frame(add_pheno)
add_pheno1 <- sapply(add_pheno, unlist)
colnames(add_pheno1) <- colnames(add_pheno)
add_pheno <- add_pheno1
add_pheno <- as.data.frame(add_pheno, stringsAsFactors = F)
add_pheno <- unique(add_pheno)
add_pheno <- add_pheno[!is.na(add_pheno$pid),]
add_pheno$`age at dateExtracted` <- as.numeric(add_pheno$`age at dateExtracted`)
add_pheno$AgeatSarcoma <- as.numeric(add_pheno$AgeatSarcoma)
add_pheno$SubjectAgeCancer <- as.numeric(add_pheno$SubjectAgeCancer)

colnames(add_pheno) <- colnames(comb_pheno)
##########
##Collate all phenotypes
comb_ALL_phen <- rbind.data.frame(comb_pheno, add_pheno)


##read sheltrin table file
# Sheltrin_tab_1 <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/Shelterin table.xlsx",
#                              sheet = 1)
# Sheltrin_tab_1 <- as.data.frame(Sheltrin_tab_1[,-c(11:13)])

##function to add attributes
get_var_attrib <- function(sheet_num){
#  xl_sheet_df <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/Shelterin_table_rect.xlsx",
#             sheet = sheet_num)
#  xl_sheet_df <- as.data.frame(xl_sheet_df[,-c(11:13)])
  xl_sheet_df <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/Shelterin_table_rect_sep182020.xlsx",
                            sheet = sheet_num)
  xl_sheet_df <- as.data.frame(xl_sheet_df)
  xl_sheet_df <- xl_sheet_df[!is.na(xl_sheet_df$ID),]
  xl_sheet_df$SAMPLE <- comb_ALL_phen[match(xl_sheet_df$ID, comb_ALL_phen$pid), 3]
  if(dim(xl_sheet_df)[2] == 10){
    xl_sheet_df$genomicclass <- comb_ALL_phen[match(xl_sheet_df$ID, comb_ALL_phen$pid), 25]
  }
  else {
    xl_sheet_df <- xl_sheet_df
  }
#  All_var1 <- fil_tab_isks_set[fil_tab_isks_set$SAMPLE %in% xl_sheet_df$SAMPLE & 
#                                 fil_tab_isks_set$gene_symbol %in% unique(xl_sheet_df$Gene),
#                               c(1:3,9,11,77,79:82)]

xl_sheet_df$mult_gene <- ifelse(grepl("/",xl_sheet_df$Gene), 1, 0)
if(sum(xl_sheet_df$mult_gene) > 0){
row_mult <- xl_sheet_df[which(xl_sheet_df$mult_gene == 1),]
row_mult_genes <- unlist(strsplit(row_mult$Gene, "/"))
#All_var2 <- fil_tab_isks_set[fil_tab_isks_set$SAMPLE %in% row_mult$SAMPLE & 
#                                   fil_tab_isks_set$gene_symbol %in% row_mult_genes,
#                                 c(1:3,9,11,77,79:82)]
#All_var12 <- rbind.data.frame(All_var1, All_var2)

add_rows <- list()
for(i in 1:length(row_mult_genes)){
  row_new <- row_mult
  row_new$Gene <- gsub(row_mult$Gene, row_mult_genes[i], row_mult$Gene)
  row_new$mult_gene <- 0
  add_rows[[i]] <- row_new
                                  }
xl_sheet_df <- xl_sheet_df[-which(xl_sheet_df$mult_gene == 1),]
xl_sheet_df <- rbind.data.frame(xl_sheet_df, do.call("rbind.data.frame", add_rows))
tab_var <- list()
for(k in 1:dim(xl_sheet_df)[1]){
  df_row <- xl_sheet_df[k,]
  tab_var[[k]] <-  fil_tab_isks_set[fil_tab_isks_set$SAMPLE %in% df_row$SAMPLE &
                                      fil_tab_isks_set$gene_symbol %in% df_row$Gene,c(1:3,9,11,77,79:82)]
  ##Add column 74 for exon number
}
tab_var_df <- do.call("rbind.data.frame", tab_var)
tab_var_df$vep_HGVSc.1 <- gsub("^ENSP.*:","", tab_var_df$vep_HGVSc.1)
xl_sheet_df <- xl_sheet_df[,-12]
xl_sheet_df_new <- cbind.data.frame(xl_sheet_df, tab_var_df)

return(xl_sheet_df_new)
                                      }
else{
xl_sheet_df <- xl_sheet_df[,-12]
#xl_sheet_df <- xl_sheet_df[,-11]
tab_var <- list()
for(k in 1:dim(xl_sheet_df)[1]){
  df_row <- xl_sheet_df[k,]
  tab_var[[k]] <-  fil_tab_isks_set[fil_tab_isks_set$SAMPLE %in% df_row$SAMPLE &
                                      fil_tab_isks_set$gene_symbol %in% df_row$Gene,c(1:3,9,11,77,79:82)]
  ##Add column 74 for exon number
}
tab_var_df <- do.call("rbind.data.frame", tab_var)
tab_var_df$vep_HGVSc.1 <- gsub("^ENSP.*:","", tab_var_df$vep_HGVSc.1)
xl_sheet_df_new <- cbind.data.frame(xl_sheet_df, tab_var_df)
return(xl_sheet_df_new)
}
}


sh1 <- get_var_attrib(sheet_num = 1) ##fix IDs for TERF1 variant in 1560; it should be 2370
sh2 <- get_var_attrib(sheet_num = 2)
sh3 <- get_var_attrib(sheet_num = 3)

write.table(sh1, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/Shelterin_tab_sh1.tsv",
            sep = "\t", row.names = F, quote = F)
write.table(sh2, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/Shelterin_tab_sh2.tsv",
            sep = "\t", row.names = F, quote = F)
write.table(sh3, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/Shelterin_tab_sh3.tsv",
            sep = "\t", row.names = F, quote = F)
