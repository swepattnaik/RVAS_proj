##Make PID file
##with DT criteria: VAF >= 0.35 & Eigenphred >= 5.6
##PID file: Apr29-2020
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

summarise_var_per_samp <- function(variant_set){
  samvar_list <- list()
  sam_id <- unique(variant_set$SAMPLE)
  for(i in 1:length(sam_id)){
    v1 <- variant_set[variant_set$SAMPLE %in% sam_id[i],]
    # v1 <- v1[v1$VAF >= 0.35 & v1$comb_score >= 3,] ##comment for unfiltered variant set
    v1 <- v1[order(v1$comb_score, decreasing = T),]
    if(!is.null(v1)){
      ntot <- dim(v1)[1]
      nC5 <- sum(ifelse(v1$auto_call %in% "C5", 1, 0))
      geneC5 <- paste(v1[v1$auto_call %in% "C5",]$gene_symbol, collapse = ";")
      nC4 <- sum(ifelse(v1$auto_call %in% "C4", 1, 0))
      geneC4 <- paste(v1[v1$auto_call %in% "C4",]$gene_symbol, collapse = ";")
      nC3 <- sum(ifelse(v1$auto_call %in% "C3", 1, 0))
      geneC3 <- paste(v1[v1$auto_call %in% "C3",]$gene_symbol, collapse = ";")
      
      
    }else{    ntot = 0
    nC5 <- 0
    geneC5 <- NULL
    nC4 <- 0
    geneC4 <- NULL
    nC3 <- 0
    geneC4 <- NULL
    }
    gene_list <- list("EXT1_2" = c("EXT1", "EXT2"), 
                      "IDH1_2" = c("IDH1", "IDH2"),
                      "TP53" = c("TP53"),
                      "BRCA1_2_PALB2" = c("BRCA1", "BRCA2", "PALB2"),
                      "MMR_APC_MUTYH" = c("MMR", "APC", "MUTYH"),
                      "ATM_ATR_CHEK2_ARF_TP53" = c("ATM", "ATR", "CHEK2", "ARF", "TP53"),
                      "NF1_SDH" = c("NF1", "SDH"),
                      "Sheltrin_main" = c("POT1", "TINF2", "TERF1", "TERF2", "TERF2IP", "ACD"),
                      "Sheltrin_extn" = c("ACD", "POT1", "TERF1", "TERF2", "TERF2IP", "TINF2", "ATM", 
                                          "BAG3", "BLM", "BRCA1", "CALD1", "CLK3", "DCLRE1B", "FANCD2", 
                                          "FBXO4", "HSPA4", "KIAA1191", "MRE11A", "NBN", "PINX1", "PRKDC", 
                                          "RAD50", "SLX4", "STUB1", "TNKS", "TNKS2", "U2AF2", "UCHL1", 
                                          "WRN", "XRCC5", "XRCC6"))
    gene_list_match <- lapply(gene_list, function(x) ifelse(x %in% v1$gene_symbol, 1, 0))
    sum_list <- unlist(lapply(gene_list_match, function(x)sum(x)))
    sum_list <- as.data.frame(t(sum_list))
    samvar_list[[i]] <- cbind.data.frame("Total.burden" = ntot, "nC5" = nC5, "C5_genes" = geneC5,
                                         "nC4" = nC4, "C4_genes" = geneC4,
                                         "nC3" = nC3, "C3_genes" = geneC3, sum_list, "SAMPLE" = sam_id[i])
  }
  return(samvar_list)
}

##filtered set
fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/all_isksmgrb_combset2020_variants_filt_all_fields_rmdup.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)
fil_tab_isks_set <- fil_tab[fil_tab$is_case == 1,]
##use this for QC testing not in final (Discovery set topSKAT300)
#top300_SKAT <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Sep05_rect_ASP_graph_final_PPI_comb_GO_cmaf_new_score_str.tsv",
#                          sep = "\t", header = T, stringsAsFactors = F)

##use this in final
##SKATO
#top300_SKAT <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/ppi_res_fil_final.tsv",
#                          sep = "\t", header = T, stringsAsFactors = F)
##SKATBinary
top300_SKAT <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/ppi_res_fil_final_SKATbin.tsv",
                          sep = "\t", header = T, stringsAsFactors = F)
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
DT_genes <- read.table("~/RVAS/DT_genes_sep2018.txt", sep = "", stringsAsFactors = F) ##includes POT1 interactors
gene_sub <- unique(c(top300_SKAT$gene[1:300], cgc_genes$Gene.Symbol, sheltrin_complex, Telo_extension, DT_genes$V1, Sheltrin_comp_extn))
comb_set_filt <- fil_tab_isks_set[fil_tab_isks_set$gene_symbol %in% gene_sub,]
comb_set_filt$gnomad_AF <- gsub("\\[|\\]", "", comb_set_filt$gnomad_AF)
comb_set_filt <- comb_set_filt[comb_set_filt$VAF >= 0.35 & comb_set_filt$comb_score >= 5.6,]
#var_per_samp_isks <- summarise_var_per_samp(comb_set_filt)
#df_sam_filt_isks_comb <- do.call("rbind.data.frame", var_per_samp_isks)
#write.table(comb_set_filt, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/CGC_topSKATBin300_repstress_potint_VARIANTS_filt_combset2020_rmdup.tsv",
#            sep = "\t", row.names = F, quote = F)

##Add MAF based filter too: DT May_12_2020
Ex_samp_id <- unique(fil_tab$SAMPLE)
p_vec <- ifelse(!is.na(as.numeric(as.character(Ex_samp_id))) | grepl("^CR|^LK",Ex_samp_id), 1, 0)
#####make genotype matrix function
make_geno_mat <- function(ftemp_file){
  ftemp_tab_var_id <- unique(ftemp_file$VARIANT)
  samp_vec <- list()
  for(m in 1:length(ftemp_tab_var_id)){
    sam_gene_gt <- ftemp_file[ftemp_file$VARIANT %in% ftemp_tab_var_id[m],][,c(1:3,9,11,82,125:126)]
    sam_gene_gt <- unique(sam_gene_gt)
    if(dim(sam_gene_gt)[1] > 1 & length(unique(sam_gene_gt$vep_consequence)) > 1){
      # if(dim(sam_gene_gt)[1] > 1 & length(unique(sam_gene_gt$vep_consequence)) > 1 & length(unique(sam_gene_gt$comb_score)) > 1){
      vep_con <- unique(sam_gene_gt$vep_consequence)
      samp_vec_con <- list()
      for(k in 1:length(vep_con)){
        sam_gene_gt_con <- sam_gene_gt[sam_gene_gt$vep_consequence %in% vep_con[k],]
        sam_gene_gt_con <- unique(sam_gene_gt_con)
        maf_vec_cont <- sum(dim(sam_gene_gt_con[grepl("^[ABZ]",sam_gene_gt_con$SAMPLE),])[1])/(2*length(p_vec))
        maf_vec_case <- sum(dim(sam_gene_gt_con[!grepl("^[ABZ]",sam_gene_gt_con$SAMPLE),])[1])/(2*length(p_vec))
        ##genotype matrix  
        sam_gene_gt_con$add_mod <- as.numeric(sam_gene_gt_con$GT)
        sam10 <- ifelse(Ex_samp_id %in% sam_gene_gt_con$SAMPLE, 1, 0)
        sam10[which(sam10 != 0)] <- sam_gene_gt_con$add_mod ##additive model
        samp_vec_con[[k]] <- c(unique(sam_gene_gt_con$VARIANT),
                               unique(sam_gene_gt_con$gene_symbol), 
                               unique(sam_gene_gt_con$vep_consequence), 
                               unique(as.character(sam_gene_gt_con$auto_call)),
                               as.numeric(unique(sam_gene_gt_con$comb_score)), as.numeric(maf_vec_case),
                               as.numeric(maf_vec_cont), sam10)
      }
      
    }
    else{
      ##compute cohort specific MAF
      maf_vec_cont <- sum(ifelse(is.na(as.numeric(sam_gene_gt$SAMPLE)), 1, 0))/(2*length(p_vec))
      maf_vec_case <- sum(ifelse(!is.na(as.numeric(sam_gene_gt$SAMPLE)) | 
                                   grepl("^CR|^LK", as.character(sam_gene_gt$SAMPLE)), 1, 0))/(2*length(p_vec))
      #maf_vec <- (maf_vec_cont + maf_vec_case)/(2*(1572 + 1110))
      ##genotype matrix  
      sam_gene_gt$add_mod <- as.numeric(sam_gene_gt$GT)
      sam10 <- ifelse(Ex_samp_id %in% sam_gene_gt$SAMPLE, 1, 0)
      sam10[which(sam10 != 0)] <- sam_gene_gt$add_mod ##additive model
      ##account for discrepancy in vep_consequence
      #vep_con <- names(sort(table(unlist(strsplit(sam_gene_gt$vep_consequence, split = "&"))),decreasing=TRUE)[1])
      # sam_gene_gt$vep_consequence <- vep_con
      
      
      samp_vec[[m]] <- c(ftemp_tab_var_id[m],
                         unique(sam_gene_gt$gene_symbol), 
                         unique(sam_gene_gt$vep_consequence), 
                         unique(as.character(sam_gene_gt$auto_call)),
                         as.numeric(unique(sam_gene_gt$comb_score)), as.numeric(maf_vec_case),
                         as.numeric(maf_vec_cont), sam10)
    }
  }
  samp_vec_mat_uni <- do.call("rbind.data.frame", samp_vec)
  colnames(samp_vec_mat_uni) <- c("VARIANT", "gene_symbol", "vep_consequence", "auto_call", "comb_score", "coh_MAF_case",
                                  "coh_MAF_cont", Ex_samp_id)
  if(exists("samp_vec_con")){
    samp_vec_mat_con <- do.call("rbind.data.frame", samp_vec_con)
    colnames(samp_vec_mat_con) <- c("VARIANT", "gene_symbol", "vep_consequence", "auto_call", "comb_score", "coh_MAF_case",
                                    "coh_MAF_cont", Ex_samp_id)
    samp_vec_mat <- rbind.data.frame(samp_vec_mat_uni, samp_vec_mat_con)
  }else{
    samp_vec_mat <- samp_vec_mat_uni
  }
  print(dim(samp_vec_mat))
  ##Intracohort filter  : equivalent to ~ 5/1661*2 for case ; ~ 5/3209*2 for control
  ##changed to uniform intra_cohort_MAF filter
  ## change to MAF filter to 3/1661*2
  #     samp_vec_mat <- samp_vec_mat[as.numeric(as.character(samp_vec_mat$coh_MAF_case)) <= 0.0015 | 
  #                                    as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) <= 0.001,]
  samp_vec_mat <- samp_vec_mat[!(as.numeric(as.character(samp_vec_mat$coh_MAF_case)) >= 0.00035 | 
                                   as.numeric(as.character(samp_vec_mat$coh_MAF_cont)) >= 0.00035),]
  
  return(samp_vec_mat)
  
}

gene_var_maf <- function(gene_sym, var_file){
  
  #cut_off <- seq(0,10)
  #  score_scale <- c((0.001:50)^2/100,0.5,1000)
 # p_val_cutoff <- list()
 # for(i in 1:length(cut_off))  {
    ftemp_tab <- var_file[var_file$gene_symbol %in% gene_sym,] 
    ftemp_tab <- ftemp_tab[ftemp_tab$VAF >= 0.35 & ftemp_tab$comb_score >= 5.6,]
  #  if(dim(ftemp_tab)[1] > 1){
      samp_vec_mat <- make_geno_mat(ftemp_tab)
      ftemp_tab_mod <- samp_vec_mat[,c(1:7)]
      ftemp_tab_mod$comb_score <- as.numeric(as.character(ftemp_tab_mod$comb_score))
      ftemp_tab_mod$auto_call <- as.character(ftemp_tab_mod$auto_call)
 #   }
 #   else {
 #     next
 #   }
    return(ftemp_tab_mod)
  }

genes <- unique(comb_set_filt$gene_symbol)
write.table(genes, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/CGC_topSKATBin300_repstress_potint_PID_unique_genes_rmdup_maf.txt",
            sep = "\t", row.names = F, quote = F)

##Parallelised SKAT  
library(doParallel)
library(doMC)
registerDoMC(30)

res20 <- list()

genes <- unique(comb_set_filt$gene_symbol)
system.time(res20 <- foreach(i=1:length(genes), .errorhandling = 'remove') %dopar% 
{gene_var_maf(genes[i], comb_set_filt)})
gene_maf_pass_df <- do.call("rbind.data.frame", res20)
comb_set_filt1 <- comb_set_filt[comb_set_filt$VARIANT %in% gene_maf_pass_df$VARIANT,]
write.table(comb_set_filt1, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/CGC_topSKATBin300_repstress_potint_VARIANTS_filt_combset2020_rmdup_maf.tsv",
            sep = "\t", row.names = F, quote = F)
var_per_samp_isks <- summarise_var_per_samp(comb_set_filt1)
df_sam_filt_isks_comb <- do.call("rbind.data.frame", var_per_samp_isks)
# 
# for(k in 1:length(genes)){
#   print(k)
#   gene_maf_pass[[k]] <- gene_var_maf(genes[[k]], comb_set_filt)
# }
# gene_maf_pass_df <- do.call("rbind.data.frame", gene_maf_pass)
# comb_set_filt1 <- comb_set_filt[comb_set_filt$VARIANT %in% gene_maf_pass_df$VARIANT,]
# var_per_samp_isks <- summarise_var_per_samp(comb_set_filt1)
# df_sam_filt_isks_comb <- do.call("rbind.data.frame", var_per_samp_isks)
#write.table(comb_set_filt1, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/CGC_topSKATBin300_repstress_potint_VARIANTS_filt_combset2020_rmdup_maf.tsv",
#            sep = "\t", row.names = F, quote = F)

#######
##Add phenotypes
##Import telseq, CHIP and clinical information file to VM 
##Also add LION's clinical data
#detach("package:gdata", unload=TRUE)
library(readxl)
##Add telomere length
telseq <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/AllCohorts_telseq_telomeres.xlsx",
                     sheet = 1)
telseq <- as.data.frame(telseq)
##telseq lions
telseq_LK <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/Lions_telseq_67samples.tsv",
                        sep = "\t", header = F, stringsAsFactors = F)
telseq_LK1 <- cbind.data.frame("Cohort" = "lions", telseq_LK) 
colnames(telseq_LK1)[2:3] <- c("Sample", "LENGTH_ESTIMATE")

telseq_all <- rbind.data.frame(telseq, telseq_LK1)
##Chip
isks_chip <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/ISKS_chip1_chip2.chip2_all_categories_nb2_depth_and_vaf_and_position_tallies.tsv",
                        sep = "\t", header = T, stringsAsFactors = F)
risc_chip <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/riscSarcoma_chip1_chip2.chip2_all_categories_nb2_depth_and_vaf_and_position_tallies.tsv",
                        sep = "\t", header = T, stringsAsFactors = F)
risc_chip$sample_id <- paste0("CR", risc_chip$sample_id)
lions_chip <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/lionsSarcoma_chip1_chip2.chip2_all_categories_nb2_depth_and_vaf_and_position_tallies.tsv",
                         sep = "\t", header = T, stringsAsFactors = F)

isks_risc_lions_chip <- rbind.data.frame(isks_chip, risc_chip, lions_chip)
#telseq$Sample <- as.character(telseq$Sample)
##Mandy phenotype
#comb_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/Copy of ISKS_RisC_PID file 170120.xlsx", sheet = 2)
comb_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/PID_Master_file_290420.xlsx", sheet = 1, col_types = c("list"))
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


##Apr29-2020

##Add QC2 fail check
QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/joint_calls_2019jan_2019nov.final_qc_output.tsv", header = T, sep = "\t", stringsAsFactors = F)

QC2_dat_fail <- QC2_dat[QC2_dat$passes_qc2 %in% "FALSE",]$new_sampleid

`%nin%` = Negate(`%in%`)
comb_pheno$QC2 <- ifelse(as.character(comb_pheno$pmn) %in% QC2_dat_fail, "Fail", 
                         ifelse(as.character(comb_pheno$pmn) %nin% QC2_dat$new_sampleid, "Unknown", "Pass"))

########David's LION phenotype (not needed for Apr29_2020)
# lion_pheno <- read_excel("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set/var_indep/Lions_ZCC_germlines_M084.xlsx",
#                          sheet = 1)
# lion_pheno <- as.data.frame(lion_pheno)
# 
# lion_pheno1 <- lion_pheno
# 
# ##For Jan17
# #lion_pheno1[,c(12:55)] <- colnames(comb_pheno)
# #colnames(lion_pheno1)[12:55] <- colnames(comb_pheno)
# #lion_pheno1[,c(12:55)] <- NA
# ##for Apr29
# lion_pheno1[,c(12:77)] <- colnames(comb_pheno)
# colnames(lion_pheno1)[12:77] <- colnames(comb_pheno)
# lion_pheno1[,c(12:77)] <- NA
# 
# lion_pheno1$pid <- lion_pheno1$matched_normal_id
# lion_pheno1$pmn <- lion_pheno1$matched_normal_id
# lion_pheno1$gender1 <- ifelse(lion_pheno1$sex %in% "Female", "F", "M")
# lion_pheno1$CaseControl <- "Proband"
# lion_pheno1$`age at dateExtracted` <- lion_pheno1$age_at_sample
# lion_pheno1$SubjectCancers <- lion_pheno1$final_diagnosis
# lion_pheno1$SarcomaTopography <- lion_pheno1$cancer_category
# lion_pheno1$sarcomatype <- lion_pheno1$cancer_type
# 
# lion_pheno1$QC2 <- ifelse(as.character(lion_pheno1$pmn) %in% QC2_dat_fail, "Fail", "Pass")
# lion_pheno1 <- lion_pheno1[,-c(1:11)]
# lion_pheno1 <- lion_pheno1[lion_pheno1$pmn %in% QC2_dat$new_sampleid,]
##########
##Collate all phenotypes

#comb_ALL_phen <- rbind.data.frame(comb_pheno, lion_pheno1)
comb_ALL_phen <- comb_pheno
##these step has to be done after running summarise_var_per_samp function 
##Added Sheltrin_extn column

comb_ALL_phen[,c(67:82)] <- df_sam_filt_isks_comb[match(as.character(comb_ALL_phen$pmn), df_sam_filt_isks_comb$SAMPLE), c(1:16)]

comb_ALL_phen$telomere_length <- telseq_all[match(comb_ALL_phen$pmn, telseq_all$Sample), 3]

#comb_ALL_phen[,c(62:67)] <- isks_risc_lions_chip[match(comb_ALL_phen$pmn, isks_risc_lions_chip$sample_id), c(1:6)]

comb_ALL_phen[,c(84:89)] <- isks_risc_lions_chip[match(comb_ALL_phen$pmn, isks_risc_lions_chip$sample_id), c(1:6)]


##Variant file merged with phenotypes(not necessary, but generate later if needed)
#comb_set_filt[,c(133:198)] <- comb_ALL_phen[match(comb_ALL_phen$SAMPLE, comb_ALL_phen$pmn),c(1:66)]
#comb_set_filt_merged <- comb_set_filt[,-c(178:191)]
#write.table(var_set_merged, 
#            "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/Variants_CGC_topskat_repstress_genes_combet2020_clin_all_unfilt.tsv",
#            sep = "\t", row.names = F, quote = F)

####Prob_NFE
p_Data_comb <- read.table("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_combset_pca.scores_clustered.tsv", 
                          header = T, sep = "\t", stringsAsFactors = F)
p_Data_comb <- p_Data_comb[p_Data_comb$superPopulation %in% c("ISKS", "RISC", "LIONS"),]
p_Data_comb <- p_Data_comb[!duplicated(p_Data_comb$sample),]


##3137 and 3468 are present in replication set and are removed from PC clustering after relatedness check.
table(comb_ALL_phen$pmn %in% p_Data_comb$sample)
comb_ALL_phen$pred.NFE <- p_Data_comb[match(comb_ALL_phen$pmn, p_Data_comb$sample),51]
comb_ALL_phen$prob.NFE <- p_Data_comb[match(comb_ALL_phen$pmn, p_Data_comb$sample),52]

comb_ALL_phen_filt <- comb_ALL_phen[comb_ALL_phen$QC2 %in% "Pass" | comb_ALL_phen$QC2 %in% "Unknown",]
##mark duplicate samples
deg2_isks <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/dup_pairs_and_deg2_isks.tsv", sep = "\t",
                        header = T, stringsAsFactors = F)

deg2_isks <- deg2_isks[deg2_isks$degree == 0,]

comb_ALL_phen_filt$duplicate <- deg2_isks[match(comb_ALL_phen_filt$pmn, deg2_isks$name2), 7]

#comb_ALL_phen_filt$duplicate <- ifelse(comb_ALL_phen_filt$pmn %in% deg2_isks$name1, deg2_isks$name2,
#                                       ifelse(comb_ALL_phen_filt$pmn %in% deg2_isks$name2, deg2_isks$name1, "none"))

#write PID file

# write.table(comb_ALL_phen_filt[,-c(52:58)], 
#             "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/PID_combset2020_CGC_skatBin_repstress_potint_predNFE_rmdup.tsv",
#             sep = "\t", row.names = F, quote = F)
#write.table(comb_ALL_phen_filt[,-c(74:80)], 
#            "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/PID_combset2020_CGC_skatBin_repstress_potint_predNFE_rmdup_Apr292020.tsv",
#            sep = "\t", row.names = F, quote = F)
write.table(comb_ALL_phen_filt[,-c(74:80)], 
            "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/PID_combset2020_CGC_skatBin_repstress_potint_predNFE_rmdup_maf_Apr292020.tsv",
            sep = "\t", row.names = F, quote = F)
