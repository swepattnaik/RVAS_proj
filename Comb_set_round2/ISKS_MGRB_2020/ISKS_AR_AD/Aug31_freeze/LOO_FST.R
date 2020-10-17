##FST package for variant effect scores and use in LOO and gene set test
#https://cran.r-project.org/web/packages/FSTpackage/FSTpackage.pdf
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
`%nin%` = Negate(`%in%`)
library(FSTpackage)
###get input files and process
##Experiments

fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_all_fields_rnd3_freeze.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)

##Additive model
Ex_samp_id <- unique(fil_tab$SAMPLE)


######get phenotype data to control for age and sex and PC's; 

QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv", 
                      header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %in% "TRUE",]
QC2_dat_pass$isFemale <- ifelse(QC2_dat_pass$f_stat < 0.2, 1, 
                                ifelse(QC2_dat_pass$f_stat > 0.8, 0, 2))
##note p_Data_PC_comb does not change
#p_Data <- read.table("~/RVAS/ISKS_MGRB_gender_age_3665.tsv", header = T, sep = "\t")
#p_Data <- read.table("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_pca_combset.scores.tsv", header = T, sep = "\t", stringsAsFactors = F)
##for future QC with prob.NFE correction use the following file and adjust the column for gender from 39 to 53
p_Data <- read.table("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_combset_pca.scores_clustered.tsv", 
                     header = T, sep = "\t", stringsAsFactors = F)
p_Data <- p_Data[p_Data$superPopulation %in% c("ISKS", "RISC", "LIONS", "MGRB"),]
p_Data_MGRB <-  p_Data[p_Data$superPopulation %in% c("MGRB"),]
#remove duplicates
dup_samp <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/fin_samp_dup_drop.txt", sep = "", header = T,
                       stringsAsFactors = F)
dup_samp <- dup_samp$x[grepl("^[ABZ]", dup_samp$x)]
p_Data_MGRB <- p_Data_MGRB[p_Data_MGRB$sample %nin% dup_samp,]
p_Data_MGRB$rect_sam <- p_Data_MGRB$sample
p_Data_MGRB <- p_Data_MGRB[p_Data_MGRB$rect_sam %in% QC2_dat_pass$new_sampleid,]
##rename p_Data
p_Data_ISKS <-  p_Data[p_Data$superPopulation %in% c("ISKS", "RISC", "LIONS"),]
rect_sam_dat <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/ISKS_RISC_LIONS_final_freeze.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
p_Data_ISKS$rect_sam <- rect_sam_dat[match(p_Data_ISKS$sample,rect_sam_dat$JCInputID),9]
p_Data_ISKS <- p_Data_ISKS[!is.na(p_Data_ISKS$rect_sam),]
##two samples missed due to mislabelling and QC2
#rect_sam_dat$JCInputRecID[rect_sam_dat$JCInputRecID %nin% p_Data_ISKS$sample]
p_Data_ISKS <- p_Data_ISKS[p_Data_ISKS$rect_sam %in% QC2_dat_pass$new_sampleid,] ##3105 is lost(QC2 fail)
p_Data_noCH <- rbind.data.frame(p_Data_ISKS, p_Data_MGRB)
p_Data_noCH <- p_Data_noCH[as.character(p_Data_noCH$rect_sam) %in% Ex_samp_id,]
##drop MGRB sample BAAUD; not in latest call
#samp_ID_match <- samp_ID_match[grep("BAAUD", samp_ID_match, invert = T)]
#p_Data_noCH <- p_Data_noCH[match(samp_ID_match, p_Data_noCH$sample),]
Ex_samp_id <- Ex_samp_id[match(p_Data_noCH$rect_sam, Ex_samp_id)]
##filter out QC fail cases
fil_tab <- fil_tab[fil_tab$SAMPLE %in% Ex_samp_id,]
##save table
#write.table(fil_tab, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/all_isksmgrb_skatinp_combset2020_clin_C3C4C5_NFE0002_AD_rm_dup_freeze.tsv", 
#            sep = "\t", row.names = F, quote = F)
##remove from gene matrix
#col_rm <- which(colnames(D_tab) ==  "BAAUD")
#D_tab <- D_tab[,-c(col_rm)]
##Add gender information

p_Data_noCH$gender <- QC2_dat_pass[match(p_Data_noCH$rect_sam, QC2_dat_pass$new_sampleid), 21]

#binary phenotype vector 
p_vec <- ifelse(!is.na(as.numeric(as.character(p_Data_noCH$rect_sam))) | grepl("^CR|^LK",as.character(p_Data_noCH$rect_sam)), 1, 0)

data(FST.example)
Y<-FST.example$Y;X<-FST.example$X;G<-FST.example$G;Z<-FST.example$Z
result.prelim<-FST.prelim(Y,X=X,out_type='D')
result<-FST.test(result.prelim,G,Z,B=5000)
Z1 <- cbind(0, Z[,2])
result1<-FST.test(result.prelim,G,Z,B=5000)

##do the calculation in SKAT script
gene_mat_comb <- as.numeric(as.character(samp_vec_mat[,5]))
gene_mat <- as.matrix(samp_vec_mat[,-c(1:7)])
G1 <- t(gene_mat)
G2 <- apply(G1, 2, function(x)x/(2*length(p_vec)))
X <- as.matrix(p_Data_noCH[,c(54,19:22)])
Z2 <- cbind(1,gene_mat_comb)
Y1 <- p_vec
res.null <- FST.prelim(Y=Y1,X=X,out_type='D')
result_TP53<-FST.test(res.null,G2,Z2,B=5000)
Z3 <- cbind(1, 10*Z2[,2])
t1 <- FST.test(res.null,G2,Z3,B=5000)

