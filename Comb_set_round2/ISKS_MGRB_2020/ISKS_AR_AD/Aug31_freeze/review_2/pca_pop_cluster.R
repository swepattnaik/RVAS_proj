#!/usr/bin/env Rscript

#.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

#setwd("~/RVAS/comb_set_2020/pop_PCA/")
setwd("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/Aug31_freeze/review_2")
`%nin%` = Negate(`%in%`)
library(mclust)

scores = read.table("~/RVAS/comb_set_2020/pop_PCA/MGRB_ISKS_1000G_pca_combset.scores.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

scores$population <- ifelse(is.na(scores$population), "ISKSMGRB", scores$population)
scores$superPopulation <- ifelse(is.na(scores$superPopulation), "ISKSMGRB", scores$superPopulation)

scores$superPopulation <- ifelse(!is.na(as.numeric(scores$sample)), "ISKS", 
                                 ifelse(grepl("^CR", scores$sample), "RISC", 
                                        ifelse(grepl("^LK", scores$sample), "LIONS", scores$superPopulation)))
scores$superPopulation <- gsub("ISKSMGRB", "MGRB", scores$superPopulation)

# scores$population = ifelse(scores$population == "ISKSMGRB" & grepl("^[ABZ]", scores$sample), "MGRB", 
#                            ifelse(scores$population == "ISKSMGRB" & !(grepl("^[ABZ]", scores$sample)), "ISKS", 
#                                   scores$population))
# scores$superPopulation = ifelse(scores$superPopulation == "ISKSMGRB" & grepl("^[ABZ]", scores$sample), "MGRB", 
#                                 ifelse(scores$superPopulation == "ISKSMGRB" & !(grepl("^[ABZ]", scores$sample)), "ISKS", 
#                                        scores$superPopulation))

pred_vars = c("PC1", "PC2", "PC3", "PC4")

`%nin%` = Negate(`%in%`)
#model_superpop = MclustDA(scores[scores$superPopulation != "ISKSMGRB",pred_vars], scores$superPopulation[scores$superPopulation != "ISKSMGRB"], modelType = "EDDA")
model_superpop = MclustDA(scores[scores$superPopulation %nin% c("ISKS","MGRB", "RISC", "LIONS"),pred_vars], 
                          scores$superPopulation[scores$superPopulation %nin% c("ISKS","MGRB", "RISC", "LIONS")], modelType = "EDDA")
summary(model_superpop)
plot(model_superpop, "scatterplot")
model_eur = MclustDA(scores[scores$superPopulation == "EUR",pred_vars], scores$population[scores$superPopulation == "EUR"], modelType = "EDDA")
summary(model_eur)
png("Eur_model_EDDA.png", height = 800, width = 800)
plot(model_eur, "scatterplot")
dev.off()

temp = predict(model_superpop, scores[,pred_vars])
scores$pred.superPop = temp$classification
temp = temp$z
colnames(temp) = paste("pred.superPop.", colnames(temp), sep = "")
scores = cbind(scores, temp)

temp = predict(model_eur, scores[,pred_vars])
scores$pred.eurPop = temp$classification
scores$pred.eurPop[scores$pred.superPop != "EUR"] = NA
temp = temp$z
temp[scores$pred.superPop != "EUR",] = NA
colnames(temp) = paste("pred.eurPop.", colnames(temp), sep = "")
scores = cbind(scores, temp)

scores$pred.NFE = FALSE
scores$pred.NFE[scores$pred.superPop == "EUR" & scores$pred.eurPop != "FIN"] = TRUE
scores$prob.NFE = 1 - scores$pred.eurPop.FIN

t1 = table(scores$pred.superPop, scores$superPopulation)
table(scores$pred.NFE, scores$superPopulation)
t2 = t1[,c(5:6,8,7)]
mat_class = matrix(t2, nrow = 5, ncol = 4)
colnames(mat_class) = colnames(t2)
rownames(mat_class) = rownames(t2)

write.table(scores, "MGRB_ISKS_1000G_combset_pca.scores_clustered.tsv", col.names = TRUE, row.names = FALSE,
            quote = F, sep = "\t")
saveRDS(mat_class, "MGRB_ISKS_1000G_pred_class.rds", compress = T)


##generate table score table after correction for duplication and relatedness

##QC pass filter
QC2_dat <- read.table("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/Rep_set_variants/joint_calls_2019jan_2019nov.final_qc_output.tsv", 
                      header = T, sep = "\t", stringsAsFactors = F)
QC2_dat_pass <- QC2_dat[QC2_dat$passes_qc2 %in% "TRUE",]
QC2_dat_pass$isFemale <- ifelse(QC2_dat_pass$f_stat < 0.2, 1, 
                                ifelse(QC2_dat_pass$f_stat > 0.8, 0, 2))

scores = read.table("MGRB_ISKS_1000G_pca_combset.scores.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

scores$population <- ifelse(is.na(scores$population), "ISKSMGRB", scores$population)
scores$superPopulation <- ifelse(is.na(scores$superPopulation), "ISKSMGRB", scores$superPopulation)

scores$superPopulation <- ifelse(!is.na(as.numeric(scores$sample)), "ISKS", 
                                 ifelse(grepl("^CR", scores$sample), "RISC", 
                                        ifelse(grepl("^LK", scores$sample), "LIONS", scores$superPopulation)))
scores$superPopulation <- gsub("ISKSMGRB", "MGRB", scores$superPopulation)

scores_1000G = scores[scores$population %nin% "ISKSMGRB", ]
scores_ISKS_MGRB = scores[scores$population %in% "ISKSMGRB", ]
scores_MGRB = scores_ISKS_MGRB[scores_ISKS_MGRB$superPopulation %in% "MGRB",]
scores_ISKS = scores_ISKS_MGRB[scores_ISKS_MGRB$superPopulation %nin% "MGRB",]
  
#remove duplicates
dup_samp <- read.delim("~/RVAS/comb_set_2020/PC_relate/ldpruned/fin_samp_dup_drop.txt", sep = "", header = T,
                         stringsAsFactors = F)
dup_samp <- dup_samp$x[grepl("^[ABZ]", dup_samp$x)]
scores_MGRB <- scores_MGRB[scores_MGRB$sample %nin% dup_samp,]
scores_MGRB$rect_sam <- scores_MGRB$sample
scores_MGRB <- scores_MGRB[scores_MGRB$rect_sam %in% QC2_dat_pass$new_sampleid,]

rect_sam_dat <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/PID/ISKS_RISC_LIONS_final_freeze.tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
scores_ISKS$rect_sam <- rect_sam_dat[match(scores_ISKS$sample,rect_sam_dat$JCInputID),9]
scores_ISKS <- scores_ISKS[!is.na(scores_ISKS$rect_sam),]
# write.table(p_Data_ISKS, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/p_Data_ISKS.tsv",
#             sep = "\t", quote = F, row.names = F)
##two samples missed due to mislabelling and QC2
#rect_sam_dat$JCInputRecID[rect_sam_dat$JCInputRecID %nin% p_Data_ISKS$sample]
scores_ISKS <- scores_ISKS[scores_ISKS$rect_sam %in% QC2_dat_pass$new_sampleid,] ##3105 is lost(QC2 fail)
scores_ISKS_MGRB <- rbind.data.frame(scores_ISKS, scores_MGRB)
scores_1000G$rect_sam = scores_1000G$sample
scores = rbind.data.frame(scores_1000G, scores_ISKS_MGRB)

##model building
pred_vars = c("PC1", "PC2", "PC3", "PC4")

`%nin%` = Negate(`%in%`)
#model_superpop = MclustDA(scores[scores$superPopulation != "ISKSMGRB",pred_vars], scores$superPopulation[scores$superPopulation != "ISKSMGRB"], modelType = "EDDA")
model_superpop = MclustDA(scores[scores$superPopulation %nin% c("ISKS","MGRB", "RISC", "LIONS"),pred_vars], 
                          scores$superPopulation[scores$superPopulation %nin% c("ISKS","MGRB", "RISC", "LIONS")], modelType = "EDDA")
summary(model_superpop)
plot(model_superpop, "scatterplot")
model_eur = MclustDA(scores[scores$superPopulation == "EUR",pred_vars], scores$population[scores$superPopulation == "EUR"], modelType = "EDDA")
summary(model_eur)
#png("Eur_model_EDDA.png", height = 800, width = 800)
#plot(model_eur, "scatterplot")
#dev.off()

temp = predict(model_superpop, scores[,pred_vars])
scores$pred.superPop = temp$classification
temp = temp$z
colnames(temp) = paste("pred.superPop.", colnames(temp), sep = "")
scores = cbind(scores, temp)

temp = predict(model_eur, scores[,pred_vars])
scores$pred.eurPop = temp$classification
scores$pred.eurPop[scores$pred.superPop != "EUR"] = NA
temp = temp$z
temp[scores$pred.superPop != "EUR",] = NA
colnames(temp) = paste("pred.eurPop.", colnames(temp), sep = "")
scores = cbind(scores, temp)

scores$pred.NFE = FALSE
scores$pred.NFE[scores$pred.superPop == "EUR" & scores$pred.eurPop != "FIN"] = TRUE
scores$prob.NFE = 1 - scores$pred.eurPop.FIN

t1 = table(scores$pred.superPop, scores$superPopulation)
table(scores$pred.NFE, scores$superPopulation)
t2 = t1[,c(5:6,8,7)]
mat_class = matrix(t2, nrow = 5, ncol = 4)
colnames(mat_class) = colnames(t2)
rownames(mat_class) = rownames(t2)

# write.table(scores, "MGRB_ISKS_1000G_combset_pca.scores_clustered_rect_sam.tsv", col.names = TRUE, row.names = FALSE,
#             quote = F, sep = "\t")
# saveRDS(mat_class, "MGRB_ISKS_1000G_pred_class_rect_sam.rds", compress = T)


