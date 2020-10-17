##C3 Variant scoring model
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

library(gdata)
library("MASS")
library(ggplot2)
library(gridExtra)
library(geoR) #for 3-D plots
library(caret)
library(pscl)
library(ISLR)
library(knitr)

var_file_isks <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/all_isks_combset2020_variants_AR_AD_all_fields_clingene_rnd3.tsv",
                       sep = "\t", header = T, stringsAsFactors = F)
var_file_mgrb <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/MGRB/all_mgrb_combset2020_variants_AR_AD_all_fields_clingene_rnd3.tsv",
                            sep = "\t", header = T, stringsAsFactors = F)
var_file_mgrb$set <- gsub("ISKS", "MGRB", var_file_mgrb$set)
var_file <- rbind.data.frame(var_file_isks, var_file_mgrb)
var_file_C3 <- var_file[var_file$auto_call %in% "C3",]
var_file_sel_C3 <- var_file_C3[,c(2,9,53,57,58,60:61,130)]
var_file_sel_C3$is_case <- as.factor(ifelse(var_file_sel_C3$set %in% "ISKS_AR_AD", 1, 0))
ggplot(var_file_sel_C3, aes(x = is_case, y = EigenPhred, fill = is_case)) + geom_boxplot()
#set.seed(12093)
#s1 <- sample(x = unique(var_file_sel_C3$gene_symbol), size = 1000)

#train_set <- var_file_sel_C3[var_file_sel_C3$gene_symbol %in% unique(c("TP53", "NF1", "BUB1B", "BRCA1", "BRCA2",
#                                                                       "MSH2", "EXT1", "EXT2", s1)),]

train_set <- var_file_sel_C3[var_file_sel_C3$gene_symbol %in% c("TP53", "NF1", "BRCA2", "ERCC2", "EXT1", "EXT2",
                                                                "SDHA", "SDHB", "SDHD"),]

#cgc_genes <- read.delim("~/RVAS/cancer_gene_census_hg37.csv", sep = ",", header = T, stringsAsFactors = F)
#train_set <- var_file_sel_C3[var_file_sel_C3$gene_symbol %in% cgc_genes$Gene.Symbol,]
dat_num <- train_set[,c(9,3:5,7)]
dat_num$prod_eig_cond <- dat_num$EigenPhred * dat_num$CONDEL_Score
dat_num$mean_eig_cadd <- (dat_num$EigenPhred + dat_num$CADD_PHRED)/2
dat_num$AM <- apply(dat_num[,2:5], 1, mean)
##based on sensitivity specificity cut-off filter the dat_num$EigenPhred >= 6
#dat_num <- dat_num[dat_num$EigenPhred >= 6,]

##logistic regression model

logistic.regression.or.ci <- function(regress.out, level=0.95)
{
  ################################################################
  #                                                              #
  #  This function takes the output from a glm                   #
  #  (logistic model) command in R and provides not              #
  #  only the usual output from the summary command, but         #
  #  adds confidence intervals for all coefficients and OR's.    #
  #                                                              #
  #  This version accommodates multiple regression parameters    #
  #                                                              #
  ################################################################
  usual.output <- summary(regress.out)
  z.quantile <- qnorm(1-(1-level)/2)
  number.vars <- length(regress.out$coefficients)
  OR <- exp(regress.out$coefficients[-1])
  temp.store.result <- matrix(rep(NA, number.vars*2), nrow=number.vars)
  for(i in 1:number.vars)
  {
    temp.store.result[i,] <- summary(regress.out)$coefficients[i] +
      c(-1, 1) * z.quantile * summary(regress.out)$coefficients[i+number.vars]
  }
  intercept.ci <- temp.store.result[1,]
  slopes.ci <- temp.store.result[-1,]
  OR.ci <- exp(slopes.ci)
  output <- list(regression.table = usual.output, intercept.ci = intercept.ci,
                 slopes.ci = slopes.ci, OR=OR, OR.ci = OR.ci)
  return(output)
}
####Function:ROC plots
library(ROCR)
plot_roc <- function(prob, resp_var){
  pr <- prediction(prob, resp_var)
  prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  #plot(prf)
  auc <- performance(pr, measure = "auc")
  auc <- auc@y.values[[1]]
  #print("AUC")
  #  return(auc)
  return(list(auc, prf))
}
####Function: Univariate model
anv <- numeric()
mod_dev <- numeric()
null_dev <- numeric()
univar_op <- list()
pse_r_sq <- numeric()  #MacFadden's pseudo R-squared test
pred <- numeric()
pred_obj <- list()
auc_roc <- list()
col_var <- colnames(dat_num)
for(i in 1:dim(dat_num)[2]){
  output <- glm(is_case ~ dat_num[,i], data=dat_num, family=binomial(link="logit"))
  pred <- predict(output,newdata = data.frame(dat_num[,i]), type = "response")
  pred_obj[[i]] <- pred
  # summary(output)
  a <- anova(output, test = "Chisq")
  anv[i] <- a$`Pr(>Chi)`[2]
  mod_dev[i]<- a$Deviance[2]
  null_dev[i] <- a$`Resid. Dev`[1]
  pse_r_sq[i] <- pR2(output)[4]
  univar_op[[i]] <- logistic.regression.or.ci(output)
  auc_roc[[i]] <- plot_roc(pred, dat_num$is_case)
}
names(auc_roc) <- colnames(dat_num)
names(univar_op) <- colnames(dat_num)
Est_cor <- numeric()
pval_cov <- numeric()
odds_rat <- numeric()
OR_CI_lower <- numeric()
OR_CI_upper <- numeric()
AUC_chk <- numeric()
for(j in 1:length(univar_op)){
  Est_cor[j] <- univar_op[[j]]$regression.table$coefficients[2,1]
  pval_cov[j] <- univar_op[[j]]$regression.table$coefficients[2,4]
  odds_rat[j] <- univar_op[[j]]$OR
  OR_CI_lower[j] <- univar_op[[j]]$OR.ci[1]
  OR_CI_upper[j] <- univar_op[[j]]$OR.ci[2]
  AUC_chk[j] <- auc_roc[[j]][[1]]
}

df_univar <- as.matrix(cbind("Estimate" = Est_cor, "P value" = pval_cov, "OR" = odds_rat, "CI_lower" = OR_CI_lower, "CI_upper" = OR_CI_upper, "Chi-sq_sig" = anv,
                             "Deviance" = mod_dev, "Null_deviance" = null_dev, "Pseudo-R.sq" = pse_r_sq, "AUC" = AUC_chk))
rownames(df_univar) <- colnames(dat_num)

df_univar_unfilt <- df_univar
df_univar_unfilt <- df_univar_unfilt[-1,]
#df_univar <- df_univar[-c(2,11:13),]
#rownames(df_univar) <- colnames(dat_num)[-c(2,11:13)]

##Plot AUC

library(RColorBrewer)
set.seed(2236789)
#colr <- sample(brewer.pal(n = 7, name = "Dark2" ), 7, replace = F)
colr <- sample(brewer.pal(n = 7, name = "Paired" ), 5, replace = F)
plot(0,0,xlim = c(0,1),ylim = c(0,1),type = "n", xlab =  "False Positive Rate", ylab = "True Positive Rate")


for(i in 2:5){
  
  plot(auc_roc[[i]][[2]], col = colr[i], avg = "vertical", spread.estimate = "stderror",
       show.spread.at = seq(0.1, 0.8, 0.1), plotCI.col = colr[i], plotCI.lwd = 2, lwd = 2, add = TRUE)
  # plot(auc_roc[[1]][[i]], col = "blue", lty = 2, lwd = 0.25, add = TRUE)
  legend("bottomright", legend = names(auc_roc)[2:5], lty = 1, lwd =2, col = colr[2:5], cex = 0.8)
}

plot(0,0,xlim = c(0,1),ylim = c(0,1),type = "n", xlab =  "False Positive Rate", ylab = "True Positive Rate")
plot(auc_roc[[2]][[2]], col = colr[2], avg = "vertical", spread.estimate = "stderror",
     show.spread.at = seq(0.1, 0.8, 0.1), plotCI.col = colr[2], plotCI.lwd = 2, lwd = 2, add = TRUE)
plot(auc_roc[[6]][[2]], col = colr[5], avg = "vertical", spread.estimate = "stderror",
     show.spread.at = seq(0.1, 0.8, 0.1), plotCI.col = colr[6], plotCI.lwd = 2, lwd = 2, add = TRUE)
legend("bottomright", legend = names(auc_roc)[c(2,6)], lty = 1, lwd =2, col = colr[c(2,7)], cex = 0.8)


##transforming cut_off prob to predictor value
#log[p/(1-p)] = b0 + b1*x
#x = (log[p/(1-p)] - b0)/b1
output_Eigen <- glm(is_case ~ EigenPhred, data=dat_num, family=binomial(link="logit"))
P_Eig <- predict(output_Eigen,newdata = data.frame(EigenPhred = dat_num$EigenPhred), type = "response")
p = mean(P_Eig)

p = 0.420
b0 = coef(output_Eigen)[1]
b1 = coef(output_Eigen)[2]

x = (log(p/(1-p)) - b0)/b1
p = 0.473
x = (log(p/(1-p)) - b0)/b1
#check Eigen cutoff of 5
x = 5

sig_val = exp(b0 + b1*x) 

p = sig_val/(1 + sig_val)

#confusion matrix
cutoff <- seq(0.4,0.55, by = 0.001)
# predicted_values <- ifelse(P_Eig >= mean(P_Eig), 1 , 0)
choose_cutoff <- function(pred_val, data){
  accuracy <- as.numeric()
  spec <- as.numeric()
  sens <- as.numeric()
  for(i in 1:length(cutoff)){
    predicted_values <- ifelse(pred_val >= cutoff[i], 1 , 0)
    conf_matrix<-table(predicted_values,data[,1])
    conf_matrix <- t(conf_matrix)
    accuracy[i]<-(conf_matrix[1,1]+conf_matrix[2,2])/(sum(conf_matrix))
    #spec = TN/(TN+FP) = TNR
    #sens = TP/(TP+FN) = TPR
    spec[i] = conf_matrix[1,1]/(conf_matrix[1,1] + conf_matrix[1,2]) #controls as controls
    sens[i] = conf_matrix[2,2]/(conf_matrix[2,2] + conf_matrix[2,1]) #cases as cases
    print(i)
  }
  # plot(accuracy ~ cutoff )
  mat_comb <- cbind.data.frame("accuracy" = accuracy, "cutoff" = cutoff, "spec" = spec, "sens" = sens)
  return(mat_comb)
}
mat_comb_train <- choose_cutoff(P_Eig, dat_num)
# mat_comb[which(mat_comb$accuracy == max(mat_comb$accuracy)),]
plot(cutoff,mat_comb_train[,1], type = "l",col=2, ylim = c(0,1))
lines(cutoff,mat_comb_train[,3],col="darkgreen",lwd=2)
lines(cutoff,mat_comb_train[,4],col="blue",lwd=2)
legend(0.5,0.55,col=c(2,"darkgreen","blue"),lwd=c(2,2,2),c("accuracy", "Specificity", "Sensitivity"), cex = 0.8, horiz=F, bty = "n")
#par(op)
abline(v= 0.473, lty=3, col = "black", lwd = 1.5)
# abline(v= 0.4611, lty=3, col = "orange")
dev.off()
####test set
P_Eig_test <- predict(output_Eigen,newdata = data.frame(EigenPhred = dat_num_test$EigenPhred), type = "response")
#confusion matrix test
predicted_values_test <- ifelse(P_Eig_test >= mean(P_Eig), 1 , 0)
conf_matrix<-table(predicted_values_test,dat_num_test[,1])
accuracy<-(conf_matrix[1,1]+conf_matrix[2,2])/(sum(conf_matrix))

mat_comb_test <- choose_cutoff(P_Eig_test, dat_num_test)
# mat_comb[which(mat_comb$accuracy == max(mat_comb$accuracy)),]
plot(cutoff,mat_comb_test[,1], type = "l",col=2, ylim = c(0,1))
lines(cutoff,mat_comb_test[,3],col="darkgreen",lwd=2)
lines(cutoff,mat_comb_test[,4],col="blue",lwd=2)
legend(0.5,0.55,col=c(2,"darkgreen","blue"),lwd=c(2,2,2),c("accuracy", "Specificity", "Sensitivity"), cex = 0.8, horiz=F, bty = "n")
#par(op)
abline(v= 0.473, lty=3, col = "black")
#  abline(v= 0.4611, lty=3, col = "orange")
