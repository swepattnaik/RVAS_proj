##select samples for telomere PCR
##buffy coat as well as whole blood this will limit the samples we can use to those collected in Australia.
##PID number beginning with 012., 013., 015., 016., or 017.

t1 <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/PID_combset2020_CGC_skatBin_repstress_potint_predNFE_rmdup_maf_clueGO_Apr292020_June1.tsv", sep = "\t", header = T, stringsAsFactors = F)
t2 <- (t1[,c(1:3,76)])
id_mat <- c("^012.", "^013.", "^015.", "^016.", "^017.")
t2$hit <- ifelse(grepl(paste(id_mat,collapse="|"), t2$pid), 1, 0)
write.table(t2[t2$hit == 1,c(1:4)] , "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/PID_combset2020_June1_telo_length_buffy_wholeblood_AUS_cases.tsv",
sep = "\t", quote = F, row.names = F)
