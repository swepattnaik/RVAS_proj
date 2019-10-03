.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

library(SKAT)
 data(SKAT.haplotypes)
 attach(SKAT.haplotypes)
 set.seed(500)
out1000.b <- Power_Logistic(Haplotype, SNPInfo$CHROM_POS, SubRegion.Length=5000, 
                        Prevalence = 0.000025,  Causal.Percent= 20, N.Sim=1000, 
                        MaxOR=7, Negative.Percent=20)

##SKAT function to plot power analysis output
#Get_RequiredSampleSize(out1000.b, Power=0.8)

##plot it
library(ggplot2)
library(tidyr)
out.b <- out1000.b
df_power <- as.data.frame(out.b$Power, stringsAsFactors = F)
df_power$samp_size <- as.numeric(rownames(df_power))
df_power_long <- gather(df_power, alpha_sig, power, 1:3, factor_key=TRUE)
df_power_long$alpha_sig <- as.factor(df_power_long$alpha_sig)
p <- ggplot(df_power_long, aes(x=samp_size, y = power)) + geom_line(aes(colour = alpha_sig))

p <- p + theme(legend.position = "bottom",
          legend.direction = "horizontal")

png("~/RVAS/shard_sub_tier3/DT_sheet/Power_calc/power_analysis_outb_1000sim.png", bg = "transparent")
p
dev.off()


