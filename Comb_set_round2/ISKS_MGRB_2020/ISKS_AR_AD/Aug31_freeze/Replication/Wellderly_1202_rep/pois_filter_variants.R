##Merge tab and vep files containing variant from Shelterin, Centrosome, NF1 and HBOC genes
##map variant info to vep files
##filter multiple annotations using coding potential, transcript length etc
`%nin%` = Negate(`%in%`)
tab_fil <- read.delim("~/RVAS/Wellderly_2012/wellderly_1202_GT.tab", sep = "\t", header = T,
                      stringsAsFactors = F)
vep_fil <- read.delim("~/RVAS/Wellderly_2012/wellderly_1202.vep", sep = "\t", header = T,
                      stringsAsFactors = F)

tab_fil$sel_id <- paste(tab_fil$ID, tab_fil$SAMPLE, sep = "_")
vep_fil$sel_id <- paste(vep_fil$Uploaded_variation, vep_fil$SAMPLE, sep = "_")


##Append tab file
vep_fil[,70:82] <- tab_fil[match(vep_fil$sel_id, tab_fil$sel_id), 1:13]

##Filter vep based on gnomad_AF_NFE
vep_fil$gnomAD_NFE_AF <- as.numeric(vep_fil$gnomAD_NFE_AF)
vep_fil$gnomAD_NFE_AF <- ifelse(is.na(vep_fil$gnomAD_NFE_AF), 0, vep_fil$gnomAD_NFE_AF)
vep_fil_filt <- vep_fil[vep_fil$gnomAD_NFE_AF <= 0.0002,]
vep_fil_filt$gnomAD_AF <- as.numeric(vep_fil_filt$gnomAD_AF)
vep_fil_filt$gnomAD_AF <- ifelse(is.na(vep_fil_filt$gnomAD_AF), 0, vep_fil_filt$gnomAD_AF)
vep_fil_filt <- vep_fil_filt[vep_fil_filt$gnomAD_AF < 0.05,] ##condition for rare variant selection in discovery set

##Select canonical transcript protein coding transcript
vep_fil_filt1 <- vep_fil_filt[vep_fil_filt$CANONICAL %in% "YES",]

##Filter vep based on coding potential; remove non-coding variants
toMatch <- c("3_prime_UTR_variant", "5_prime_UTR_variant", "intron_variant", 
             "synonymous_variant", "upstream_gene_variant", "non_coding_transcript_exon_variant",
             "non_coding_transcript_variant", "downstream_gene_variant", "NMD_transcript_variant",
             "intergenic_variant", "mature_miRNA_variant", "inframe_deletion", "inframe_insertion")

vep_fil_filt2 <- vep_fil_filt1[!(grepl(paste(toMatch,collapse="|"), 
                                      vep_fil_filt1$Consequence)),]


vep_fil_filt2$AD <- as.numeric(gsub(",.*$", "", vep_fil_filt2$AD))
vep_fil_filt2$VAF <- ifelse(!is.na(vep_fil_filt2$DP), vep_fil_filt2$AD/vep_fil_filt2$DP, vep_fil_filt2$AD/vep_fil_filt2$DPI)

##Poisson variable-depth somatic filter (similar to Hartwig data)

pois_filt <- function(AD, DP){
  filt_var <- ifelse(ppois(AD, DP/2, lower.tail = TRUE) > 0.002 & AD > 0.3*DP, "PASS", "FAIL")
  return(filt_var)
}

vep_fil_filt2$filt_pois <- ifelse(!is.na(vep_fil_filt2$DP),pois_filt(vep_fil_filt2$AD, vep_fil_filt2$DP),
                                  pois_filt(vep_fil_filt2$AD, vep_fil_filt2$DPI))
vep_fil_filt3 <- vep_fil_filt2[vep_fil_filt2$filt_pois %in% "PASS",]


##conservative cutoff for genotype quality using GQX; which is the lower bound of GQ score.
##setting GQX to 40 is ideal?
vep_fil_filt4 <- vep_fil_filt3[vep_fil_filt3$GQX > 40,]

##Selection of longest coding transcript
##populate accurate protein length by searching for protein length from gene_symbol column..
prot_len <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/UCSC_tables_GRCh37_RefSeq_genes_20190716_gene_length_for_canonical_and_longest_transcripts.txt",
                       sep = "\t", header = T, stringsAsFactors = F)
vep_fil_filt4$Protein_canonical_length <- prot_len[match(vep_fil_filt4$SYMBOL, prot_len$gene), 2]
vep_fil_filt4$last_five_perc <- ifelse((as.numeric(vep_fil_filt4$Protein_position)/vep_fil_filt4$Protein_canonical_length) <= 0.05, 1, 0)

##retain only C4/C5
toMatch1 <- c("Likely_pathogenic", "Pathogenic", "Pathogenic/Likely_pathogenic", 
              "Pathogenic/Likely_pathogenic/drug_response", "Pathogenic/Likely_pathogenic/other",
              "Pathogenic/Affects", "Pathogenic/risk_factor", "Pathogenic/Likely_pathogenic/risk_factor",
              "risk_factor")
toMatch2 <- c("frameshift_variant", "start_lost", "stop_gained")
vep_fil_filt4_C4 <- vep_fil_filt4[(grepl(paste(toMatch2,collapse="|"), 
                                         vep_fil_filt4$Consequence)),]
vep_fil_filt4_C5 <- vep_fil_filt4[(grepl(paste(toMatch1,collapse="|"), 
                                         vep_fil_filt4$CLIN_SIG)),]
vep_fil_filt4_C45 <- rbind.data.frame(vep_fil_filt4_C4, vep_fil_filt4_C5)
vep_fil_filt4_C45 <- unique(vep_fil_filt4_C45)
vep_fil_filt4_C45[,c(2,7,18:19,54,57,59,73:81)]
vep_fil_filt4_C45_fin <- vep_fil_filt4_C45
#vep_fil_filt4_C45_fin <- vep_fil_filt4_C45_fin[vep_fil_filt4_C45_fin$GQX > 40,]
write.table(vep_fil_filt4_C45_fin, "~/RVAS/Wellderly_2012/Wellderly_vep_pois_filt4_C45_GQX40.tsv", sep = "\t",
            row.names = F, quote = F)
##inp: with CHIP correction: remove 84 cases
inp = c(7, 832, 2, 1200)
sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
colnames(sim_mat) <- c("case", "cont")
rownames(sim_mat) <- c("hits", "no_hits")
fisher.test(sim_mat, conf.int = T, conf.level = 0.95)
