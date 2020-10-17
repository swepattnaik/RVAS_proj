##Non-coding variants
##filtered set
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

##This function currently works for the HAIL processed variant annotation output

library(g3viz)

fil_tab <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS/round2/all_isksmgrb_combset2020_variants_filt_all_fields_rmdup.tsv",
                      sep = "\t", stringsAsFactors = F, header = T)
fil_tab <- fil_tab[!is.na(fil_tab$SAMPLE),]
dim(fil_tab)
#fil_tab_isks_set <- fil_tab[fil_tab$is_case == 1,]
gen_g3viz_inp <- function(var_file, gene){
  fil_tab_gene <- var_file[var_file$gene_symbol %in% gene,]
  
  fil_tab_gene_fmt <- fil_tab_gene[,c(9, 131, 2, 4, 11, 126, 17, 1, 129)]
  
  var_cols <- do.call("rbind.data.frame", strsplit(as.character(fil_tab_gene_fmt$VARIANT), ":")) 
  colnames(var_cols) <- c("Chromosome", "Start_Position", 
                          "Reference_Allele", "Tumor_Seq_Allele2")                   
  
  var_cols$End_position <-  var_cols$Start_Position            
  var_cols$Tumor_Seq_Allele1 <- var_cols$Reference_Allele
  var_cols$Chromosome <- gsub("^", "chr", var_cols$Chromosome)
  var_cols <- var_cols[,c(1:2,5,3,6,4)]
  
  fil_tab_gene_fmt$vep_hgvsp <- gsub("^.*:", "", fil_tab_gene_fmt$vep_hgvsp)
  AA_pos <- strsplit(fil_tab_gene_fmt$vep_hgvsp, split = "\\_|fs")
  fil_tab_gene_fmt$AA_Position <- unlist(lapply(AA_pos, function(x) as.numeric(gsub("[^0-9]", "",x[1]))))
  
  
  fil_tab_gene_fmt <- fil_tab_gene_fmt[,-c(3)]
  fil_tab_gene_fmt$HGVSp_Short <- fil_tab_gene_fmt$vep_hgvsp
  colnames(fil_tab_gene_fmt) <- c("Hugo_Symbol", "cohort", "VAF", "Variant_Classification", "comb_score",
                                  "HGVSp", "Sample", "MAF_coh","AA_Position", "HGVSp_Short")
  #fil_tab_gene_fmt$HGVSp <- gsub("^.*_\\(", "",fil_tab_gene_fmt$HGVSp)
  #fil_tab_gene_fmt$HGVSp <- gsub("\\)$", "",fil_tab_gene_fmt$HGVSp)
  ##combine dataframes
  fil_tab_gene_fin <- cbind.data.frame(fil_tab_gene_fmt, var_cols)
  fil_tab_gene_fin$Variant_Type <- ifelse(fil_tab_gene_fin$Variant_Classification == "missense_variant",
                                          "SNV", "Other")
  fil_tab_gene_fin$Mutation_Class <- fil_tab_gene_fmt$Variant_Classification
  

  return(fil_tab_gene_fin)
}


##Plot function

plot_lollipop <- function(inp_fil, gene, covariate){
#  chart.options <- g3Lollipop.theme(theme.name = "ggplot2",
#                                    title.text = gene) 
 # gen_g3viz_inp(var_file = fil_tab, gene = gene)
  #   g3Lollipop.options(highlight.text.angle = "120")
  
 # if(covariate == "EigenPass"){ ##specific to eigenphred score
 #   i_fil <- readRDS("~/RVAS/shard_sub_tier3/DT_sheet/filt3_var_CIM12Sep2018.rds")
 #   c_score <- i_fil[i_fil$gene_symbol %in% gene,]$comb_score
    inp_fil1 <- gen_g3viz_inp(inp_fil, gene = gene)
 #   inp_fil1$comb_score <- c_score
    inp_fil1 <- inp_fil1[inp_fil1$comb_score >= 5.6,]
 #   inp_fil$EigenPass <- ifelse(inp_fil$Start_Position %in% inp_fil1$Start_Position, "Pass", "Fail")
    #  inp_fil$comb_score <- c_score
    #  inp_fil$EigenPass <- ifelse(inp_fil$comb_score >= 5, "Pass", "Fail")
 # }
  
  g3Lollipop(inp_fil1,
             gene.symbol = gene,
             plot.options = g3Lollipop.options( chart.type = "pie", 
                                                lollipop.color.scheme = "pie2",
                                                highlight.text.angle = "60",
                                                title.text = gene),
             btn.style = "blue",
             factor.col = covariate, save.png.btn = T, save.svg.btn = F,
             output.filename = paste(gene, "_lolli", ".png", sep = "")) 
}

#tp53_fin <- gen_g3viz_inp(var_file = fil_tab, gene = "TP53")
#tp53_fin$EigenPhred <- ifelse(is.na(tp53_fin$EigenPhred), 20, tp53_fin$EigenPhred)


plot_lollipop(inp_fil = fil_tab, gene = "TP53", covariate = "Mutation_Class")
plot_lollipop(inp_fil = fil_tab, gene = "TP53", covariate = "cohort")

chart.options = g3Lollipop.theme(
  theme.name = "default",
  title.text = "default theme title",
  y.axis.label = "y-label",
  legend.title = "legend-title")

g3Lollipop(tp53_fin,
           plot.options = chart.options,
           btn.style = "blue",
           gene.symbol = "TP53")


##genetrack
library(trackViewer)
library(rtracklayer)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(VariantAnnotation)
library(httr)

getGeneIDsFromTxDb <- function(gr, txdb){
  stopifnot(is(gr, "GRanges"))
  stopifnot(length(gr)>0)
  stopifnot(is(txdb, "TxDb"))
  if(length(gr)>1){
    warning("The length of gr is greater than 1. Only first genomic location will be used.")
    gr <- gr[1]
  }
  genes <- genes(txdb, columns="gene_id")
  genes <- subsetByOverlaps(genes, gr)
  return(genes$gene_id)
}
grW <- parse2GRanges("chr11:122,830,799-123,116,707")
ids <- getGeneIDsFromTxDb(grW, TxDb.Hsapiens.UCSC.hg19.knownGene)
symbols <- mget(ids, org.Hs.egSYMBOL)
genes <- geneTrack(ids, TxDb.Hsapiens.UCSC.hg19.knownGene, 
                   symbols, asList=FALSE)
optSty <- optimizeStyle(trackList(repA, fox2, genes), theme="safe")
trackListW <- optSty$tracks
viewerStyleW <- optSty$style
viewTracks(trackListW, gr=grW, viewerStyle=viewerStyleW)
