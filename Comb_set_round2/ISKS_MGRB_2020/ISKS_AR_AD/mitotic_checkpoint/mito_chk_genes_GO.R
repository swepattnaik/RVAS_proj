.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )

`%nin%` = Negate(`%in%`)

library(org.Hs.eg.db)
library(GO.db)

#go_id_1 = GOID( GOTERM[ Term(GOTERM) == "mitotic G2 DNA damage checkpoint"])
#go_id = GOID( GOTERM[ Term(GOTERM) == "mitotic G2/M transition checkpoint"])

go_id_1 = GOID( GOTERM[ Term(GOTERM) == "mitotic cytokinesis checkpoint"])
go_id_2 = GOID( GOTERM[ Term(GOTERM) == "mitotic DNA integrity checkpoint"])
go_id_3 = GOID( GOTERM[ Term(GOTERM) == "mitotic G1/S transition checkpoint"])
go_id_4 = GOID( GOTERM[ Term(GOTERM) == "mitotic G2/M transition checkpoint"])
go_id_5 = GOID( GOTERM[ Term(GOTERM) == "signal transduction involved in mitotic cell cycle checkpoint"])
go_id_6 = GOID( GOTERM[ Term(GOTERM) == "mitotic spindle checkpoint"])

#get(go_id, org.Hs.egGO2ALLEGS)

allegs = get(go_id_3, org.Hs.egGO2ALLEGS)

genes = unlist(mget(allegs,org.Hs.egSYMBOL))

#mitotic G2 DNA damage checkpoint

Go_terms <- c("mitotic cytokinesis checkpoint", "mitotic DNA integrity checkpoint", "mitotic G1/S transition checkpoint",
              "mitotic G2/M transition checkpoint", "signal transduction involved in mitotic cell cycle checkpoint",
              "mitotic spindle checkpoint")

mit_chkpoint_genes <- list()
for(i in 1:length(Go_terms)){
  g1 <- GOID( GOTERM[ Term(GOTERM) == Go_terms[i]])
  allegs = get(g1, org.Hs.egGO2ALLEGS)
  mit_chkpoint_genes[[i]] = unique(unlist(mget(allegs,org.Hs.egSYMBOL)))
  print(i)
                      
}

length(unique(unlist(mit_chkpoint_genes)))

mito_chk_pt_genes <- unique(unlist(mit_chkpoint_genes))

PPI_col <- read.delim("~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/comb_set_2020/ISKS_AR_AD/SKAT/scrambled_null_ppi_fisher_res_fil_final.tsv", 
                      sep = "\t", header = T, stringsAsFactors = F)
#Mitotic checkpoint genes preset in enriched SKATbin set
dput(mito_chk_pt_genes[mito_chk_pt_genes %in% PPI_col$gene])

write.table(mito_chk_pt_genes, "~/RVAS/shard_sub_tier3/DT_sheet/EXOME_isks_risc/test/git_proj_rvas/Comb_set_round2/ISKS_MGRB_2020/ISKS_AR_AD/mitotic_checkpoint/mitocheck_point.txt",
            sep = "", row.names = F, quote = F)
