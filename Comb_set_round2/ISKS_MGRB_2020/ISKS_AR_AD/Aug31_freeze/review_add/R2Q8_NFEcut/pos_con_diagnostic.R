##R2Q8; NFE <= 0.0002


cpx_OR_fisher_sarc <- function(case, control, case_coh_size, cont_coh_size, cpx_name){
  
  #ppi_res_tab <- ppi_res[ppi_res$gene %in% cpx_list,]
  # ppi_res_tab[,2] <- ifelse(ppi_res_tab[,2] == 0, 1, ppi_res_tab[,2])
  inp <- c(case, case_coh_size - case , 
           control, cont_coh_size - control)
  sim_mat <- matrix(inp ,nrow = 2, ncol = 2)
  colnames(sim_mat) <- c("case", "cont")
  rownames(sim_mat) <- c("hits", "no_hits")
  #ft <- fisher.test(sim_mat, alternative = "greater")
  #ft <- fisher.test(sim_mat)
  ft <- fisher.test(sim_mat, conf.int = T, conf.level = 0.99)
  #return(ft)
  ft_df <- cbind.data.frame("Region" = cpx_name ,"Cases" = sim_mat[1,1],
                            "Controls" = sim_mat[1,2],
                            "Fish_pval" = ft$p.value,"CI_lower" = ft$conf.int[1],
                            "CI_upper" = ft$conf.int[2],
                            "OR_Fish" = ft$estimate, "case_coh_size" = sum(sim_mat[,1]),
                            "control_coh_size" = sum(sim_mat[,2]))
  return(ft_df)
}
