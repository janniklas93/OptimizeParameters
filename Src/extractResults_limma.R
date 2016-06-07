# abbreviations
# li = Limma
# g = GCRMA, m = MAS, l = LESN2, q = quantile, s = scaling, mp = median.polish, t = tukey.biweight, lm = log.median, rlm = rlm
# p005 = p_value 0.05, p001 = p-value 0.01
# l05 = lfc 0.5, l1 = lfc 1, l15 = lfc 1.5
targetPathway = as.character(target_pathways$target[target_pathways$DataSet %in% dataSet])

all_ranks_li = c()
all_p_li = c()
all_q_li = c()
for(bm in 1:length(back_methods)){
  for(nm in 1:length(normalize_methods)){
    for(sm in 1:length(summary_methods)){
      for(pc in 1:length(p_change)){
        for(lc in 1:length(lfc_change)){
          current = read.table(paste(outputPath, paste(back_methods[bm], normalize_methods[nm], summary_methods[sm], sep = "_"), "limma_Results", paste(paste("pVal", str_replace(as.character(p_change[pc]), "\\.", "_"), sep = ""), paste("lFc", str_replace(as.character(lfc_change[lc]), "\\.", "_"), sep = ""), sep = "_"), "diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
          current_rank = which(current$GS %in% targetPathway)
          all_ranks_li = c(all_ranks_li, current_rank)
          all_p_li     = c(all_p_li, current$p_value[current_rank])
          all_q_li     = c(all_q_li, current$q_value[current_rank])
        }
      }
    }
  }
}

parameterSet_li = c("m_q_mp_p005_l05", "m_q_mp_p005_l1", "m_q_mp_p005_l15", "m_q_mp_p001_l05", "m_q_mp_p001_l1", "m_q_mp_p001_l15",
                    "m_q_t_p005_l05", "m_q_t_p005_l1", "m_q_t_p005_l15", "m_q_t_p001_l05", "m_q_t_p001_l1", "m_q_t_p001_l15",
                    "m_q_lm_p005_l05", "m_q_lm_p005_l1", "m_q_lm_p005_l15", "m_q_lm_p001_l05", "m_q_lm_p001_l1", "m_q_lm_p001_l15",
                    "m_q_rlm_p005_l05", "m_q_rlm_p005_l1", "m_q_rlm_p005_l15", "m_q_rlm_p001_l05", "m_q_rlm_p001_l1", "m_q_rlm_p001_l15",  
                    "m_s_mp_p005_l05", "m_s_mp_p005_l1", "m_s_mp_p005_l15", "m_s_mp_p001_l05", "m_s_mp_p001_l1", "m_s_mp_p001_l15",
                    "m_s_t_p005_l05", "m_s_t_p005_l1", "m_s_t_p005_l15", "m_s_t_p001_l05", "m_s_t_p001_l1", "m_s_t_p001_l15",
                    "m_s_lm_p005_l05", "m_s_lm_p005_l1", "m_s_lm_p005_l15", "m_s_lm_p001_l05", "m_s_lm_p001_l1", "m_s_lm_p001_l15",
                    "m_s_rlm_p005_l05", "m_s_rlm_p005_l1", "m_s_rlm_p005_l15", "m_s_rlm_p001_l05", "m_s_rlm_p001_l1", "m_s_rlm_p001_l15",
                    "g_q_mp_p005_l05", "g_q_mp_p005_l1", "g_q_mp_p005_l15", "g_q_mp_p001_l05", "g_q_mp_p001_l1", "g_q_mp_p001_l15",
                    "g_q_t_p005_l05", "g_q_t_p005_l1", "g_q_t_p005_l15", "g_q_t_p001_l05", "g_q_t_p001_l1", "g_q_t_p001_l15",
                    "g_q_lm_p005_l05", "g_q_lm_p005_l1", "g_q_lm_p005_l15", "g_q_lm_p001_l05", "g_q_lm_p001_l1", "g_q_lm_p001_l15",
                    "g_q_rlm_p005_l05", "g_q_rlm_p005_l1", "g_q_rlm_p005_l15", "g_q_rlm_p001_l05", "g_q_rlm_p001_l1", "g_q_rlm_p001_l15",  
                    "g_s_mp_p005_l05", "g_s_mp_p005_l1", "g_s_mp_p005_l15", "g_s_mp_p001_l05", "g_s_mp_p001_l1", "g_s_mp_p001_l15",
                    "g_s_t_p005_l05", "g_s_t_p005_l1", "g_s_t_p005_l15", "g_s_t_p001_l05", "g_s_t_p001_l1", "g_s_t_p001_l15",
                    "g_s_lm_p005_l05", "g_s_lm_p005_l1", "g_s_lm_p005_l15", "g_s_lm_p001_l05", "g_s_lm_p001_l1", "g_s_lm_p001_l15",
                    "g_s_rlm_p005_l05", "g_s_rlm_p005_l1", "g_s_rlm_p005_l15", "g_s_rlm_p001_l05", "g_s_rlm_p001_l1", "g_s_rlm_p001_l15",
                    "l_q_mp_p005_l05", "l_q_mp_p005_l1", "l_q_mp_p005_l15", "l_q_mp_p001_l05", "l_q_mp_p001_l1", "l_q_mp_p001_l15",
                    "l_q_t_p005_l05", "l_q_t_p005_l1", "l_q_t_p005_l15", "l_q_t_p001_l05", "l_q_t_p001_l1", "l_q_t_p001_l15",
                    "l_q_lm_p005_l05", "l_q_lm_p005_l1", "l_q_lm_p005_l15", "l_q_lm_p001_l05", "l_q_lm_p001_l1", "l_q_lm_p001_l15",
                    "l_q_rlm_p005_l05", "l_q_rlm_p005_l1", "l_q_rlm_p005_l15", "l_q_rlm_p001_l05", "l_q_rlm_p001_l1", "l_q_rlm_p001_l15",  
                    "l_s_mp_p005_l05", "l_s_mp_p005_l1", "l_s_mp_p005_l15", "l_s_mp_p001_l05", "l_s_mp_p001_l1", "l_s_mp_p001_l15",
                    "l_s_t_p005_l05", "l_s_t_p005_l1", "l_s_t_p005_l15", "l_s_t_p001_l05", "l_s_t_p001_l1", "l_s_t_p001_l15",
                    "l_s_lm_p005_l05", "l_s_lm_p005_l1", "l_s_lm_p005_l15", "l_s_lm_p001_l05", "l_s_lm_p001_l1", "l_s_lm_p001_l15",
                    "l_s_rlm_p005_l05", "l_s_rlm_p005_l1", "l_s_rlm_p005_l15", "l_s_rlm_p001_l05", "l_s_rlm_p001_l1", "l_s_rlm_p001_l15")


result_out_li = data.frame("parameter_setting" = parameterSet_li, "rank" = all_ranks_li, "p_value" = all_p_li, "q_value" = all_q_li)
result_out_li = result_out_li[order(result_out_li$rank, result_out_li$q_value, decreasing = FALSE), ]
resultFile_li = paste(outputPath, paste("overall_results_", dataSet, "_limma", ".txt", sep = ""), sep = "/")
write.table(result_out_li, file = resultFile_li, sep = "\t", quote = FALSE, row.names = FALSE)

# write best/worst ranks, p.values and p-values for limma to file
best_rank_li = min(all_ranks_li)
best_rank_li_setting = which(all_ranks_li %in% best_rank_li)
best_rank_li_setting = parameterSet_li[best_rank_li_setting]
worst_rank_li = max(all_ranks_li)
worst_rank_li_setting = which(all_ranks_li %in% worst_rank_li)
worst_rank_li_setting = parameterSet_li[worst_rank_li_setting]

best_p_li = min(all_p_li)
best_p_li_setting = which(all_p_li %in% best_p_li)
best_p_li_setting = parameterSet_li[best_p_li_setting]
worst_p_li = max(all_p_li)
worst_p_li_setting = which(all_p_li %in% worst_p_li)
worst_p_li_setting = parameterSet_li[worst_p_li_setting]

best_q_li = min(all_q_li)
best_q_li_setting = which(all_q_li %in% best_q_li)
best_q_li_setting = parameterSet_li[best_q_li_setting]
worst_q_li = max(all_q_li)
worst_q_li_setting = which(all_q_li %in% worst_q_li)
worst_q_li_setting = parameterSet_li[worst_q_li_setting]

resultFile_best_li = paste(outputPath, paste("best_worst_results_", dataSet, "_limma", ".txt", sep = ""), sep = "/")
write(paste(paste("best rank: ", best_rank_li, sep = "  "), paste("parameter setting: ", best_rank_li_setting, sep = "  "), sep = "  "), file = resultFile_best_li, append = TRUE)
write(paste(paste("worst rank: ", worst_rank_li, sep = "  "), paste("parameter setting: ", worst_rank_li_setting, sep = "  "), sep = "  "), file = resultFile_best_li, append = TRUE)
write(paste(paste("best p-value: ", best_p_li, sep = "  "), paste("parameter setting: ", best_p_li_setting, sep = "  "), sep = "  "), file = resultFile_best_li, append = TRUE)
write(paste(paste("worst p-value: ", worst_p_li, sep = "  "), paste("parameter setting: ", worst_p_li_setting, sep = "  "), sep = "  "), file = resultFile_best_li, append = TRUE)
write(paste(paste("best q-value: ", best_q_li, sep = "  "), paste("parameter setting: ", best_q_li_setting, sep = "  "), sep = "  "), file = resultFile_best_li, append = TRUE)
write(paste(paste("worst q-value: ", worst_q_li, sep = "  "), paste("parameter setting: ", worst_q_li_setting, sep = "  "), sep = "  "), file = resultFile_best_li, append = TRUE)

if(FALSE){
######## Limma results ######### -------------------------------------------------------------------------------------------------------------------------------------------------------
# abbreviations
# li = Limma, gsea = GSEA
# g = GCRMA, m = MAS, l = LESN2, q = quantile, s = scaling, mp = median.polish, t = tukey.biweight, lm = log.median, rlm = rlm
#### GCRMA quantile median.polish --------------------------------------------------
li_g_q_mp_0.05  = read.table(paste(outputPath, "GCRMA_quantile_median.polish/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_q_mp_0.1   = read.table(paste(outputPath, "GCRMA_quantile_median.polish/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_q_mp_0.15  = read.table(paste(outputPath, "GCRMA_quantile_median.polish/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_q_mp_0.2   = read.table(paste(outputPath, "GCRMA_quantile_median.polish/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### GCRMA quantile tukey.biweight
li_g_q_t_0.05  = read.table(paste(outputPath, "GCRMA_quantile_tukey.biweight/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_q_t_0.1   = read.table(paste(outputPath, "GCRMA_quantile_tukey.biweight/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_q_t_0.15  = read.table(paste(outputPath, "GCRMA_quantile_tukey.biweight/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_q_t_0.2   = read.table(paste(outputPath, "GCRMA_quantile_tukey.biweight/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### GCRMA quantile log.median 
li_g_q_lm_0.05 = read.table(paste(outputPath, "GCRMA_quantile_log.median/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_q_lm_0.1  = read.table(paste(outputPath, "GCRMA_quantile_log.median/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_q_lm_0.15 = read.table(paste(outputPath, "GCRMA_quantile_log.median/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_q_lm_0.2  = read.table(paste(outputPath, "GCRMA_quantile_log.median/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### GCRMA quantile rlm
li_g_q_rlm_0.05 = read.table(paste(outputPath, "GCRMA_quantile_rlm/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_q_rlm_0.1  = read.table(paste(outputPath, "GCRMA_quantile_rlm/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_q_rlm_0.15 = read.table(paste(outputPath, "GCRMA_quantile_rlm/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_q_rlm_0.2  = read.table(paste(outputPath, "GCRMA_quantile_rlm/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### GCRMA scaling median.polish
li_g_s_mp_0.05   = read.table(paste(outputPath, "GCRMA_scaling_median.polish/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_s_mp_0.1    = read.table(paste(outputPath, "GCRMA_scaling_median.polish/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_s_mp_0.15   = read.table(paste(outputPath, "GCRMA_scaling_median.polish/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_s_mp_0.2    = read.table(paste(outputPath, "GCRMA_scaling_median.polish/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### GCRMA scaling tukey.biweight
li_g_s_t_0.05   = read.table(paste(outputPath, "GCRMA_scaling_tukey.biweight/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_s_t_0.1    = read.table(paste(outputPath, "GCRMA_scaling_tukey.biweight/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_s_t_0.15   = read.table(paste(outputPath, "GCRMA_scaling_tukey.biweight/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_s_t_0.2    = read.table(paste(outputPath, "GCRMA_scaling_tukey.biweight/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### GCRMA scaling log.median 
li_g_s_lm_0.05  = read.table(paste(outputPath, "GCRMA_scaling_log.median/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_s_lm_0.1   = read.table(paste(outputPath, "GCRMA_scaling_log.median/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_s_lm_0.15  = read.table(paste(outputPath, "GCRMA_scaling_log.median/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_s_lm_0.2   = read.table(paste(outputPath, "GCRMA_scaling_log.median/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### GCRMA scaling rlm
li_g_s_rlm_0.05 = read.table(paste(outputPath, "GCRMA_scaling_rlm/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_s_rlm_0.1 = read.table(paste(outputPath, "GCRMA_scaling_rlm/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_s_rlm_0.15 = read.table(paste(outputPath, "GCRMA_scaling_rlm/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_g_s_rlm_0.2 = read.table(paste(outputPath, "GCRMA_scaling_rlm/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### MAS quantile median.polish --------------------------------------------
li_m_q_mp_0.05  = read.table(paste(outputPath, "MAS_quantile_median.polish/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_q_mp_0.1   = read.table(paste(outputPath, "MAS_quantile_median.polish/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_q_mp_0.15  = read.table(paste(outputPath, "MAS_quantile_median.polish/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_q_mp_0.2   = read.table(paste(outputPath, "MAS_quantile_median.polish/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### MAS quantile tukey.biweight
li_m_q_t_0.05  = read.table(paste(outputPath, "MAS_quantile_tukey.biweight/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_q_t_0.1   = read.table(paste(outputPath, "MAS_quantile_tukey.biweight/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_q_t_0.15  = read.table(paste(outputPath, "MAS_quantile_tukey.biweight/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_q_t_0.2   = read.table(paste(outputPath, "MAS_quantile_tukey.biweight/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### MAS quantile log.median 
li_m_q_lm_0.05 = read.table(paste(outputPath, "MAS_quantile_log.median/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_q_lm_0.1  = read.table(paste(outputPath, "MAS_quantile_log.median/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_q_lm_0.15 = read.table(paste(outputPath, "MAS_quantile_log.median/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_q_lm_0.2  = read.table(paste(outputPath, "MAS_quantile_log.median/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### MAS quantile rlm
li_m_q_rlm_0.05 = read.table(paste(outputPath, "MAS_quantile_rlm/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_q_rlm_0.1  = read.table(paste(outputPath, "MAS_quantile_rlm/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_q_rlm_0.15 = read.table(paste(outputPath, "MAS_quantile_rlm/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_q_rlm_0.2  = read.table(paste(outputPath, "MAS_quantile_rlm/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### MAS scaling median.polish
li_m_s_mp_0.05   = read.table(paste(outputPath, "MAS_scaling_median.polish/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_mp_0.1    = read.table(paste(outputPath, "MAS_scaling_median.polish/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_mp_0.15   = read.table(paste(outputPath, "MAS_scaling_median.polish/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_mp_0.2    = read.table(paste(outputPath, "MAS_scaling_median.polish/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### MAS scaling tukey.biweight
li_m_s_t_0.05   = read.table(paste(outputPath, "MAS_scaling_tukey.biweight/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_t_0.1    = read.table(paste(outputPath, "MAS_scaling_tukey.biweight/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_t_0.15   = read.table(paste(outputPath, "MAS_scaling_tukey.biweight/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_t_0.2    = read.table(paste(outputPath, "MAS_scaling_tukey.biweight/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### MAS scaling log.median 
li_m_s_lm_0.05  = read.table(paste(outputPath, "MAS_scaling_log.median/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_lm_0.1   = read.table(paste(outputPath, "MAS_scaling_log.median/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_lm_0.15  = read.table(paste(outputPath, "MAS_scaling_log.median/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_lm_0.2   = read.table(paste(outputPath, "MAS_scaling_log.median/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### MAS scaling rlm
li_m_s_rlm_0.05 = read.table(paste(outputPath, "MAS_scaling_rlm/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_rlm_0.1  = read.table(paste(outputPath, "MAS_scaling_rlm/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_rlm_0.15 = read.table(paste(outputPath, "MAS_scaling_rlm/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_rlm_0.2  = read.table(paste(outputPath, "MAS_scaling_rlm/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### LESN2 quantile median.polish --------------------------------------------------
li_l_q_mp_0.05  = read.table(paste(outputPath, "LESN2_quantile_median.polish/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_q_mp_0.1   = read.table(paste(outputPath, "LESN2_quantile_median.polish/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_q_mp_0.15  = read.table(paste(outputPath, "LESN2_quantile_median.polish/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_q_mp_0.2   = read.table(paste(outputPath, "LESN2_quantile_median.polish/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### LESN2 quantile tukey.biweight
li_l_q_t_0.05  = read.table(paste(outputPath, "LESN2_quantile_tukey.biweight/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_q_t_0.1   = read.table(paste(outputPath, "LESN2_quantile_tukey.biweight/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_q_t_0.15  = read.table(paste(outputPath, "LESN2_quantile_tukey.biweight/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_q_t_0.2   = read.table(paste(outputPath, "LESN2_quantile_tukey.biweight/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### LESN2 quantile log.median 
li_l_q_lm_0.05 = read.table(paste(outputPath, "LESN2_quantile_log.median/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_q_lm_0.1  = read.table(paste(outputPath, "LESN2_quantile_log.median/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_q_lm_0.15 = read.table(paste(outputPath, "LESN2_quantile_log.median/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_q_lm_0.2  = read.table(paste(outputPath, "LESN2_quantile_log.median/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### LESN2 quantile rlm
li_l_q_rlm_0.05 = read.table(paste(outputPath, "LESN2_quantile_rlm/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_q_rlm_0.1  = read.table(paste(outputPath, "LESN2_quantile_rlm/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_q_rlm_0.15 = read.table(paste(outputPath, "LESN2_quantile_rlm/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_q_rlm_0.2  = read.table(paste(outputPath, "LESN2_quantile_rlm/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### LESN2 scaling median.polish
li_l_s_mp_0.05   = read.table(paste(outputPath, "LESN2_scaling_median.polish/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_s_mp_0.1    = read.table(paste(outputPath, "LESN2_scaling_median.polish/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_s_mp_0.15   = read.table(paste(outputPath, "LESN2_scaling_median.polish/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_s_mp_0.2    = read.table(paste(outputPath, "LESN2_scaling_median.polish/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### LESN2 scaling tukey.biweight
li_l_s_t_0.05   = read.table(paste(outputPath, "LESN2_scaling_tukey.biweight/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_s_t_0.1    = read.table(paste(outputPath, "LESN2_scaling_tukey.biweight/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_s_t_0.15   = read.table(paste(outputPath, "LESN2_scaling_tukey.biweight/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_s_t_0.2    = read.table(paste(outputPath, "LESN2_scaling_tukey.biweight/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### LESN2 scaling log.median 
li_l_s_lm_0.05  = read.table(paste(outputPath, "LESN2_scaling_log.median/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_s_lm_0.1   = read.table(paste(outputPath, "LESN2_scaling_log.median/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_s_lm_0.15  = read.table(paste(outputPath, "LESN2_scaling_log.median/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_l_s_lm_0.2   = read.table(paste(outputPath, "LESN2_scaling_log.median/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
#### LESN2 scaling rlm
li_m_s_rlm_0.05 = read.table(paste(outputPath, "LESN2_scaling_rlm/limma_Results/0_05/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_rlm_0.1  = read.table(paste(outputPath, "LESN2_scaling_rlm/limma_Results/0_1/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_rlm_0.15 = read.table(paste(outputPath, "LESN2_scaling_rlm/limma_Results/0_15/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
li_m_s_rlm_0.2  = read.table(paste(outputPath, "LESN2_scaling_rlm/limma_Results/0_2/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")

#### GSEA results --------------
#result_gsea_GCRMA_Normal  = read.table(paste(outputPath, "GCRMA/gseaResults", paste(dataSet, "_GCRMA.SUMMARY.RESULTS.REPORT.Normal.txt", sep = ""), sep = "/"), sep = "\t", header = TRUE)
#result_gsea_GCRMA_Disease = read.table(paste(outputPath, "GCRMA/gseaResults", paste(dataSet, "_GCRMA.SUMMARY.RESULTS.REPORT.Disease.txt", sep = ""), sep = "/"), sep = "\t", header = TRUE)
#result_gsea_MAS_Normal    = read.table(paste(outputPath, "MAS/gseaResults", paste(dataSet, "_MAS.SUMMARY.RESULTS.REPORT.Normal.txt", sep = ""), sep = "/"), sep = "\t", header = TRUE)
#result_gsea_MAS_Disease   = read.table(paste(outputPath, "MAS/gseaResults", paste(dataSet, "_MAS.SUMMARY.RESULTS.REPORT.Disease.txt", sep = ""), sep = "/"), sep = "\t", header = TRUE)
#result_gsea_LESN2_Normal  = read.table(paste(outputPath, "LESN2/gseaResults", paste(dataSet, "_LESN2.SUMMARY.RESULTS.REPORT.Normal.txt", sep = ""), sep = "/"), sep = "\t", header = TRUE)
#result_gsea_LESN2_Disease = read.table(paste(outputPath, "LESN2/gseaResults", paste(dataSet, "_LESN2.SUMMARY.RESULTS.REPORT.Normal.txt", sep = ""), sep = "/"), sep = "\t", header = TRUE)

#result_gsea_GCRMA = rbind(result_gsea_GCRMA_Normal, result_gsea_GCRMA_Disease)
#result_gsea_GCRMA = result_gsea_GCRMA[order(result_gsea_GCRMA$FDR.q.val, result_gsea_GCRMA$NOM.p.val, decreasing = FALSE), ]
#result_gsea_MAS = rbind(result_gsea_MAS_Normal, result_gsea_MAS_Disease)
#result_gsea_MAS = result_gsea_MAS[order(result_gsea_MAS$FDR.q.val, result_gsea_MAS$NOM.p.val, decreasing = FALSE), ]
#result_gsea_LESN2 = rbind(result_gsea_LESN2_Normal, result_gsea_LESN2_Disease)
#result_gsea_LESN2 = result_gsea_GCRMA[order(result_gsea_GCRMA$FDR.q.val, result_gsea_LESN2$NOM.p.val, decreasing = FALSE), ]

rm(result_gsea_GCRMA_Normal, result_gsea_GCRMA_Disease, result_gsea_MAS_Normal, result_gsea_MAS_Disease)

####### Benchmark ######## --------------------------------------------------------------------------------------------------------------------------------------------------
targetPathway = as.character(target_pathways$target[target_pathways$DataSet %in% dataSet])

#### extract rank of target geneSet --------------------------------------------------------
# GCRMA quantile median.polish ----------------------------------
rank_li_g_q_mp_0.05 = which(li_g_q_mp_0.05$GS %in% targetPathway)
rank_li_g_q_mp_0.1 = which(li_g_q_mp_0.1$GS %in% targetPathway)
rank_li_g_q_mp_0.15 = which(li_g_q_mp_0.15$GS %in% targetPathway)
rank_li_g_q_mp_0.2 = which(li_g_q_mp_0.2$GS %in% targetPathway)
# GCRMA quantile tukey.biweight
rank_li_g_q_t_0.05 = which(li_g_q_t_0.05$GS %in% targetPathway)
rank_li_g_q_t_0.1 = which(li_g_q_t_0.1$GS %in% targetPathway)
rank_li_g_q_t_0.15 = which(li_g_q_t_0.15$GS %in% targetPathway)
rank_li_g_q_t_0.2 = which(li_g_q_t_0.2$GS %in% targetPathway)
# GCRMA quantile log.median
rank_li_g_q_lm_0.05 = which(li_g_q_lm_0.05$GS %in% targetPathway)
rank_li_g_q_lm_0.1  = which(li_g_q_lm_0.1$GS %in% targetPathway)
rank_li_g_q_lm_0.15 = which(li_g_q_lm_0.15$GS %in% targetPathway)
rank_li_g_q_lm_0.2  = which(li_g_q_lm_0.2$GS %in% targetPathway)
# GCRMA quantile rlm
rank_li_g_q_rlm_0.05 = which(li_g_q_rlm_0.05$GS %in% targetPathway)
rank_li_g_q_rlm_0.1 = which(li_g_q_rlm_0.1$GS %in% targetPathway)
rank_li_g_q_rlm_0.15 = which(li_g_q_rlm_0.15$GS %in% targetPathway)
rank_li_g_q_rlm_0.2 = which(li_g_q_rlm_0.2$GS %in% targetPathway)
# GCRMA scaling median.polish
rank_li_g_s_mp_0.05 = which(li_g_s_mp_0.05$GS %in% targetPathway)
rank_li_g_s_mp_0.1 = which(li_g_s_mp_0.1$GS %in% targetPathway)
rank_li_g_s_mp_0.15 = which(li_g_s_mp_0.15$GS %in% targetPathway)
rank_li_g_s_mp_0.2 = which(li_g_s_mp_0.2$GS %in% targetPathway)
# GCRMA scaling tukey.biweight
rank_li_g_s_t_0.05 = which(li_g_s_t_0.05$GS %in% targetPathway)
rank_li_g_s_t_0.1 = which(li_g_s_t_0.1$GS %in% targetPathway)
rank_li_g_s_t_0.15 = which(li_g_s_t_0.15$GS %in% targetPathway)
rank_li_g_s_t_0.2 = which(li_g_s_t_0.2$GS %in% targetPathway)
# GCRMA scaling log.median
rank_li_g_s_lm_0.05 = which(li_g_s_lm_0.05$GS %in% targetPathway)
rank_li_g_s_lm_0.1 = which(li_g_s_lm_0.1$GS %in% targetPathway)
rank_li_g_s_lm_0.15 = which(li_g_s_lm_0.15$GS %in% targetPathway)
rank_li_g_s_lm_0.2 = which(li_g_s_lm_0.2$GS %in% targetPathway)
# GCRMA scaling rlm
rank_li_g_s_rlm_0.05 = which(li_g_s_rlm_0.05$GS %in% targetPathway)
rank_li_g_s_rlm_0.05 = which(li_g_s_rlm_0.05$GS %in% targetPathway)
rank_li_g_s_rlm_0.05 = which(li_g_s_rlm_0.05$GS %in% targetPathway)
rank_li_g_s_rlm_0.05 = which(li_g_s_rlm_0.05$GS %in% targetPathway)
# MAS quantile median.polish --------------------------------------
rank_li_m_q_mp_0.05 = which(li_m_q_mp_0.05$GS %in% targetPathway)
rank_li_m_q_mp_0.1 = which(li_m_q_mp_0.1$GS %in% targetPathway)
rank_li_m_q_mp_0.15 = which(li_m_q_mp_0.15$GS %in% targetPathway)
rank_li_m_q_mp_0.2 = which(li_m_q_mp_0.2$GS %in% targetPathway)
# MAS quantile tukey.biweight
rank_li_m_q_t_0.05 = which(li_m_q_t_0.05$GS %in% targetPathway)
rank_li_m_q_t_0.1 = which(li_m_q_t_0.1$GS %in% targetPathway)
rank_li_m_q_t_0.15 = which(li_m_q_t_0.15$GS %in% targetPathway)
rank_li_m_q_t_0.2 = which(li_m_q_t_0.2$GS %in% targetPathway)
# MAS quantile log.median
rank_li_m_q_lm_0.05 = which(li_m_q_lm_0.05$GS %in% targetPathway)
rank_li_m_q_lm_0.1  = which(li_m_q_lm_0.1$GS %in% targetPathway)
rank_li_m_q_lm_0.15 = which(li_m_q_lm_0.15$GS %in% targetPathway)
rank_li_m_q_lm_0.2  = which(li_m_q_lm_0.2$GS %in% targetPathway)
# MAS quantile rlm
rank_li_m_q_rlm_0.05 = which(li_m_q_rlm_0.05$GS %in% targetPathway)
rank_li_m_q_rlm_0.1  = which(li_m_q_rlm_0.1$GS %in% targetPathway)
rank_li_m_q_rlm_0.15 = which(li_m_q_rlm_0.15$GS %in% targetPathway)
rank_li_m_q_rlm_0.2  = which(li_m_q_rlm_0.2$GS %in% targetPathway)
# MAS scaling median.polish
rank_li_m_s_mp_0.05 = which(li_m_s_mp_0.05$GS %in% targetPathway)
rank_li_m_s_mp_0.1  = which(li_m_s_mp_0.1$GS %in% targetPathway)
rank_li_m_s_mp_0.15 = which(li_m_s_mp_0.15$GS %in% targetPathway)
rank_li_m_s_mp_0.2  = which(li_m_s_mp_0.2$GS %in% targetPathway)
# MAS scaling tukey.biweight
rank_li_m_s_t_0.05 = which(li_m_s_t_0.05$GS %in% targetPathway)
rank_li_m_s_t_0.1 = which(li_m_s_t_0.1$GS %in% targetPathway)
rank_li_m_s_t_0.15 = which(li_m_s_t_0.15$GS %in% targetPathway)
rank_li_m_s_t_0.2 = which(li_m_s_t_0.2$GS %in% targetPathway)
# MAS scaling log.median
rank_li_m_s_lm_0.05 = which(li_m_s_lm_0.05$GS %in% targetPathway)
rank_li_m_s_lm_0.1 = which(li_m_s_lm_0.1$GS %in% targetPathway)
rank_li_m_s_lm_0.15 = which(li_m_s_lm_0.15$GS %in% targetPathway)
rank_li_m_s_lm_0.2 = which(li_m_s_lm_0.2$GS %in% targetPathway)
# MAS scaling rlm
rank_li_m_s_rlm_0.05 = which(li_m_s_rlm_0.05$GS %in% targetPathway)
rank_li_m_s_rlm_0.05 = which(li_m_s_rlm_0.05$GS %in% targetPathway)
rank_li_m_s_rlm_0.05 = which(li_m_s_rlm_0.05$GS %in% targetPathway)
rank_li_m_s_rlm_0.05 = which(li_m_s_rlm_0.05$GS %in% targetPathway)
# LESN2 quantile median.polish ------------------------------------
rank_li_l_q_mp_0.05 = which(li_l_q_mp_0.05$GS %in% targetPathway)
rank_li_l_q_mp_0.1 = which(li_l_q_mp_0.1$GS %in% targetPathway)
rank_li_l_q_mp_0.15 = which(li_l_q_mp_0.15$GS %in% targetPathway)
rank_li_l_q_mp_0.2 = which(li_l_q_mp_0.2$GS %in% targetPathway)
# LESN2 quantile tukey.biweight
rank_li_l_q_t_0.05 = which(li_l_q_t_0.05$GS %in% targetPathway)
rank_li_l_q_t_0.1 = which(li_l_q_t_0.1$GS %in% targetPathway)
rank_li_l_q_t_0.15 = which(li_l_q_t_0.15$GS %in% targetPathway)
rank_li_l_q_t_0.2 = which(li_l_q_t_0.2$GS %in% targetPathway)
# LESN2 quantile log.median
rank_li_l_q_lm_0.05 = which(li_l_q_lm_0.05$GS %in% targetPathway)
rank_li_l_q_lm_0.1  = which(li_l_q_lm_0.1$GS %in% targetPathway)
rank_li_l_q_lm_0.15 = which(li_l_q_lm_0.15$GS %in% targetPathway)
rank_li_l_q_lm_0.2  = which(li_l_q_lm_0.2$GS %in% targetPathway)
# LESN2 quantile rlm
rank_li_l_q_rlm_0.05 = which(li_l_q_rlm_0.05$GS %in% targetPathway)
rank_li_l_q_rlm_0.1  = which(li_l_q_rlm_0.1$GS %in% targetPathway)
rank_li_l_q_rlm_0.15 = which(li_l_q_rlm_0.15$GS %in% targetPathway)
rank_li_l_q_rlm_0.2  = which(li_l_q_rlm_0.2$GS %in% targetPathway)
# LESN2 scaling median.polish
rank_li_l_s_mp_0.05 = which(li_l_s_mp_0.05$GS %in% targetPathway)
rank_li_l_s_mp_0.1  = which(li_l_s_mp_0.1$GS %in% targetPathway)
rank_li_l_s_mp_0.15 = which(li_l_s_mp_0.15$GS %in% targetPathway)
rank_li_l_s_mp_0.2  = which(li_l_s_mp_0.2$GS %in% targetPathway)
# LESN2 scaling tukey.biweight
rank_li_m_s_t_0.05 = which(li_m_s_t_0.05$GS %in% targetPathway)
rank_li_m_s_t_0.1 = which(li_m_s_t_0.1$GS %in% targetPathway)
rank_li_m_s_t_0.15 = which(li_m_s_t_0.15$GS %in% targetPathway)
rank_li_m_s_t_0.2 = which(li_m_s_t_0.2$GS %in% targetPathway)
# LESN2 scaling log.median
rank_li_l_s_lm_0.05 = which(li_l_s_lm_0.05$GS %in% targetPathway)
rank_li_l_s_lm_0.1 = which(li_l_s_lm_0.1$GS %in% targetPathway)
rank_li_l_s_lm_0.15 = which(li_l_s_lm_0.15$GS %in% targetPathway)
rank_li_l_s_lm_0.2 = which(li_l_s_lm_0.2$GS %in% targetPathway)
# LESN2 scaling rlm
rank_li_l_s_rlm_0.05 = which(li_l_s_rlm_0.05$GS %in% targetPathway)
rank_li_l_s_rlm_0.05 = which(li_l_s_rlm_0.05$GS %in% targetPathway)
rank_li_l_s_rlm_0.05 = which(li_l_s_rlm_0.05$GS %in% targetPathway)
rank_li_l_s_rlm_0.05 = which(li_l_s_rlm_0.05$GS %in% targetPathway)


#### extract p-values of target gene sets ---------------------------------------------------------------
# p GCRMA quantile median.polish ----------------------------------
p_li_g_q_mp_0.05 = which(li_g_q_mp_0.05$p_value %in% targetPathway)
p_li_g_q_mp_0.1 = which(li_g_q_mp_0.1$p_value %in% targetPathway)
p_li_g_q_mp_0.15 = which(li_g_q_mp_0.15$p_value %in% targetPathway)
p_li_g_q_mp_0.2 = which(li_g_q_mp_0.2$p_value %in% targetPathway)
# p GCRMA quantile tukey.biweight
p_li_g_q_t_0.05 = which(li_g_q_t_0.05$p_value %in% targetPathway)
p_li_g_q_t_0.1 = which(li_g_q_t_0.1$p_value %in% targetPathway)
p_li_g_q_t_0.15 = which(li_g_q_t_0.15$p_value %in% targetPathway)
p_li_g_q_t_0.2 = which(li_g_q_t_0.2$p_value %in% targetPathway)
# p GCRMA quantile log.median
p_li_g_q_lm_0.05 = which(li_g_q_lm_0.05$p_value %in% targetPathway)
p_li_g_q_lm_0.1  = which(li_g_q_lm_0.1$p_value %in% targetPathway)
p_li_g_q_lm_0.15 = which(li_g_q_lm_0.15$p_value %in% targetPathway)
p_li_g_q_lm_0.2  = which(li_g_q_lm_0.2$p_value %in% targetPathway)
# p GCRMA quantile rlm
p_li_g_q_rlm_0.05 = which(li_g_q_rlm_0.05$p_value %in% targetPathway)
p_li_g_q_rlm_0.1 = which(li_g_q_rlm_0.1$p_value %in% targetPathway)
p_li_g_q_rlm_0.15 = which(li_g_q_rlm_0.15$p_value %in% targetPathway)
p_li_g_q_rlm_0.2 = which(li_g_q_rlm_0.2$p_value %in% targetPathway)
# p GCRMA scaling median.polish
p_li_g_s_mp_0.05 = which(li_g_s_mp_0.05$p_value %in% targetPathway)
p_li_g_s_mp_0.1 = which(li_g_s_mp_0.1$p_value %in% targetPathway)
p_li_g_s_mp_0.15 = which(li_g_s_mp_0.15$p_value %in% targetPathway)
p_li_g_s_mp_0.2 = which(li_g_s_mp_0.2$p_value %in% targetPathway)
# p GCRMA scaling tukey.biweight
p_li_g_s_t_0.05 = which(li_g_s_t_0.05$p_value %in% targetPathway)
p_li_g_s_t_0.1 = which(li_g_s_t_0.1$p_value %in% targetPathway)
p_li_g_s_t_0.15 = which(li_g_s_t_0.15$p_value %in% targetPathway)
p_li_g_s_t_0.2 = which(li_g_s_t_0.2$p_value %in% targetPathway)
# p GCRMA scaling log.median
p_li_g_s_lm_0.05 = which(li_g_s_lm_0.05$p_value %in% targetPathway)
p_li_g_s_lm_0.1 = which(li_g_s_lm_0.1$p_value %in% targetPathway)
p_li_g_s_lm_0.15 = which(li_g_s_lm_0.15$p_value %in% targetPathway)
p_li_g_s_lm_0.2 = which(li_g_s_lm_0.2$p_value %in% targetPathway)
# p GCRMA scaling rlm
p_li_g_s_rlm_0.05 = which(li_g_s_rlm_0.05$p_value %in% targetPathway)
p_li_g_s_rlm_0.05 = which(li_g_s_rlm_0.05$p_value %in% targetPathway)
p_li_g_s_rlm_0.05 = which(li_g_s_rlm_0.05$p_value %in% targetPathway)
p_li_g_s_rlm_0.05 = which(li_g_s_rlm_0.05$p_value %in% targetPathway)
# p MAS quantile median.polish --------------------------------------
p_li_m_q_mp_0.05 = which(li_m_q_mp_0.05$p_value %in% targetPathway)
p_li_m_q_mp_0.1 = which(li_m_q_mp_0.1$p_value %in% targetPathway)
p_li_m_q_mp_0.15 = which(li_m_q_mp_0.15$p_value %in% targetPathway)
p_li_m_q_mp_0.2 = which(li_m_q_mp_0.2$p_value %in% targetPathway)
# p MAS quantile tukey.biweight
p_li_m_q_t_0.05 = which(li_m_q_t_0.05$p_value %in% targetPathway)
p_li_m_q_t_0.1 = which(li_m_q_t_0.1$p_value %in% targetPathway)
p_li_m_q_t_0.15 = which(li_m_q_t_0.15$p_value %in% targetPathway)
p_li_m_q_t_0.2 = which(li_m_q_t_0.2$p_value %in% targetPathway)
# p MAS quantile log.median
p_li_m_q_lm_0.05 = which(li_m_q_lm_0.05$p_value %in% targetPathway)
p_li_m_q_lm_0.1  = which(li_m_q_lm_0.1$p_value %in% targetPathway)
p_li_m_q_lm_0.15 = which(li_m_q_lm_0.15$p_value %in% targetPathway)
p_li_m_q_lm_0.2  = which(li_m_q_lm_0.2$p_value %in% targetPathway)
# p MAS quantile rlm
p_li_m_q_rlm_0.05 = which(li_m_q_rlm_0.05$p_value %in% targetPathway)
p_li_m_q_rlm_0.1  = which(li_m_q_rlm_0.1$p_value %in% targetPathway)
p_li_m_q_rlm_0.15 = which(li_m_q_rlm_0.15$p_value %in% targetPathway)
p_li_m_q_rlm_0.2  = which(li_m_q_rlm_0.2$p_value %in% targetPathway)
# p MAS scaling median.polish
p_li_m_s_mp_0.05 = which(li_m_s_mp_0.05$p_value %in% targetPathway)
p_li_m_s_mp_0.1  = which(li_m_s_mp_0.1$p_value %in% targetPathway)
p_li_m_s_mp_0.15 = which(li_m_s_mp_0.15$p_value %in% targetPathway)
p_li_m_s_mp_0.2  = which(li_m_s_mp_0.2$p_value %in% targetPathway)
# p MAS scaling tukey.biweight
p_li_m_s_t_0.05 = which(li_m_s_t_0.05$p_value %in% targetPathway)
p_li_m_s_t_0.1 = which(li_m_s_t_0.1$p_value %in% targetPathway)
p_li_m_s_t_0.15 = which(li_m_s_t_0.15$p_value %in% targetPathway)
p_li_m_s_t_0.2 = which(li_m_s_t_0.2$p_value %in% targetPathway)
# p MAS scaling log.median
p_li_m_s_lm_0.05 = which(li_m_s_lm_0.05$p_value %in% targetPathway)
p_li_m_s_lm_0.1 = which(li_m_s_lm_0.1$p_value %in% targetPathway)
p_li_m_s_lm_0.15 = which(li_m_s_lm_0.15$p_value %in% targetPathway)
p_li_m_s_lm_0.2 = which(li_m_s_lm_0.2$p_value %in% targetPathway)
# p MAS scaling rlm
p_li_m_s_rlm_0.05 = which(li_m_s_rlm_0.05$p_value %in% targetPathway)
p_li_m_s_rlm_0.05 = which(li_m_s_rlm_0.05$p_value %in% targetPathway)
p_li_m_s_rlm_0.05 = which(li_m_s_rlm_0.05$p_value %in% targetPathway)
p_li_m_s_rlm_0.05 = which(li_m_s_rlm_0.05$p_value %in% targetPathway)
# p LESN2 quantile median.polish ------------------------------------
p_li_l_q_mp_0.05 = which(li_l_q_mp_0.05$p_value %in% targetPathway)
p_li_l_q_mp_0.1 = which(li_l_q_mp_0.1$p_value %in% targetPathway)
p_li_l_q_mp_0.15 = which(li_l_q_mp_0.15$p_value %in% targetPathway)
p_li_l_q_mp_0.2 = which(li_l_q_mp_0.2$p_value %in% targetPathway)
# p LESN2 quantile tukey.biweight
p_li_l_q_t_0.05 = which(li_l_q_t_0.05$p_value %in% targetPathway)
p_li_l_q_t_0.1 = which(li_l_q_t_0.1$p_value %in% targetPathway)
p_li_l_q_t_0.15 = which(li_l_q_t_0.15$p_value %in% targetPathway)
p_li_l_q_t_0.2 = which(li_l_q_t_0.2$p_value %in% targetPathway)
# p LESN2 quantile log.median
p_li_l_q_lm_0.05 = which(li_l_q_lm_0.05$p_value %in% targetPathway)
p_li_l_q_lm_0.1  = which(li_l_q_lm_0.1$p_value %in% targetPathway)
p_li_l_q_lm_0.15 = which(li_l_q_lm_0.15$p_value %in% targetPathway)
p_li_l_q_lm_0.2  = which(li_l_q_lm_0.2$p_value %in% targetPathway)
# p LESN2 quantile rlm
p_li_l_q_rlm_0.05 = which(li_l_q_rlm_0.05$p_value %in% targetPathway)
p_li_l_q_rlm_0.1  = which(li_l_q_rlm_0.1$p_value %in% targetPathway)
p_li_l_q_rlm_0.15 = which(li_l_q_rlm_0.15$p_value %in% targetPathway)
p_li_l_q_rlm_0.2  = which(li_l_q_rlm_0.2$p_value %in% targetPathway)
# p LESN2 scaling median.polish
p_li_l_s_mp_0.05 = which(li_l_s_mp_0.05$p_value %in% targetPathway)
p_li_l_s_mp_0.1  = which(li_l_s_mp_0.1$p_value %in% targetPathway)
p_li_l_s_mp_0.15 = which(li_l_s_mp_0.15$p_value %in% targetPathway)
p_li_l_s_mp_0.2  = which(li_l_s_mp_0.2$p_value %in% targetPathway)
# p LESN2 scaling tukey.biweight
p_li_m_s_t_0.05 = which(li_m_s_t_0.05$p_value %in% targetPathway)
p_li_m_s_t_0.1 = which(li_m_s_t_0.1$p_value %in% targetPathway)
p_li_m_s_t_0.15 = which(li_m_s_t_0.15$p_value %in% targetPathway)
p_li_m_s_t_0.2 = which(li_m_s_t_0.2$p_value %in% targetPathway)
# p LESN2 scaling log.median
p_li_l_s_lm_0.05 = which(li_l_s_lm_0.05$p_value %in% targetPathway)
p_li_l_s_lm_0.1 = which(li_l_s_lm_0.1$p_value %in% targetPathway)
p_li_l_s_lm_0.15 = which(li_l_s_lm_0.15$p_value %in% targetPathway)
p_li_l_s_lm_0.2 = which(li_l_s_lm_0.2$p_value %in% targetPathway)
# p LESN2 scaling rlm
p_li_l_s_rlm_0.05 = which(li_l_s_rlm_0.05$p_value %in% targetPathway)
p_li_l_s_rlm_0.05 = which(li_l_s_rlm_0.05$p_value %in% targetPathway)
p_li_l_s_rlm_0.05 = which(li_l_s_rlm_0.05$p_value %in% targetPathway)
p_li_l_s_rlm_0.05 = which(li_l_s_rlm_0.05$p_value %in% targetPathway)

#### extract q-values of target gene sets ---------------------------------------------------------------
# q GCRMA quantile median.polish ----------------------------------
q_li_g_q_mp_0.05 = which(li_g_q_mp_0.05$p_value %in% targetPathway)
q_li_g_q_mp_0.1 = which(li_g_q_mp_0.1$p_value %in% targetPathway)
q_li_g_q_mp_0.15 = which(li_g_q_mp_0.15$p_value %in% targetPathway)
q_li_g_q_mp_0.2 = which(li_g_q_mp_0.2$p_value %in% targetPathway)
# q GCRMA quantile tukey.biweight
q_li_g_q_t_0.05 = which(li_g_q_t_0.05$p_value %in% targetPathway)
q_li_g_q_t_0.1 = which(li_g_q_t_0.1$p_value %in% targetPathway)
q_li_g_q_t_0.15 = which(li_g_q_t_0.15$p_value %in% targetPathway)
q_li_g_q_t_0.2 = which(li_g_q_t_0.2$p_value %in% targetPathway)
# q GCRMA quantile log.median
q_li_g_q_lm_0.05 = which(li_g_q_lm_0.05$p_value %in% targetPathway)
q_li_g_q_lm_0.1  = which(li_g_q_lm_0.1$p_value %in% targetPathway)
q_li_g_q_lm_0.15 = which(li_g_q_lm_0.15$p_value %in% targetPathway)
q_li_g_q_lm_0.2  = which(li_g_q_lm_0.2$p_value %in% targetPathway)
# q GCRMA quantile rlm
q_li_g_q_rlm_0.05 = which(li_g_q_rlm_0.05$p_value %in% targetPathway)
q_li_g_q_rlm_0.1 = which(li_g_q_rlm_0.1$p_value %in% targetPathway)
q_li_g_q_rlm_0.15 = which(li_g_q_rlm_0.15$p_value %in% targetPathway)
q_li_g_q_rlm_0.2 = which(li_g_q_rlm_0.2$p_value %in% targetPathway)
# q GCRMA scaling median.polish
q_li_g_s_mp_0.05 = which(li_g_s_mp_0.05$p_value %in% targetPathway)
q_li_g_s_mp_0.1 = which(li_g_s_mp_0.1$p_value %in% targetPathway)
q_li_g_s_mp_0.15 = which(li_g_s_mp_0.15$p_value %in% targetPathway)
q_li_g_s_mp_0.2 = which(li_g_s_mp_0.2$p_value %in% targetPathway)
# q GCRMA scaling tukey.biweight
q_li_g_s_t_0.05 = which(li_g_s_t_0.05$p_value %in% targetPathway)
q_li_g_s_t_0.1 = which(li_g_s_t_0.1$p_value %in% targetPathway)
q_li_g_s_t_0.15 = which(li_g_s_t_0.15$p_value %in% targetPathway)
q_li_g_s_t_0.2 = which(li_g_s_t_0.2$p_value %in% targetPathway)
# q GCRMA scaling log.median
q_li_g_s_lm_0.05 = which(li_g_s_lm_0.05$p_value %in% targetPathway)
q_li_g_s_lm_0.1 = which(li_g_s_lm_0.1$p_value %in% targetPathway)
q_li_g_s_lm_0.15 = which(li_g_s_lm_0.15$p_value %in% targetPathway)
q_li_g_s_lm_0.2 = which(li_g_s_lm_0.2$p_value %in% targetPathway)
# q GCRMA scaling rlm
q_li_g_s_rlm_0.05 = which(li_g_s_rlm_0.05$p_value %in% targetPathway)
q_li_g_s_rlm_0.05 = which(li_g_s_rlm_0.05$p_value %in% targetPathway)
q_li_g_s_rlm_0.05 = which(li_g_s_rlm_0.05$p_value %in% targetPathway)
q_li_g_s_rlm_0.05 = which(li_g_s_rlm_0.05$p_value %in% targetPathway)
# q MAS quantile median.polish --------------------------------------
q_li_m_q_mp_0.05 = which(li_m_q_mp_0.05$p_value %in% targetPathway)
q_li_m_q_mp_0.1 = which(li_m_q_mp_0.1$p_value %in% targetPathway)
q_li_m_q_mp_0.15 = which(li_m_q_mp_0.15$p_value %in% targetPathway)
q_li_m_q_mp_0.2 = which(li_m_q_mp_0.2$p_value %in% targetPathway)
# q MAS quantile tukey.biweight
q_li_m_q_t_0.05 = which(li_m_q_t_0.05$p_value %in% targetPathway)
q_li_m_q_t_0.1 = which(li_m_q_t_0.1$p_value %in% targetPathway)
q_li_m_q_t_0.15 = which(li_m_q_t_0.15$p_value %in% targetPathway)
q_li_m_q_t_0.2 = which(li_m_q_t_0.2$p_value %in% targetPathway)
# q MAS quantile log.median
q_li_m_q_lm_0.05 = which(li_m_q_lm_0.05$p_value %in% targetPathway)
q_li_m_q_lm_0.1  = which(li_m_q_lm_0.1$p_value %in% targetPathway)
q_li_m_q_lm_0.15 = which(li_m_q_lm_0.15$p_value %in% targetPathway)
q_li_m_q_lm_0.2  = which(li_m_q_lm_0.2$p_value %in% targetPathway)
# q MAS quantile rlm
q_li_m_q_rlm_0.05 = which(li_m_q_rlm_0.05$p_value %in% targetPathway)
q_li_m_q_rlm_0.1  = which(li_m_q_rlm_0.1$p_value %in% targetPathway)
q_li_m_q_rlm_0.15 = which(li_m_q_rlm_0.15$p_value %in% targetPathway)
q_li_m_q_rlm_0.2  = which(li_m_q_rlm_0.2$p_value %in% targetPathway)
# q MAS scaling median.polish
q_li_m_s_mp_0.05 = which(li_m_s_mp_0.05$p_value %in% targetPathway)
q_li_m_s_mp_0.1  = which(li_m_s_mp_0.1$p_value %in% targetPathway)
q_li_m_s_mp_0.15 = which(li_m_s_mp_0.15$p_value %in% targetPathway)
q_li_m_s_mp_0.2  = which(li_m_s_mp_0.2$p_value %in% targetPathway)
# q MAS scaling tukey.biweight
q_li_m_s_t_0.05 = which(li_m_s_t_0.05$p_value %in% targetPathway)
q_li_m_s_t_0.1 = which(li_m_s_t_0.1$p_value %in% targetPathway)
q_li_m_s_t_0.15 = which(li_m_s_t_0.15$p_value %in% targetPathway)
q_li_m_s_t_0.2 = which(li_m_s_t_0.2$p_value %in% targetPathway)
# q MAS scaling log.median
q_li_m_s_lm_0.05 = which(li_m_s_lm_0.05$p_value %in% targetPathway)
q_li_m_s_lm_0.1 = which(li_m_s_lm_0.1$p_value %in% targetPathway)
q_li_m_s_lm_0.15 = which(li_m_s_lm_0.15$p_value %in% targetPathway)
q_li_m_s_lm_0.2 = which(li_m_s_lm_0.2$p_value %in% targetPathway)
# q MAS scaling rlm
q_li_m_s_rlm_0.05 = which(li_m_s_rlm_0.05$p_value %in% targetPathway)
q_li_m_s_rlm_0.05 = which(li_m_s_rlm_0.05$p_value %in% targetPathway)
q_li_m_s_rlm_0.05 = which(li_m_s_rlm_0.05$p_value %in% targetPathway)
q_li_m_s_rlm_0.05 = which(li_m_s_rlm_0.05$p_value %in% targetPathway)
# q LESN2 quantile median.polish ------------------------------------
q_li_l_q_mp_0.05 = which(li_l_q_mp_0.05$p_value %in% targetPathway)
q_li_l_q_mp_0.1 = which(li_l_q_mp_0.1$p_value %in% targetPathway)
q_li_l_q_mp_0.15 = which(li_l_q_mp_0.15$p_value %in% targetPathway)
q_li_l_q_mp_0.2 = which(li_l_q_mp_0.2$p_value %in% targetPathway)
# q LESN2 quantile tukey.biweight
q_li_l_q_t_0.05 = which(li_l_q_t_0.05$p_value %in% targetPathway)
q_li_l_q_t_0.1 = which(li_l_q_t_0.1$p_value %in% targetPathway)
q_li_l_q_t_0.15 = which(li_l_q_t_0.15$p_value %in% targetPathway)
q_li_l_q_t_0.2 = which(li_l_q_t_0.2$p_value %in% targetPathway)
# q LESN2 quantile log.median
q_li_l_q_lm_0.05 = which(li_l_q_lm_0.05$p_value %in% targetPathway)
q_li_l_q_lm_0.1  = which(li_l_q_lm_0.1$p_value %in% targetPathway)
q_li_l_q_lm_0.15 = which(li_l_q_lm_0.15$p_value %in% targetPathway)
q_li_l_q_lm_0.2  = which(li_l_q_lm_0.2$p_value %in% targetPathway)
# q LESN2 quantile rlm
q_li_l_q_rlm_0.05 = which(li_l_q_rlm_0.05$p_value %in% targetPathway)
q_li_l_q_rlm_0.1  = which(li_l_q_rlm_0.1$p_value %in% targetPathway)
q_li_l_q_rlm_0.15 = which(li_l_q_rlm_0.15$p_value %in% targetPathway)
q_li_l_q_rlm_0.2  = which(li_l_q_rlm_0.2$p_value %in% targetPathway)
# q LESN2 scaling median.polish
q_li_l_s_mp_0.05 = which(li_l_s_mp_0.05$p_value %in% targetPathway)
q_li_l_s_mp_0.1  = which(li_l_s_mp_0.1$p_value %in% targetPathway)
q_li_l_s_mp_0.15 = which(li_l_s_mp_0.15$p_value %in% targetPathway)
q_li_l_s_mp_0.2  = which(li_l_s_mp_0.2$p_value %in% targetPathway)
# q LESN2 scaling tukey.biweight
q_li_m_s_t_0.05 = which(li_m_s_t_0.05$p_value %in% targetPathway)
q_li_m_s_t_0.1 = which(li_m_s_t_0.1$p_value %in% targetPathway)
q_li_m_s_t_0.15 = which(li_m_s_t_0.15$p_value %in% targetPathway)
q_li_m_s_t_0.2 = which(li_m_s_t_0.2$p_value %in% targetPathway)
# q LESN2 scaling log.median
q_li_l_s_lm_0.05 = which(li_l_s_lm_0.05$p_value %in% targetPathway)
q_li_l_s_lm_0.1 = which(li_l_s_lm_0.1$p_value %in% targetPathway)
q_li_l_s_lm_0.15 = which(li_l_s_lm_0.15$p_value %in% targetPathway)
q_li_l_s_lm_0.2 = which(li_l_s_lm_0.2$p_value %in% targetPathway)
# q LESN2 scaling rlm
q_li_l_s_rlm_0.05 = which(li_l_s_rlm_0.05$p_value %in% targetPathway)
q_li_l_s_rlm_0.05 = which(li_l_s_rlm_0.05$p_value %in% targetPathway)
q_li_l_s_rlm_0.05 = which(li_l_s_rlm_0.05$p_value %in% targetPathway)
q_li_l_s_rlm_0.05 = which(li_l_s_rlm_0.05$p_value %in% targetPathway)

}