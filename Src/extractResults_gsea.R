#### extract results gsea ---------------
# abbreviations
# g = GCRMA, m = MAS, l = LESN2, q = quantile, s = scaling, mp = median.polish, t = tukey.biweight, lm = log.median, rlm = rlm
# p = phenotype, gs = gene_set
# s2n = Signal2Noise, ts = tTest, lrc = log2_Ratio_of_Classes

targetPathway = as.character(target_pathways$target[target_pathways$DataSet %in% dataSet])

all_ranks_gs = c()
all_p_gs = c()
all_q_gs = c()
for(bm in 1:length(back_methods)){
  for(nm in 1:length(normalize_methods)){
    for(sm in 1:length(summary_methods)){
      for(ls in 1:length(local_statistics)){
        for(pt in 1:length(permutation_type)){
          current_path    = paste(outputPath, paste(back_methods[bm], normalize_methods[nm], summary_methods[sm], sep = "_"), paste("gseaResults", local_statistics[ls], permutation_type[pt], sep = "_"), sep = "/")
          results_path_gs = list.files(current_path, pattern = c(back_methods[bm], normalize_methods[nm], summary_methods[sm], local_statistics[ls], permutation_type[pt]))
          results_path_gs = paste(current_path, results_path_gs, sep = "/")
          file_disease    = list.files(results_path_gs, pattern = c("gsea_report_for_Disease_.*\\.xls$"))
          path_disease    = paste(results_path_gs, file_disease, sep = "/")
          disease_table   = read.table(path_disease, header = TRUE, sep = "\t")
          file_normal     = list.files(results_path_gs, pattern = c("gsea_report_for_Normal_.*\\.xls$"))
          path_normal     = paste(results_path_gs, file_normal, sep = "/")
          normal_table    = read.table(path_normal, header = TRUE, sep = "\t")
          normal_disease  = rbind(disease_table, normal_table)
          normal_disease  = normal_disease[order(normal_disease$FDR.q.val, normal_disease$NOM.p.val, decreasing = FALSE), ]
          current_rank    = which(normal_disease$NAME %in% toupper(targetPathway))
          all_ranks_gs    = c(all_ranks_gs, current_rank)
          all_p_gs        = c(all_p_gs, normal_disease$NOM.p.val[current_rank])
          all_q_gs        = c(all_q_gs, normal_disease$FDR.q.val[current_rank])
        }
      }
    }
  }
}

parameterSet_gs = c("m_q_mp_s2n_p", "m_q_mp_s2n_gs", "m_q_mp_ts_p", "m_q_mp_ts_gs", "m_q_mp_lrc_p", "m_q_mp_lrc_gs", 
                    "m_q_t_s2n_p", "m_q_t_s2n_gs", "m_q_t_ts_p", "m_q_t_ts_gs", "m_q_t_lrc_p", "m_q_t_lrc_gs",
                    "m_q_lm_s2n_p", "m_q_lm_s2n_gs", "m_q_lm_ts_p", "m_q_lm_ts_gs", "m_q_lm_lrc_p", "m_q_lm_lrc_gs",
                    "m_q_rlm_s2n_p", "m_q_rlm_s2n_gs", "m_q_rlm_ts_p", "m_q_rlm_ts_gs", "m_q_rlm_lrc_p", "m_q_rlm_lrc_gs",  
                    "m_s_mp_s2n_p", "m_s_mp_s2n_gs", "m_s_mp_ts_p", "m_s_mp_ts_gs", "m_s_mp_lrc_p", "m_s_mp_lrc_gs", 
                    "m_s_t_s2n_p", "m_s_t_s2n_gs", "m_s_t_ts_p", "m_s_t_ts_gs", "m_s_t_lrc_p", "m_s_t_lrc_gs",
                    "m_s_lm_s2n_p", "m_s_lm_s2n_gs", "m_s_lm_ts_p", "m_s_lm_ts_gs", "m_s_lm_lrc_p", "m_s_lm_lrc_gs",
                    "m_s_rlm_s2n_p", "m_s_rlm_s2n_gs", "m_s_rlm_ts_p", "m_s_rlm_ts_gs", "m_s_rlm_lrc_p", "m_s_rlm_lrc_gs",
                    "g_q_mp_s2n_p", "g_q_mp_s2n_gs", "g_q_mp_ts_p", "g_q_mp_ts_gs", "g_q_mp_lrc_p", "g_q_mp_lrc_gs", 
                    "g_q_t_s2n_p", "g_q_t_s2n_gs", "g_q_t_ts_p", "g_q_t_ts_gs", "g_q_t_lrc_p", "g_q_t_lrc_gs",
                    "g_q_lm_s2n_p", "g_q_lm_s2n_gs", "g_q_lm_ts_p", "g_q_lm_ts_gs", "g_q_lm_lrc_p", "g_q_lm_lrc_gs",
                    "g_q_rlm_s2n_p", "g_q_rlm_s2n_gs", "g_q_rlm_ts_p", "g_q_rlm_ts_gs", "g_q_rlm_lrc_p", "g_q_rlm_lrc_gs",  
                    "g_s_mp_s2n_p", "g_s_mp_s2n_gs", "g_s_mp_ts_p", "g_s_mp_ts_gs", "g_s_mp_lrc_p", "g_s_mp_lrc_gs", 
                    "g_s_t_s2n_p", "g_s_t_s2n_gs", "g_s_t_ts_p", "g_s_t_ts_gs", "g_s_t_lrc_p", "g_s_t_lrc_gs",
                    "g_s_lm_s2n_p", "g_s_lm_s2n_gs", "g_s_lm_ts_p", "g_s_lm_ts_gs", "g_s_lm_lrc_p", "g_s_lm_lrc_gs",
                    "g_s_rlm_s2n_p", "g_s_rlm_s2n_gs", "g_s_rlm_ts_p", "g_s_rlm_ts_gs", "g_s_rlm_lrc_p", "g_s_rlm_lrc_gs",  
                    "l_q_mp_s2n_p", "l_q_mp_s2n_gs", "l_q_mp_ts_p", "l_q_mp_ts_gs", "l_q_mp_lrc_p", "l_q_mp_lrc_gs", 
                    "l_q_t_s2n_p", "l_q_t_s2n_gs", "l_q_t_ts_p", "l_q_t_ts_gs", "l_q_t_lrc_p", "l_q_t_lrc_gs",
                    "l_q_lm_s2n_p", "l_q_lm_s2n_gs", "l_q_lm_ts_p", "l_q_lm_ts_gs", "l_q_lm_lrc_p", "l_q_lm_lrc_gs",
                    "l_q_rlm_s2n_p", "l_q_rlm_s2n_gs", "l_q_rlm_ts_p", "l_q_rlm_ts_gs", "l_q_rlm_lrc_p", "l_q_rlm_lrc_gs",  
                    "l_s_mp_s2n_p", "l_s_mp_s2n_gs", "l_s_mp_ts_p", "l_s_mp_ts_gs", "l_s_mp_lrc_p", "l_s_mp_lrc_gs", 
                    "l_s_t_s2n_p", "l_s_t_s2n_gs", "l_s_t_ts_p", "l_s_t_ts_gs", "l_s_t_lrc_p", "l_s_t_lrc_gs",
                    "l_s_lm_s2n_p", "l_s_lm_s2n_gs", "l_s_lm_ts_p", "l_s_lm_ts_gs", "l_s_lm_lrc_p", "l_s_lm_lrc_gs",
                    "l_s_rlm_s2n_p", "l_s_rlm_s2n_gs", "l_s_rlm_ts_p", "l_s_rlm_ts_gs", "l_s_rlm_lrc_p", "l_s_rlm_lrc_gs")

result_out_gs = data.frame("parameter_setting" = parameterSet_gs, "rank" = all_ranks_gs, "p_value" = all_p_gs, "q_value" = all_q_gs)
result_out_gs = result_out_gs[order(result_out_gs$rank, result_out_gs$q_value, decreasing = FALSE), ]
resultFile_gs = paste(outputPath, paste("overall_results_", dataSet, "_gsea", ".txt", sep = ""), sep = "/")
write.table(result_out_gs, file = resultFile_gs, sep = "\t", quote = FALSE, row.names = FALSE)

# write best/worst ranks, p.values and p-values for limma to file
best_rank_gs = min(all_ranks_gs)
best_rank_gs_setting = which(all_ranks_gs %in% best_rank_gs)
best_rank_gs_setting = parameterSet_gs[best_rank_gs_setting]
worst_rank_gs = max(all_ranks_gs)
worst_rank_gs_setting = which(all_ranks_gs %in% worst_rank_gs)
worst_rank_gs_setting = parameterSet_gs[worst_rank_gs_setting]

best_p_gs = min(all_p_gs)
best_p_gs_setting = which(all_p_gs %in% best_p_gs)
best_p_gs_setting = parameterSet_gs[best_p_gs_setting]
worst_p_gs = max(all_p_gs)
worst_p_gs_setting = which(all_p_gs %in% worst_p_gs)
worst_p_gs_setting = parameterSet_gs[worst_p_gs_setting]

best_q_gs = min(all_q_gs)
best_q_gs_setting = which(all_q_gs %in% best_q_gs)
best_q_gs_setting = parameterSet_gs[best_q_gs_setting]
worst_q_gs = max(all_q_gs)
worst_q_gs_setting = which(all_q_gs %in% worst_q_gs)
worst_q_gs_setting = parameterSet_gs[worst_q_gs_setting]

resultFile_best_gs = paste(outputPath, paste("best_worst_results_", dataSet, "_gsea", ".txt", sep = ""), sep = "/")
write(paste(paste("best rank: ", best_rank_gs, sep = "  "), paste("parameter setting: ", best_rank_gs_setting, sep = "  "), sep = "  "), file = resultFile_best_gs, append = TRUE)
write(paste(paste("worst rank: ", worst_rank_gs, sep = "  "), paste("parameter setting: ", worst_rank_gs_setting, sep = "  "), sep = "  "), file = resultFile_best_gs, append = TRUE)
write(paste(paste("best p-value: ", best_p_gs, sep = "  "), paste("parameter setting: ", best_p_gs_setting, sep = "  "), sep = "  "), file = resultFile_best_gs, append = TRUE)
write(paste(paste("worst p-value: ", worst_p_gs, sep = "  "), paste("parameter setting: ", worst_p_gs_setting, sep = "  "), sep = "  "), file = resultFile_best_gs, append = TRUE)
write(paste(paste("best q-value: ", best_q_gs, sep = "  "), paste("parameter setting: ", best_q_gs_setting, sep = "  "), sep = "  "), file = resultFile_best_gs, append = TRUE)
write(paste(paste("worst q-value: ", worst_q_gs, sep = "  "), paste("parameter setting: ", worst_q_gs_setting, sep = "  "), sep = "  "), file = resultFile_best_gs, append = TRUE)