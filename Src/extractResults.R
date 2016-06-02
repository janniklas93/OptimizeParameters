result_limma_GCRMA = read.table(paste(outputPath, "GCRMA/limma_Results/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
result_limma_MAS   = read.table(paste(outputPath, "MAS/limma_Results/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")
result_limma_LESN2 = read.table(paste(outputPath, "LESN2/limma_Results/diffGenes_pathway_enrichment.csv", sep = "/"), header = TRUE, sep = ",")

result_gsea_GCRMA_Normal  = read.table(paste(outputPath, "GCRMA/gseaResults", paste(dataSet, "_GCRMA.SUMMARY.RESULTS.REPORT.Normal.txt", sep = ""), sep = "/"), sep = "\t", header = TRUE)
result_gsea_GCRMA_Disease = read.table(paste(outputPath, "GCRMA/gseaResults", paste(dataSet, "_GCRMA.SUMMARY.RESULTS.REPORT.Disease.txt", sep = ""), sep = "/"), sep = "\t", header = TRUE)
result_gsea_MAS_Normal    = read.table(paste(outputPath, "MAS/gseaResults", paste(dataSet, "_MAS.SUMMARY.RESULTS.REPORT.Normal.txt", sep = ""), sep = "/"), sep = "\t", header = TRUE)
result_gsea_MAS_Disease   = read.table(paste(outputPath, "MAS/gseaResults", paste(dataSet, "_MAS.SUMMARY.RESULTS.REPORT.Disease.txt", sep = ""), sep = "/"), sep = "\t", header = TRUE)
result_gsea_LESN2_Normal  = read.table(paste(outputPath, "LESN2/gseaResults", paste(dataSet, "_LESN2.SUMMARY.RESULTS.REPORT.Normal.txt", sep = ""), sep = "/"), sep = "\t", header = TRUE)
result_gsea_LESN2_Disease = read.table(paste(outputPath, "LESN2/gseaResults", paste(dataSet, "_LESN2.SUMMARY.RESULTS.REPORT.Normal.txt", sep = ""), sep = "/"), sep = "\t", header = TRUE)

result_gsea_GCRMA = rbind(result_gsea_GCRMA_Normal, result_gsea_GCRMA_Disease)
result_gsea_GCRMA = result_gsea_GCRMA[order(result_gsea_GCRMA$FDR.q.val, result_gsea_GCRMA$NOM.p.val, decreasing = FALSE), ]
result_gsea_MAS = rbind(result_gsea_MAS_Normal, result_gsea_MAS_Disease)
result_gsea_MAS = result_gsea_MAS[order(result_gsea_MAS$FDR.q.val, result_gsea_MAS$NOM.p.val, decreasing = FALSE), ]
result_gsea_LESN2 = rbind(result_gsea_LESN2_Normal, result_gsea_LESN2_Disease)
result_gsea_LESN2 = result_gsea_GCRMA[order(result_gsea_GCRMA$FDR.q.val, result_gsea_LESN2$NOM.p.val, decreasing = FALSE), ]

rm(result_gsea_GCRMA_Normal, result_gsea_GCRMA_Disease, result_gsea_MAS_Normal, result_gsea_MAS_Disease)

targetPathway = as.character(target_pathways$target[target_pathways$DataSet %in% dataSet])

# extract rank of target geneSet 
rank_GCRMA_lm = which(result_limma_GCRMA$GS %in% targetPathway)
rank_MAS_lm   = which(result_limma_MAS$GS %in% targetPathway)
rank_LESN2_lm = which(result_limma_LESN2$GS %in% targetPathway)
rank_GCRMA_gs = which(result_gsea_GCRMA$GS %in% targetPathway)
rank_MAS_gs   = which(result_gsea_MAS$GS %in% targetPathway)
rank_LESN2_gs = which(result_gsea_LESN2$GS %in% targetPathway)

# extract p-values of target gene sets 
q_GCRMA_lm = result_limma_GCRMA$q_value[rank_GCRMA_lm]
q_MAS_lm   = result_limma_MAS$q_value[rank_MAS_lm]
q_LESN2_lm = result_limma_LESN2$q_value[rank_LESN2_lm]
q_GCRMA_gs = result_gsea_GCRMA$FDR.q.val[rank_GCRMA_gs]
q_MAS_gs   = result_gsea_MAS$FDR.q.val[rank_MAS_gs]
q_LESN2_gs = result_gsea_LESN2$FDR.q.val[rank_LESN2_gs]

# extract p-values of target gene sets 
p_GCRMA_lm = result_limma_GCRMA$p_value[rank_GCRMA_lm]
p_MAS_lm   = result_limma_MAS$p_value[rank_MAS_lm]
p_LESN2_lm = result_limma_LESN2$p_value[rank_LESN2_lm]
p_GCRMA_gs = result_gsea_GCRMA$NOM.p.val[rank_GCRMA_gs]
p_MAS_gs   = result_gsea_MAS$NOM.p.val[rank_MAS_gs]
p_LESN2_gs = result_gsea_LESN2$NOM.p.val[rank_LESN2_gs]

# prepare output to file
rank_all_lm = c(rank_GCRMA_lm, rank_MAS_lm, rank_LESN2_lm)
q_all_lm = c(q_GCRMA_lm, q_MAS_lm, q_LESN2_lm)
p_all_lm = c(p_GCRMA_lm, p_MAS_lm, p_LESN2_lm))

rank_all_gs = c(rank_GCRMA_gs, rank_MAS_gs, rank_LESN2_lm)
q_all_gs = c(q_GCRMA_gs, q_MAS_gs, q_LESN2_lm)
p_all_gs = c(p_GCRMA_gs, p_MAS_gs, p_LESN2_lm))

resultFile = paste(outputPath, paste("Optimization_Results_", dataSet, ".txt", sep = ""), sep = "/")

result_limma = data.frame(
  "Parameters" = back_methods,
  "Rank"       = rank_all_lm,
  "q_value"    = q_all_lm,
  "p_value"    = p_all_lm)

result_gsea = data.frame(
  "Parameters" = back_methods[c(1:2)],
  "Rank"       = rank_all_gs,
  "q_value"    = q_all_gs,
  "p_value"    = p_all_gs
)

# write result report
write(c(paste("Limma Results ", dataSet , sep = "")), file = resultFile , append = TRUE)
write.table(result_gsea, file = resultFile, sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE)
write(c("\n"), file = resultFile , append = TRUE)
write(c(paste("GSEA Results ", dataSet , sep = "")), file = resultFile , append = TRUE)
write.table(result_gsea, file = resultFile, sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE)

# write best ranks, p.values and p-values for limma and gsea to file
best_rank_lm = min(rank_all_lm)
best_rank_gs = min(rank_all_gs)
write(paste(dataSet, best_rank_lm, best_rank_gs, sep = "\t"), file = all_ranks_file, append = TRUE)
