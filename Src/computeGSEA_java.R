gseaOutputPath = paste(outputPath, paste(backMethod, normalizeMethod, summaryMethod, sep = "_"), paste("gseaResults", local_statistic, permutation, sep = "_"), sep = "/")
dir.create(gseaOutputPath, recursive = TRUE, showWarnings = FALSE)
project_ID     = paste(dataSet, backMethod, normalizeMethod, summaryMethod, permutation, local_statistic, sep = "_")

# create parameter file
para_ID = c("res", "cls", "gmx", "out", "permute", "metric", "rpt_label", "collapse", "gui")
para_values = c(
  paste(celFilesPath, paste("ExpressionSet_", paste(backMethod, normalizeMethod, summaryMethod, sep = "_"), ".gct", sep = ""), sep = "/"),
  paste(celFilesPath, "phenotypes_GSEA.cls", sep = "/"),
  geneSetDB_path,
  gseaOutputPath,
  permutation,
  local_statistic,
  project_ID,
  "false",
  "false")

parameters = data.frame(para_ID, para_values)
para_path = paste(gseaOutputPath, "Run_Parameters.txt", sep = "/")
write.table(parameters, file = para_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

system(paste("sudo java -Xmx1024m -cp", paste(pipelineLoc, "gsea2-2.2.2.jar", sep = "/"), "xtools.gsea.Gsea -param_file", para_path, sep = " "), intern = TRUE)
#system("java -Xmx1024m -cp /Users/Shared/OptimizeParameters/gsea2-2.2.2.jar xtools.gsea.Gsea -param_file /Users/Shared/OptimizeParameters/DataSets/GSE1297/Output/Analysis_2016_06_03/GCRMA_scaling_median.polish/gseaResults_tTest_phenotype/Run_Parameters.txt")