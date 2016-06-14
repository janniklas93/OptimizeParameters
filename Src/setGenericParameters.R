date = Sys.Date()
date = str_replace_all(date, "-", "_")
date = "2016_06_13"
project_name = paste("Analysis_", date, sep = "")
celFilesPath = paste(pipelineLoc, "DataSets", dataSet, "Input", sep = "/")
outputPath = paste(dirname(celFilesPath), "Output", project_name, sep = "/")
dir.create(outputPath)

#mas5_path  = paste(outputPath, "MAS", sep = "/"); dir.create(mas5_path)
#lesn2_path = paste(outputPath, "LESN2", sep = "/"); dir.create(lesn2_path)
#gcrma_path = paste(outputPath, "GCRMA", sep = "/"); dir.create(gcrma_path)

final_results_path = paste(pipelineLoc, "DataSets", dataSet, "Output/Final_Results" , paste("Analysis_", date, sep = ""), sep = "/")
dir.create(final_results_path, showWarnings = FALSE, recursive = TRUE)
best_ranks_file      = paste(final_results_path, "best_Ranks.tab", sep = "/")
best_q.values_file   = paste(final_results_path, "best_Q.Values.tab", sep = "/")
best_p.values_file   = paste(final_results_path, "best_P.Values.tab", sep = "/")
worst_ranks_file     = paste(final_results_path, "worst_Ranks.tab", sep = "/")
worst_q.values_file  = paste(final_results_path, "worst_Q.Values.tab", sep = "/")
worst_p.values_file  = paste(final_results_path, "worst_P.Values.tab", sep = "/")
# initialize final result file
write(paste("Dataset", "limma", "gsea", sep = "\t"), file = best_ranks_file)
write(paste("Dataset", "limma", "gsea", sep = "\t"), file = best_q.values_file)
write(paste("Dataset", "limma", "gsea", sep = "\t"), file = best_p.values_file)
write(paste("Dataset", "limma", "gsea", sep = "\t"), file = worst_ranks_file)
write(paste("Dataset", "limma", "gsea", sep = "\t"), file = worst_q.values_file)
write(paste("Dataset", "limma", "gsea", sep = "\t"), file = worst_p.values_file)
  
cohortsPath    = paste(celFilesPath, "cohorts.tab", sep = "/")
qcPath         = dirname(outputPath)
geneSetDB_path = paste(pipelineLoc, "Misc", "kegg_gene_sets.gmt", sep = "/")
user_folder    = as.character(system("echo $HOME", intern = TRUE))
#dir.create(limmaPath, showWarnings = FALSE)

phenodata = read.table(cohortsPath, header = TRUE, sep = "\t" )
  
strEndsWith <- function(haystack, needle)
{
  hl <- nchar(haystack)
  nl <- nchar(needle)
  if(nl>hl)
  {
    return(F)
  } else
  {
    return(substr(haystack, hl-nl+1, hl) == needle)
  }
}