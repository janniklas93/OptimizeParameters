date = Sys.Date()
project_name = paste("Analysis_", date, sep = "")
celFilesPath = paste(pipelineLoc, "DataSets", dataSet, "Input", sep = "/")
outputPath = paste(dirname(celFilesPath), "Output", project_name, sep = "/")
dir.create(outputPath)

mas5_path  = paste(outputPath, "MAS", sep = "/"); dir.create(mas5_path)
lesn2_path = paste(outputPath, "LESN2", sep = "/"); dir.create(lesn2_path)
gcrma_path = paste(outputPath, "GCRMA", sep = "/"); dir.create(gcrma_path)

final_results_path = paste(pipelineLoc, "DataSets/Final_Results" , paste("Analysis_", date, sep = ""), sep = "/")
dir.create(final_results_path, showWarnings = FALSE, recursive = TRUE)
all_ranks_file     = paste(final_results_path, "allRanks.tab", sep = "/")
all_q.values_file  = paste(final_results_path, "allQ.Values.tab", sep = "/")
all_p.values_file  = paste(final_results_path, "allP.Values.tab", sep = "/")
# initialize final result file
write(paste("Dataset", "limma", "gsea", sep = "\t"), file = all_ranks_file)
write(paste("Dataset", "limma", "gsea", sep = "\t"), file = all_q.values_file)
write(paste("Dataset", "limma", "gsea", sep = "\t"), file = all_p.values_file)
  
cohortsPath    = paste(celFilesPath, "cohorts.tab", sep = "/")
qcPath         = dirname(outputPath)
geneSetDB_path = paste(pipelineLoc, "Misc", "kegg_biocarta.gmt", sep = "/")
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