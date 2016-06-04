kocent = (system('uname -n', intern = TRUE) == "kocent")
if (kocent){
  pipelineLoc = "/usr/Generic_mRNA_Expression_Pipeline" #server path
}else if(system('uname -n', intern = T) == "MacBook-Jan-Niklas.local" ){
  pipelineLoc = paste(system("echo $HOME", intern = T), "OptimizeParameters", sep = "/")
}

pipelineLoc = "/Users/jan-niklas/OptimizeParameters"
setwd(pipelineLoc)
source("Src/loadPackages.R")

set_ctrl = "Normal"
set_case = "Disease"
p_val = 0.05

data_sets = c("GSE1297", "GSE14762", "GSE14924_CD4", "GSE14924_CD8", "GSE15932_Dia", "GSE15932_Panc", "GSE19420", "GSE19728",
              "GSE20153", "GSE20164", "GSE20291", "GSE21354", "GSE24250", "GSE30153", "GSE32676", "GSE3585", "GSE4107", "GSE4183",
"GSE5281_EC", "GSE5281_HIP", "GSE5281_VCX", "GSE781", "GSE8762", "GSE9348", "GSE9476")
data_sets = c("GSE1297")#, "GSE14762")
pack_foreach = c("stringr", "hgu133plus2.db", "affy", "simpleaffy", "affyPLM", "affycoretools", "affyQCReport", "annaffy", "limma", "xlsx")

source("Src/parametersToOptimize.R")

numCores = detectCores() - 1
registerDoMC(numCores)
#clusters = makeCluster(numCores)
#clusterExport(clusters, c("data_sets", "pipelineLoc", "set_ctrl", "set_case", "p_val", "lfc_exp"))
#registerDoParallel(clusters)

source("Src/processGeneSetDB.R")

target_pathways = read.table("~/OptimizeParameters/Misc/target_pathways.tab", sep = "\t", header = TRUE)

foreach(i = 1:length(data_sets), .packages = pack_foreach) %dopar% {

  dataSet = data_sets[i]
  #### set generic initial parameters ####
  source("Src/setGenericParameters.R", local = TRUE)
  
  #### read cel files from input directory ####
  source("Src/readCEL.R", local = TRUE)
  
  #### create cohorts for data ####
  source("Src/createCohorts.R", local = TRUE)
  
  #### quality control ####
  if(! file.exists(paste(outputPath, "QC_report.pdf"))){
    QCReport(rawData, file = paste(qcPath, "QC_report.pdf", sep = "/"))
  }
  
  #### perform normalization with changing parameters ####
  foreach(b = 1:length(back_methods)) %dopar% {
    backMethod = back_methods[b]
    foreach(n = 1:length(normalize_methods)) %dopar% {
      normalizeMethod = normalize_methods[n]
      foreach(s = 1:length(summary_methods)) %dopar% {
        summaryMethod = summary_methods[s]
      
        source("Src/normalize.r", local = TRUE)
  
        #### annotate data ####
        source("Src/annotation.R", local = TRUE)
    
        for(p in 1:length(p_change)){
          p_val = p_change[p]
          for(l in 1:length(lfc_change)){
            lfc_exp = lfc_change[l]
            #### limma analysis ####
            source("Src/computeLimma.R", local = TRUE)
  
            #### limma gene set analysis ####
            source("Src/limmaGSOA.R", local = TRUE)
          }
        }
        
        #### create expressionSet.gct if doesn't exist ####
        if(! file.exists(paste(celFilesPath, paste("ExpressionSet_", paste(backMethod, normalizeMethod, summaryMethod, sep = "_"), ".gct", sep =""), sep = "/"))){
          source("Src/createExpression.gct.R", local = TRUE)
        }
        #### create phenotype label cls file if doesn't exist ####
        if(! file.exists(paste(celFilesPath, "phenotypes_GSEA.cls", sep = "/"))){
          source("Src/createPhenoLabels.cls.R", local = TRUE)  
        }
    
        for(g in 1:length(permutation_type)) %dopar% {
          permutation = permutation_type[g]
          for(h in 1:length(local_statistics)) %dopar% {
            local_statistic = local_statistics[h]
            #### gsea analysis ####
            source("Src/computeGSEA.R", local = TRUE)
          }
        }
      }
    }
  }
  #### performance limma - extract ranks, p-values, q-values
  source("Src/extractResults_limma.R", local = TRUE)
  
  #### performance gsea - extract ranks, p-values, q-values
  source("Src/extractResults_gsea.R", local = TRUE)
  
  #### write final results
  source("Src/writeFinalResults.R", local = TRUE)
}
stopCluster(clusters)