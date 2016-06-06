geneSetDB_path = paste(pipelineLoc, "Misc", "kegg_gene_sets.gmt", sep = "/")
gene_set_DB = readLines(geneSetDB_path)
gene_set_DB = lapply(gene_set_DB, FUN = strsplit, "\t")

counti = c()
gene_sets = list()
pathwayID = c()
numberGenes = c()

#### process gene set db ####
for(i in 1:length(gene_set_DB)){
  temp        = unlist(gene_set_DB[[i]])
  temp_path   = temp[1]
  pathwayID   = c(pathwayID, temp_path)
  temp        = temp[- c(1, 2)]
  temp_number = length(temp)
  numberGenes = c(numberGenes, temp_number)
  counti      = c(counti, temp)
  singleSet   = as.list(temp)
  gene_sets[[length(gene_sets) + 1]] = singleSet
}