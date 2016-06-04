p_value = function(pathwayGenes, diffGenes){
  pathwayGenes = unlist(pathwayGenes)
  quantity = sum(diffGenes %in% pathwayGenes)
  p_value = phyper(q = (quantity - 1), m = length(pathwayGenes), n = (length(counti) - length(pathwayGenes)), k = all_diffGenes, lower.tail = FALSE)
}

number = function(pathwayGenes, diffGenes){
  pathwayGenes = unlist(pathwayGenes)
  number = sum(diffGenes %in% pathwayGenes)
}

diffGenes = unique(as.character(topall_res$HGNC_symb))
diffGenes = diffGenes[! diffGenes %in% ""]
all_diffGenes = sum(diffGenes %in% counti)
all_non_diffGenes = length(counti) - all_diffGenes

p_values = unlist(lapply(gene_sets, FUN = p_value, diffGenes))
q_values = round(p.adjust(p_values, "BH"), 10)
diffGenes_set = unlist(lapply(gene_sets, FUN = number, diffGenes))

resultData = data.frame("GS" = pathwayID, "genes_in_set" = numberGenes, "diffGenes_in_set" = diffGenes_set, "percentage" = round((diffGenes_set / numberGenes) * 100, 3), "p_value" = p_values, "q_value" = q_values)
resultData = resultData[order(resultData$q_value, decreasing = FALSE), ]

GSOA_out = paste(outputPath, paste(backMethod, normalizeMethod, summaryMethod, sep = "_"), "limma_Results", paste(paste("pVal", str_replace(as.character(p_val), "\\.", "_"), sep = ""), paste("lFc", str_replace(as.character(lfc_exp), "\\.", "_"), sep = ""), sep = "_"), sep = "/")

write.table(resultData, file = paste(GSOA_out, "diffGenes_pathway_enrichment.csv", sep = "/"), row.names = FALSE, sep = ",")