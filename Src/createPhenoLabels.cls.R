fileConn = file(paste(celFilesPath, "phenotypes_GSEA.cls", sep = "/"))
phenoLabels = rep(0, length(expression_out) - 2)
phenoLabels[index_case] = set_case
phenoLabels[index_ctrl] = set_ctrl
if (phenoLabels[1] == set_case){
  pheno_order = paste("#", set_case, set_ctrl, sep = " ")
} else {
  pheno_order = paste( "#", set_ctrl, set_case, sep = " ")
}
phenoLabels = paste(phenoLabels, collapse = " ")
writeLines(c(paste(length(expression_out) - 2 , "2", "1", sep = " "), pheno_order, phenoLabels), fileConn)
close(fileConn)
#write.table(phenoLabels, file = paste( cel_files_path, "phenotypes_GSEA.cls", sep ="/" ), sep = " ", col.names = FALSE, row.names = F, append=TRUE)