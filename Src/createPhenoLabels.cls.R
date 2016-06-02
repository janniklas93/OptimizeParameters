fileConn = file(paste(celFilesPath, "phenotypes_GSEA.cls", sep = "/"))
phenoLabels = rep(0, length(expression_out) - 2)
case_GSEA = unlist(strsplit(set_case, " "))
case_GSEA = case_GSEA[1]
ctrl_GSEA = unlist(strsplit(set_ctrl, " "))
ctrl_GSEA = ctrl_GSEA[1]
phenoLabels[index_case] = case_GSEA
phenoLabels[index_ctrl] = ctrl_GSEA
phenoLabels = paste(phenoLabels, collapse = " ")
if (phenoLabels[1] == set_case){
  pheno_order = paste("#", case_GSEA, ctrl_GSEA, sep = " ")
} else {
  pheno_order = paste( "#", ctrl_GSEA, case_GSEA, sep = " ")
}
writeLines(c(paste(length(expression_out) - 2 , "2", "1", sep = " "), pheno_order, phenoLabels), fileConn)
close(fileConn)
#write.table(phenoLabels, file = paste( cel_files_path, "phenotypes_GSEA.cls", sep ="/" ), sep = " ", col.names = FALSE, row.names = F, append=TRUE)