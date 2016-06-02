probe_ids      = unlist(hgnc_symbols)
descr          = rep("NA", length(probe_ids))
data           = as.data.frame(exprs(eset))
expression_out = data.frame("NAMES" = probe_ids, "DESCRIPTION" = descr, data)

index_case_gsea = sapply(index_case, function(x) x + 2)
index_ctrl_gsea = sapply(index_ctrl, function(x) x + 2)
expression_out = expression_out[expression_out$NAMES != "", ]
names_dupli = expression_out$NAMES[which(duplicated(expression_out$NAMES))]
names_dupli = unique(names_dupli)
reduced = list()
  
if(length(names_dupli) > 0){
    
  # collapsing duplicate gene symbols by logFC
  index_dupli = which(expression_out$NAMES %in% names_dupli)
  for (i in 1:length(names_dupli)){
    index = which(expression_out$NAMES == names_dupli[i])
    current = expression_out[index, ]
    exprs_case_gsea = rowMeans(current[, index_case_gsea ])
    exprs_ctrl_gsea = rowMeans(current[, index_ctrl_gsea ])
    dif_exp_gsea    = exprs_case_gsea - exprs_ctrl_gsea
    index_max = which.max(abs(dif_exp_gsea))
    temp      = as.data.frame(current[index_max,])
    reduced   = c(reduced, list(temp))
  }
  reduced = do.call("rbind", reduced)
  expression_out = expression_out[- index_dupli, ]
  expression_out = rbind(expression_out, reduced)
}
  
# write ExpressionSet.gct file in Input directory
fileConn = file(paste(celFilesPath, paste("ExpressionSet_", backMethod, ".gct", sep = ""), sep = "/"))
writeLines(c("#1.2", paste(length(expression_out$NAMES), length(expression_out) - 2, sep = "\t")), fileConn)
close(fileConn)
write.table(expression_out, file = paste(celFilesPath, paste("ExpressionSet_", backMethod, ".gct", sep = ""), sep ="/" ), sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)