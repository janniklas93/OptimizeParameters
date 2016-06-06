probe_ids      = unlist(hgnc_symbols)
descr          = rep("NA", length(probe_ids))
data           = as.data.frame(exprs(eset))
expression_out = data.frame("NAMES" = probe_ids, "DESCRIPTION" = descr, data)

index_case_gsea = sapply(index_case, function(x) x + 2)
index_ctrl_gsea = sapply(index_ctrl, function(x) x + 2)
expression_out = expression_out[expression_out$NAMES != "", ]