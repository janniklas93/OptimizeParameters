#topall = volc_all[-which(hgnc_symbols %in% ""), ]
#hgnc_sym_limma = hgnc_symbols[-which(hgnc_symbols %in% "")]
#topall = cbind(topall, hgnc_sym_limma)
#topall = topall[topall$adj.P.Val <= 0.05, ]
#topall = topall[order(abs(topall$logFC), decreasing = TRUE), ]
#topall = topall[1:(anteil_diff * length(topall[, 1])), ]

topall = topTable(fit, coef = "contrast", number  = nrow(eset), adjust = "none", p.value = p_val, lfc = lfc_exp)

if(dim(topall)[1] == 0){
  topall_res = data.frame(
    "Probe_ids"           = "no_diff_exp",
    "expr_ctrl"           = 999,
    "expr_case"           = 999,
    "logFC"               = 999,
    "P_Value"             = 999,
    "HGNC_symb"           = "NNNNNNN",
    "HGNC_name"           = "NNNNNNN"
  )
} else{
  probe_ids        = rownames(topall)
  index_topall     = match(probe_ids, rownames(eset), nomatch = 0)
  if ((dim(topall)[1] == 1)){
    exprs_case = mean(exprs(eset)[index_topall, index_case])
    exprs_ctrl = mean(exprs(eset)[index_topall, index_ctrl])
  }else{
    exprs_case = rowMeans(exprs(eset)[index_topall, index_case])
    exprs_ctrl = rowMeans(exprs(eset)[index_topall, index_ctrl])
  }
  topall$logFC     = round(topall$logFC, 2)
  hgnc_sym_limma   = hgnc_symbols[index_topall]
  hgnc_names_limma = hgnc_names[index_topall]
  ensembl_limma    = ensembl_genes[index_topall]

  topall_res = data.frame(
    "Probe_ids"           = probe_ids,
    "expr_ctrl"           = round(exprs_ctrl, 2),
    "expr_case"           = round(exprs_case, 2),
    "logFC"               = topall$logFC,
    "P_Value"             = topall$P.Val,
    "HGNC_symb"           = hgnc_sym_limma,
    "HGNC_name"           = stringr::str_replace_all(hgnc_names_limma, ",", ";")
  )
}

topall_res = topall_res[order(topall_res$logFC, decreasing = TRUE), ]
dir.create(paste(outputPath, paste(backMethod, normalizeMethod, summaryMethod, sep = "_"),  "limma_Results", paste(paste("pVal", str_replace(as.character(p_val), "\\.", "_"), sep = ""), paste("lFc", str_replace(as.character(lfc_exp), "\\.", "_"), sep = ""), sep = "_"), sep = "/"), showWarnings = FALSE, recursive = TRUE)
write.table(topall_res, file = paste(outputPath, paste(backMethod, normalizeMethod, summaryMethod, sep = "_"), "limma_Results", paste(paste("pVal", str_replace(as.character(p_val), "\\.", "_"), sep = ""), paste("lFc", str_replace(as.character(lfc_exp), "\\.", "_"), sep = ""), sep = "_"), "dif_exp_results.csv", sep = "/"), row.names = FALSE, sep = ",")