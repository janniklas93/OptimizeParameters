eset = eset[, which(colnames(eset) %in% phenodata$ID)]
pData(eset)$Group = phenodata[match(colnames(eset), phenodata$ID, nomatch = 0), 2]

fit = lmFit(eset, design)
cont.matrix = makeContrasts(contrast = CASE - CTRL,  levels = design)
fit = contrasts.fit(fit, cont.matrix)
fit = eBayes(fit)
volc_all = topTable(fit, coef = "contrast", number  = nrow(eset), adjust  = "none", p.value = 1, lfc = 0)
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
  index_topall     = which(rownames(eset) %in% rownames(topall))
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
dir.create(paste(outputPath, backMethod, "limma_Results", paste("lfc_", lfc_exp, sep = ""), sep = "/"), showWarnings = FALSE, recursive = TRUE)
write.table(topall_res, file = paste(outputPath, backMethod, "limma_Results", paste("lfc_", lfc_exp, sep = ""), "dif_exp_results.csv", sep = "/"), row.names = FALSE, sep = ",")