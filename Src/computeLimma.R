eset = eset[, which(colnames(eset) %in% phenodata$ID)]
pData(eset)$Group = phenodata[match(colnames(eset), phenodata$ID, nomatch = 0), 2]

fit = lmFit(eset, design)
cont.matrix = makeContrasts(contrast = CASE - CTRL,  levels = design)
fit = contrasts.fit(fit, cont.matrix)
fit = eBayes(fit)
volc_all = topTable(fit, coef = "contrast", number  = nrow(eset), adjust  = "none", p.value = 1, lfc = 0)