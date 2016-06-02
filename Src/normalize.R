eset = threestep(rawData, background.method = backMethod, normalize.method = "quantile", summary.method = "median.polish")

rownames(pData(eset)) = base::gsub(c(".gz|.CEL|.cel|.GZ"), "", rownames(pData(eset)))
mapping_cohort_p      = match(base::gsub(c(".gz|.CEL|.cel|.GZ"), "", rownames(pData(eset))), base::gsub(c(".gz|.CEL|.cel|.GZ"), "", names(cohortsVec)), nomatch = 0)
mapping_cohort_c      = match(base::gsub(c(".gz|.CEL|.cel|.GZ"), "", names(cohortsVec)) , base::gsub(c(".gz|.CEL|.cel|.GZ"), "", rownames(pData(eset))) , nomatch = 0)
eset = eset[, rownames(pData(eset)) %in% base::gsub(c(".gz|.CEL|.cel|.GZ"), "", names(cohortsVec))]
pData(eset)$Cohort = cohortsVec[mapping_cohort_p]