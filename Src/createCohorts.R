rownames(pData(rawData)) = base::gsub(c(".gz|.CEL|.cel|.GZ"), "", rownames(pData(rawData)))
sampleNames(rawData)     = base::gsub( c(".gz|.CEL|.cel|.GZ"), "", sampleNames(rawData)) 
phenodata$ID = base::gsub(c(".gz|.CEL|.cel|.GZ"), "", phenodata$ID) # the reason is to assure, that no .gz ending is present when we add it
phenodata = phenodata[phenodata$ID %in% rownames(pData(rawData)), ]

cohortsVec = as.vector(phenodata$Group)
cohortsVec[cohortsVec %in% set_ctrl] = "CTRL"
cohortsVec[cohortsVec %in% set_case] = "CASE"
index_cohorts_vec = which(cohortsVec %in% c("CTRL", "CASE") == TRUE)
cohortsVec = cohortsVec[cohortsVec %in% c("CTRL", "CASE")]
names(cohortsVec) = as.character(phenodata$ID[index_cohorts_vec])
nmbr_samples = sum(cohortsVec %in% c("CTRL", "CASE"))

design = model.matrix( ~ 0 + cohortsVec )
colnames(design)[colnames(design) == paste0("cohortsVec", "CTRL")] = "CTRL"
colnames(design)[colnames(design) == paste0("cohortsVec", "CASE")] = "CASE"  

index_ctrl = as.integer(which(design[, colnames(design) == "CTRL" ] == 1))
index_case = as.integer(which(design[, colnames(design) == "CASE" ] == 1))

#### maybe remove ?? ###
raw_data_group_vec = rep("", dim(pData(rawData))[1])
raw_data_group_vec[index_ctrl] = "CTRL"
raw_data_group_vec[index_case] = "CASE"
pData(rawData) = cbind(pData(rawData), raw_data_group_vec)
raw_data_group_vec = raw_data_group_vec[which(raw_data_group_vec != "")]
colnames(pData(rawData))[-1] = "Cohort"
  