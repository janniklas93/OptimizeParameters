chip_type = eset@annotation

if (chip_type == "hgu133plus2"){  
  hgnc_symbols  = unlist(BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2SYMBOL))
  hgnc_names    = unlist(BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2GENENAME))
  
  hgnc_symbols[is.na(hgnc_symbols)] = ""
  hgnc_names[is.na(hgnc_names)] = ""
  
} else if (chip_type == "hgu133a"){ 
  hgnc_symbols  = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aSYMBOL))
  hgnc_names    = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aGENENAME))
  
  hgnc_symbols[is.na(hgnc_symbols)] = ""
  hgnc_names[is.na(hgnc_names)  ] = ""
}