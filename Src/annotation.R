chip_type = eset@annotation

if (chip_type == "hgu133plus2"){  
  hgnc_symbols  = unlist(BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2SYMBOL))
  hgnc_names    = unlist(BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2GENENAME))
  ensembl_genes = unlist(BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2ENSEMBL))
  entrez_genes  = unlist(BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2ENTREZID))
  
  hgnc_symbols[is.na(hgnc_symbols)] = ""
  hgnc_names[is.na(hgnc_names)] = ""
  ensembl_genes[is.na(ensembl_genes)] = ""
  entrez_genes[is.na(entrez_genes)] = ""
  
} else if (chip_type == "hgu133a"){ 
  hgnc_symbols  = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aSYMBOL))
  hgnc_names    = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aGENENAME))
  ensembl_genes = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aENSEMBL))
  entrez_genes  = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aENTREZID))
  
  hgnc_symbols[is.na(hgnc_symbols)] = ""
  hgnc_names[is.na(hgnc_names)  ] = ""
  ensembl_genes[is.na(ensembl_genes)  ] = ""
  entrez_genes[is.na(entrez_genes)] = ""
}