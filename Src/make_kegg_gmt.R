for(i in 2:length(pathway_id)){
  pathway_id[i] = unlist(strsplit(as.character(pathway_id[i]), split = "-"))[1]
}
for(i in 1:length(pathway_id)){
  pathway_id[i] = str_replace(pathway_id[i], " - Homo sapiens \\(human\\)", "")
}

for(i in 1:length(pathway_id)){
  pathway_id[i] = paste("KEGG_", pathway_id[i], sep = "")
}

for(i in 1:length(pathway_id)){
  pathway_id[i] = str_replace_all(pathway_id[i], ",", "")
}

for(i in 1:length(pathway_id)){
  pathway_id[i] = str_replace_all(pathway_id[i], " ", "_")
}

file = "~/OptimizeParameters/Misc/kegg_gene_sets.gmt"
for(i in 1:length(pathway_id)){
  stringout = paste(pathway_id[i], source_path[i], genes[i], sep = "\t")
  stringout = str_replace_all(stringout, ",", "\t")
  write(stringout, file = file, append = TRUE)
}
