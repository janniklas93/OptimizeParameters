celFiles = affy::list.celfiles(celFilesPath, full = TRUE)
rawData = suppressWarnings(read.affybatch(filenames = celFiles))