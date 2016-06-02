gseaOutputPath = paste(outputPath, backMethod, "gseaResults", paste(reshuffling_type, "weighted", weighted_score, sep = "_"), sep = "/")
dir.create(gseaOutputPath, recursive = TRUE, showWarnings = FALSE)
project_ID     = paste(dataSet, backMethod, sep = "_")

GSEA(                                                                      # Input/Output Files :-------------------------------------------
                                                                           input.ds              = paste(celFilesPath, paste("ExpressionSet_", backMethod, ".gct", sep = ""), sep = "/"), # Input gene expression Affy dataset file in RES or GCT format
                                                                           input.cls             = paste(celFilesPath, "phenotypes_GSEA.cls", sep = "/"),                      # Input class vector (phenotype) file in CLS format
                                                                           gs.db                 = geneSetDB_path,                                                             # Gene set database in GMT format
                                                                           output.directory      = gseaOutputPath,                                                             # Directory where to store output and results (default: "")
                                                                           #  Program parameters :----------------------------------------------------------------------------------------------------------------------------
                                                                           doc.string            = project_ID,    # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
                                                                           non.interactive.run   = FALSE,           # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
                                                                           reshuffling.type      = reshuffling_type, # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
                                                                           nperm                 = 1000,            # Number of random permutations (default: 1000)
                                                                           weighted.score.type   =  weighted_score, # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
                                                                           nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
                                                                           fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
                                                                           fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
                                                                           topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
                                                                           adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
                                                                           gs.size.threshold.min = 15,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
                                                                           gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
                                                                           reverse.sign          = FALSE,           # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
                                                                           preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
                                                                           random.seed           = 111,             # Random number generator seed. (default: 123456)
                                                                           perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
                                                                           fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
                                                                           replace               = FALSE,           # For experts only, Resampling mode (replacement or not replacement) (default: F)
                                                                           save.intermediate.results = FALSE,       # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
                                                                           OLD.GSEA              = FALSE,           # Use original (old) version of GSEA (default: F)
                                                                           use.fast.enrichment.routine = TRUE       # Use faster routine to compute enrichment for random permutations (default: T)
)
#--------------------------------------------------------------------------------------------------------------------------------------------------

# Overlap and leading gene subset assignment analysis of the GSEA results                                                                           
GSEA.Analyze.Sets(
  directory = gseaOutputPath,     # Directory where to store output and results (default: "")
  topgs     = 20,                 # number of top scoring gene sets used for analysis
  height    = 16,
  width     = 16
)