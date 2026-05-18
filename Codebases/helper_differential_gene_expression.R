# Function for running pseudobulk differential expression testing on cell types/groups
# Brad Rosenberg
# 2022-01-17
# Analysis strategy adapted from Amezquita et al
# Multi-Sample Single-Cell Analyses with Bioconductor
# https://bioconductor.org/books/3.14/OSCA.multisample/multi-sample-comparisons.html

prep_pseudobulk <- function(
  seurat, #Seurat object
  ident, # Seurat object ident from which to select cellgroup
  populations = NULL, # the cell groups to include in differential testing
  condition = NULL, # metadata entry within which to contrast
  block = NULL, # blocking factor (e.g. donor series)
  mincells = 10, # Minimum number of cells per sample per cellgroup to include
  min_prop = 0 # Gene expresed by minimum proportion of cells for given population / condition
){
  # Packages required
  require(edgeR)
  require(Seurat)
  require(scran)
  require(scuttle)
  require(dplyr)
  require(tidyverse)

  # Set identity class to the level at which Seurat object will be subset
  Idents(seurat) <- ident

  # Subset to only populations being evaluated
  if(!is.null(populations)){
    seurat <- subset(seurat, idents = populations)
  }

  # Drop any leftover levels
  seurat@meta.data[,ident] <-
    droplevels(seurat@meta.data[,ident])

  # Ensure cell populations names are syntactically valid for edgeR
  seurat@meta.data[,ident] <-
    factor(make.names(seurat@meta.data[,ident]))
  # Set populations object
  populations <- levels(seurat@meta.data[,ident])

  # normalize count data for seurat object
  DefaultAssay(seurat) <- "RNA"
  seurat <- NormalizeData(seurat)

  # Convert seurat to singlecellexperiment object
  DefaultAssay(seurat) <- "RNA"
  sce <- as.SingleCellExperiment(seurat, assay = "RNA")
  colData(sce) <- droplevels(colData(sce))

  # Identify genes to filter by proportion detected
  prop_exp <- summarizeAssayByGroup(
    sce,
    ids = colData(sce)[,c(block, ident, condition)],
    statistics =  "prop.detected",
    store.number = "ncells",
    threshold = 0
  )

  ident_to_test <- levels(factor(colData(prop_exp)[,ident]))
  condition_to_test <- levels(factor(colData(prop_exp)[,condition]))
  genes_to_keep <- lapply(ident_to_test, function(x){
    se <- prop_exp[, prop_exp[[ident]] == x]
    expressed <- lapply(condition_to_test, function(y){
      se <- se[, se[[condition]] == y]
      expressed <- rowSums(assay(se) > min_prop) == ncol(assay(se))
      return(expressed)
    })
    df_expressed <- do.call(cbind, expressed)
    genes_to_keep <- rownames(df_expressed[rowSums(df_expressed) >= 1,])
    return(genes_to_keep)
  })
  genes_to_keep <- unlist(genes_to_keep) %>% unique()

  # Create pseudobulk samples
  summed <- aggregateAcrossCells(
    sce,
    id = colData(sce)[,c(block, ident, condition)]
  )

  # Filter genes at below proporiton expressed threshold per population/condition
  summed <- summed[genes_to_keep,]

  # Filter out low abundance celltype groups
  discarded <- summed$ncells < mincells
  summed <- summed[, !discarded]
  colData(summed) <- droplevels(colData(summed))
  populations <- populations[populations %in%
                               levels(colData(summed)[,ident])
  ]
  colData(summed)[,block] <- factor(colData(summed)[,block])
  colData(summed)[,condition] <- factor(colData(summed)[,condition])
  return(summed)
}

run_standard_edgeR_fit <- function(counts, group){
  y <- DGEList(counts = counts)

  # Filter by expression
  #keep <- filterByExpr(y, group = group)
  #y <- y[keep, , keep.lib.sizes=FALSE]

  # Normalize
  y <- calcNormFactors(y)
  # Dispersion
  y <- estimateDisp(y, design)
  # Fit
  fit <- glmQLFit(y, design)
  return(fit)
}

# "Clip" fold change values for clustering, such that high values don't dominate
clip <- function(x, lower, upper) {
  base::pmax(base::pmin(x, upper), lower)
}

