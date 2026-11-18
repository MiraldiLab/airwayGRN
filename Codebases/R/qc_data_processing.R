read_cellranger_data <- function(runID, sample_metadata, cellranger_path){
  require(Seurat)
  require(dplyr)

  # Read in 10X GEX and HTO output
  ## per_sample_outs/ contains FILTERED barcodes
  mat_10X <- Read10X(
    file.path(cellranger_path, runID, "outs/count/filtered_feature_bc_matrix")
  )

  ## IAV features to remove
  features_to_remove <-
    grep(rownames(mat_10X), pattern = "^\\+|^\\-", value = TRUE)
  mat_10X$`Gene Expression` <-
    mat_10X$`Gene Expression`[!(rownames(mat_10X$`Gene Expression`) %in%
                                  features_to_remove
    ), ]

  # Create Seurat object with GEX data
  seurat <- CreateSeuratObject(
    counts = mat_10X$`Gene Expression`,
    project = runID
  )

  # Remove unused MULTISeq HTO barcodes for each pool
  HTO_counts <- mat_10X$`Antibody Capture`
  HashIDs <- sample_metadata %>%
    filter(library == runID) %>%
    dplyr::pull("MS_barcode")

  ## Remove count data for HTOs not represented in the library pool
  HTO_counts <- HTO_counts[HashIDs, , drop = FALSE]

  # Add HTO data as additional assay in Seurat object
  seurat[["HTO"]] <- CreateAssayObject(counts = HTO_counts)

  # Return Seurat object
  return(seurat)
}

read_cellranger_metrics <- function(runIDs, cellranger_path){
  cellranger_metrics_list <- lapply(runIDs, function(runID){
    # Import data from 10X files
    path_cellranger_metrics <- file.path(
      cellranger_path, runID, "outs/count/metrics_summary.csv"
      )
    cellranger_metrics <- read.csv(file = path_cellranger_metrics, header = TRUE, check.names = TRUE)
    names_metrics <- colnames(cellranger_metrics)
    # names_metrics <- str_remove_all(colnames(cellranger_metrics), "[.]")
    cellranger_metrics <- cellranger_metrics %>%
    gsub(pattern = "%|,",  replacement = "") %>%
    as.numeric()
    names(cellranger_metrics) <- names_metrics
    return(cellranger_metrics)
  })
  df_cellranger_metrics <- as.data.frame(do.call(rbind, cellranger_metrics_list))
  df_cellranger_metrics$pool <- runIDs
  return(df_cellranger_metrics)
}

plot_cellranger_metrics <- function(df, metric){
  ggplot(df, aes(x = pool, y = get(metric))) +
    geom_bar(stat = "identity") +
    ylab(
      gsub(x = metric, pattern = "Antibody.", replacement = "") %>%
      gsub(x = ., pattern = "..summed.over.all.recognized.antibody.barcodes.",
           replacement = "") %>%
      gsub(x = ., pattern = "\\.", replacement = " ") %>%
      str_wrap(., width = 15)
      ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size =8),
      aspect.ratio = 1
    )
  }

plot_violin_percell_metrics <- function(seurat = seurat, metric, group = "orig.ident"){
  VlnPlot(seurat, features = metric, group.by = group, pt.size = 0) +
    geom_point(size = 1, alpha = 0.05, position = "jitter") +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(size = 6),
      title = element_text(size = 8),
      aspect.ratio = 1
    )
  }

demultiplex_by_HTO <- function(seurat, sample_metadata){
  library <- seurat@project.name
  quantile <- 0.9 # Set empirically for this project
  # Demultiplex
  ## TODO: allow for selection of different demultiplex functions
  seurat <- HTODemux(seurat, assay = "HTO",
                positive.quantile = quantile, verbose = FALSE)

  ## Assemble merged metadata table outside of Seurat object
  ## Subset metadata to current library pool
  sample_metadata <- filter(sample_metadata, library == seurat@project.name) %>%
    droplevels()

  seurat$cell_barcode <- rownames(seurat@meta.data)
  df_metadata <- seurat@meta.data

  ## Assign sample metadata based on HTO
  df_metadata <- left_join(
    x = df_metadata,
    y = sample_metadata,
    by = c("hash.ID" = "MS_barcode")
  )
  rownames(df_metadata) <- df_metadata$cell_barcode

  # Add to Seurat object
  seurat <- AddMetaData(seurat, metadata = df_metadata)

  # Reorder classification factor for downstream plotting
  seurat$HTO_classification.global <- factor(seurat$HTO_classification.global,
                                             levels = c("Negative",
                                                        "Singlet",
                                                        "Doublet"))

  # Return single run Seurat object with demultiplex assignments and metadata
  return(seurat)
}

run_tSNE_HTO <- function(seurat, remove.negatives = TRUE){
  Idents(seurat) <- "hash.ID"
  hto_tsne <- seurat

  # Remove droplets designated "Negative" by demultiplexing?
  if(remove.negatives == TRUE){
    hto_tsne <- subset(hto_tsne, idents = "Negative", invert = TRUE)
    }

  # Calculate a distance matrix using HTO
  hto.dist.mtx <- as.matrix(
    dist(t(GetAssayData(object = hto_tsne, assay = "HTO")))
  )

  # Calculate tSNE embeddings with a distance matrix
  hto_tsne <- RunTSNE(hto_tsne, distance.matrix = hto.dist.mtx,
                      perplexity = 100
  )
  return(hto_tsne)
}

detect_doublets <- function(seurat, BPPARAM = BiocParallel::bpparam()){
  # Convert to sce object for scDblFinder workflow
  DefaultAssay(seurat) <- "RNA"
  ##  Add metadata column for HTO classification (logical)
  seurat$Doublet <- seurat$HTO_classification.global == "Doublet"
  sce <- as.SingleCellExperiment(seurat, assay = "RNA")

  # Run scDblFinder
  ## Primary use defaults
  sce <- scDblFinder(
    sce = sce,
    clusters = FALSE, # Generate artificial doublets
    #dbr = dblrate,
    knownDoublets = "Doublet", # Specify HTO classifications
    BPPARAM = BPPARAM
  )
  # Add scDblFinder results as metadata back to Seurat object
  df_metadata <- as.data.frame(
    colData(sce)[,c("ident", "scDblFinder.class", "scDblFinder.score",
                    "scDblFinder.weighted", "scDblFinder.cxds_score")
    ]
  )
  seurat <- AddMetaData(seurat, metadata = df_metadata)
  return(seurat)
}
