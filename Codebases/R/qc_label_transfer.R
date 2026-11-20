label_transfer_seurat <- function(seurat,
                                  reference,
                                  query_assay,
                                  ref_assay,
                                  norm_method = "pca",
                                  labels = NULL,
                                  name = "ref"
                                  ) {
  DefaultAssay(reference) <- ref_assay
  anchors <- FindTransferAnchors(
    reference = reference,
    query = seurat,
    dims = 1:30,
    reference.reduction = "pca",
    normalization.method = norm_method,
    query.assay = query_assay,
    reference.assay = ref_assay,
    recompute.residuals = FALSE
  )
  predictions <- TransferData(anchorset = anchors,
                              refdata = reference@meta.data[[labels]],
                              dims = 1:30)
  colnames(predictions) <- paste(name, colnames(predictions), sep = "_")
  seurat <- AddMetaData(seurat, metadata = predictions)
  return(seurat)
}

# Assign major group
# For each cluster, determine the cell type label with the most cells
# Assign this label as the "major group"

assign_major_group <- function(seurat, reference_name, clusters = "seurat_clusters") {
  mat_predictions <- table(unlist(seurat[[clusters]]),
                           unlist(seurat[[reference_name]])) %>%
    as.matrix()
  assignments <- colnames(mat_predictions)[apply(mat_predictions,1,which.max)]
  df_assignments <- data.frame(
    seurat_clusters = as.character(rownames(mat_predictions)),
    group_major = assignments
  )
  # names(df_assignments)[names(df_assignments) == "group_major"] <-
  #   paste("group_major", reference_name, sep = "_")

  # Add annotations to Seurat obj metadata
  df_metadata <- seurat@meta.data
  df_metadata$cell_barcode <- rownames(df_metadata)
  df_metadata <- left_join(
    x = df_metadata,
    y = df_assignments,
    by = "seurat_clusters"
  )
  rownames(df_metadata) <- df_metadata$cell_barcode

  # Add to Seurat object
  seurat <- AddMetaData(seurat, metadata = df_metadata)
  return(seurat)
}
