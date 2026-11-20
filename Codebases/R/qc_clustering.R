assign_cluster_resolution <- function(seurat, res, name = "seurat_clusters") {
  require(Seurat)
  seurat <- SetIdent(
    seurat,
    value = paste0("integrated_snn_res.", res)
    )
  seurat <- AddMetaData(
    object = seurat,
    metadata = Idents(seurat),
    col.name = name
  )
  return(seurat)
}

# Output top DGE markers per cluster
top_markers <- function(list_markers, cluster, n = 10){
  df <- list_markers[[cluster]]
  cols <- grep(colnames(df), pattern = "avg_log2FC", value = TRUE)
  markers <- df %>%
    filter(if_all(cols, ~ . > 0)) %>%
    rownames() %>%
    head(n) %>%
    str_flatten(collapse = " ")
  return(markers)
}
