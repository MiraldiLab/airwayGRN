# ==============================================================================
# scRNA-Seq quality control, demultiplexing, clustering, and cell type annotation
# Alec P. Pankow & Brad R. Rosenberg
# Icahn School of Medicine at Mount Sinai, New York, NY USA
# ==============================================================================

# ==============================================================================
# 1. SETUP
# ==============================================================================

# ---- R version check ---------------------------------------------------------

r_version <- "4.0.3"
if (r_version == paste0(R.Version()$major,
  sep = ".",
  R.Version()$minor
)) {
  print(paste(
    "Running R Version",
    r_version
  ))
} else {
  print(paste(
    "To help ensure reproducibility, please use R Version",
    r_version
  ))
}

# ---- Packages ----------------------------------------------------------------
library(bookdown)
library(here)
library(knitr)
library(future)
library(future.apply)
library(devtools)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(magrittr)
library(git2r)
library(patchwork)
library(pals)
library(clustree)
library(glmGamPoi)
library(multtest)
library(metap)
library(SingleR)
library(rgeos)
library(Seurat)
library(viridis)

# Bioconductor packages
options(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = FALSE)
library(BiocManager)
library(scDblFinder)

# ---- RNG seed ----------------------------------------------------------------
set.seed(42)

# ---- Parallelization (future / BiocParallel) ---------------------------------
# Set number of cores on compute instance here
ncores <- 5 # One per donor series
plan("multisession", workers = ncores)
# Set memory caps here
mem_max <- 32 # 32 GB
options(future.globals.maxSize = mem_max * 1024 ^ 3)

# BiocParallel
bpparam <- BiocParallel::bpparam()
bpparam$workers <- ncores 
bpparam$RNGseed <- 42

# ---- Paths -------------------------------------------------------------------
# UPDATE THIS PATH TO POINT TO CELLRANGER MULTI OUTPUT
# SUBDIRS SHOULD BE RUNIDS IN `metadata/scrna.csv` AND CONTAIN
# `/outs/count/filtered_feature_bc_matrix/`
cellranger_multi_path <- "<path-to-cellranger-multi-output>"

source(here::here("R/qc_data_processing.R"))
source(here::here("R/qc_clustering.R"))
source(here::here("R/qc_label_transfer.R"))

# Generic violin plots
vln_theme <- theme(
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.text = element_text(size = 6),
  title = element_text(size = 8)
)

# ==============================================================================
# 2. READ IN DATASETS
# ==============================================================================

# ---- Sample metadata ---------------------------------------------------------
sample_metadata <- read.csv(
  here::here("metadata/scrna.csv"),
  header = TRUE
)
# Set a sampleID column
sample_metadata$sampleID <- 
  paste(sample_metadata$treatment,
        sample_metadata$sample,
        sample_metadata$donor,
        sep = "_")

# Set a compound factor for timepoint_donor
sample_metadata$donor_timepoint <-
  paste(sample_metadata$donor, sample_metadata$timepoint, sep = "_")

# Set a compound factor for timepoint_donor_series
sample_metadata$donor_timepoint_series <-
  paste(sample_metadata$donor_timepoint, sample_metadata$series, sep = "_")

# Convert all columns to factors
sample_metadata %<>% mutate_at(colnames(sample_metadata), factor)

# Output table
knitr::kable(sample_metadata)

# ---- 10X gene expression + HTO data ------------------------------------------
runs_10X <- list.dirs(path = cellranger_multi_path, 
                      recursive = FALSE, 
                      full.names = FALSE) %>%
  grep(pattern = "^IFN_|^BRRO04_", value = TRUE)

# Read in each 10X run to Seurat object in list
seurat_list <- future_lapply(runs_10X, read_cellranger_data,
                             sample_metadata = sample_metadata,
                             cellranger_path = cellranger_multi_path,
                             future.seed = TRUE)
names(seurat_list) <- runs_10X

# ==============================================================================
# 3. CELLRANGER METRICS
# ==============================================================================

# Read in cellranger metrics
df_cellranger_metrics <- read_cellranger_metrics(
  runIDs = runs_10X,
  cellranger_path = cellranger_multi_path
)

# Select GEX metrics to plot
metrics_to_plot_gex <- c(
  "Mean.Reads.per.Cell",
  "Median.Genes.per.Cell",
  "Sequencing.Saturation",
  "Q30.Bases.in.RNA.Read",
  "Reads.Mapped.Confidently.to.Transcriptome",
  "Median.UMI.Counts.per.Cell"
)

p_gex_list <- lapply(metrics_to_plot_gex,
                     plot_cellranger_metrics,
                     df = df_cellranger_metrics)

# Select HTO metrics to plot
metrics_to_plot_hto <- c( 
  "Antibody..Mean.Reads.per.Cell",
  "Antibody..Sequencing.Saturation",
  "Antibody..Q30.Bases.in.Antibody.Read",
  "Antibody..Median.UMIs.per.Cell..summed.over.all.recognized.antibody.barcodes."
)

p_hto_list <- lapply(metrics_to_plot_hto,
                     plot_cellranger_metrics,
                     df = df_cellranger_metrics)
# Plot GEX metrics
wrap_plots(p_gex_list, nrow = 2) +
  plot_annotation(title = "CellRanger Gene Expression run metrics")

# Plot MULTISeq HTO metrics
wrap_plots(p_hto_list, ncol = 3) +
  plot_annotation(title = "CellRanger MULTISeq HTO run metrics")

# ==============================================================================
# 4. INITIAL QUALITY ASSESSMENT AND FILTERING
# ==============================================================================

# Select cell metrics to plot
per_cell_metrics_to_plot <- c("nFeature_RNA", "nCount_RNA", "nCount_HTO")

# Merge seurat objects for plotting purposes
seurat_merged <- merge(seurat_list[[1]], y = seurat_list[c(-1)])

vln_list <- lapply(per_cell_metrics_to_plot, 
                   plot_violin_percell_metrics, seurat = seurat_merged)
names(vln_list) <- per_cell_metrics_to_plot

# Plot specific adjustments      
min_umi <- 4000
vln_list[["nCount_RNA"]] <- vln_list[["nCount_RNA"]] +
  scale_y_continuous(trans = "log10") +
  geom_hline(yintercept = min_umi, linetype = 2, color = "red")

vln_list[["nCount_HTO"]] <- vln_list[["nCount_HTO"]] +
  scale_y_continuous(trans = "log10")

# Plot
wrap_plots(vln_list, nrow = 1, guides = "collect") +
  plot_annotation(title = "Initial QC and filtering prior to demultiplexing")

# ==============================================================================
# 5. HTO DEMULTIPLEXING
# ==============================================================================

# HTO labeling was performed using MULTI-Seq (McGinnis et al. 2019)

# ---- Initial UMI filter ------------------------------------------------------
# Perform initial filtering on nCount_RNA
seurat_list <- lapply(seurat_list, subset, nCount_RNA < min_umi, invert = TRUE)

# ---- Normalize HTO data ------------------------------------------------------
seurat_list <- future_lapply(seurat_list, future.seed = TRUE, NormalizeData,
                             assay = "HTO", normalization.method = "CLR",
                             margin = 2)

# ---- HTODemux ----------------------------------------------------------------
# Demultiplex with HTOdemux
seurat_list <- future_lapply(seurat_list, future.seed = TRUE,
                             demultiplex_by_HTO,
                             sample_metadata = sample_metadata)

# ---- tSNE on HTO data (per library) ------------------------------------------
# This can take a long time, but can be informative as a per-sample summary
hto_tsne_list <- future_lapply(seurat_list, future.seed = TRUE,
                               run_tSNE_HTO, remove.negatives = TRUE)

hto_tsne_plots_list <- lapply(hto_tsne_list, function(x){
  DimPlot(x, group.by = "hash.ID", pt.size = 0.00001) +
    ggtitle(x@project.name) +
    scale_color_manual(values=cols25()) +
    theme(
      axis.text = element_text(size = 4),
      title = element_text(size = 6),
      legend.text = element_text(size = 4),
      aspect.ratio = 1
      )
  })
  
wrap_plots(hto_tsne_plots_list) +
  plot_annotation(title = "Initial HTO demultiplexing by tSNE")

# ---- HTODemux margin filter --------------------------------------------------
# Optimal margin cutoff previously determined by Emma DeGrace analyses in lab
margin_cutoff <- 0.25
hto_tsne_margin_plots_list <- lapply(hto_tsne_list, function(x){
  FeaturePlot(x, features = "HTO_margin", max.cutoff = margin_cutoff) +
    ggtitle(x@project.name) +
    scale_color_viridis() +
    theme(
      axis.text = element_text(size = 4),
      title = element_text(size = 6),
      legend.text = element_text(size = 4),
      aspect.ratio = 1
      )
})
wrap_plots(hto_tsne_margin_plots_list, guides = "collect") +
  plot_annotation(title = paste0("Minimum HTODemux margin set to ", margin_cutoff))

# Remove cells by HTODemux margin cutoff
seurat_list <- lapply(seurat_list, subset, HTO_margin > margin_cutoff)

# ---- Post-filter tSNE --------------------------------------------------------
# Regenerate tSNE embeddings
hto_tsne_list <- future_lapply(seurat_list, future.seed = TRUE,
                               run_tSNE_HTO, remove.negatives = TRUE)

hto_tsne_plots_list <- lapply(hto_tsne_list, function(x){
  DimPlot(x, group.by = "hash.ID", pt.size = 0.00001) +
    ggtitle(x@project.name) +
    scale_color_manual(values=cols25()) +
    theme(
      axis.text = element_text(size = 4),
      title = element_text(size = 6),
      legend.text = element_text(size = 4),
      aspect.ratio = 1
      )
  })
  
wrap_plots(hto_tsne_plots_list) +
  plot_annotation(
    title = "HTO demultiplexing (post HTODemux margin filter) by tSNE")

# ---- HTO heatmaps ------------------------------------------------------------
hm_list <- lapply(seurat_list, function(x) {
  Idents(x) <- "hash.ID"
  hm <- HTOHeatmap(x, assay = "HTO", classification = "hash.ID") +
    scale_fill_viridis(option = "magma", direction = 1) +
    ggtitle(x@project.name) +
    theme(
      axis.text = element_text(size = 6),
      title = element_text(size = 8), 
    ) + NoLegend()
  return(hm)
})

wrap_plots(hm_list, ncol = 2)

# ==============================================================================
# 6. POST-DEMUX QC AND FILTERING
# ==============================================================================

# Add percent mitochondrial metadata
seurat_list <- lapply(seurat_list, function(x){
  DefaultAssay(x) <- "RNA"
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  return(x)
})

# Select cell metrics to plot
per_cell_metrics_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

# Generate violin plots
vln_list <- lapply(per_cell_metrics_to_plot, function(x){
  p_list <- lapply(seurat_list, function(y){
    p <- plot_violin_percell_metrics(seurat = y, metric = x, 
                                     group = "HTO_classification.global") +
      ggtitle(y@project.name) +
      ylab(x)
  })
  return(p_list)
})
names(vln_list) <- per_cell_metrics_to_plot
                   
# Plot specific adjustments 
vln_list[["nFeature_RNA"]] <- lapply(vln_list[["nFeature_RNA"]], function(p){
  p + scale_y_continuous(limits = c(0, 8000)) +
    theme(axis.text.x = element_blank())
})


min_umi_postdemux <- 5500
vln_list[["nCount_RNA"]] <- lapply(vln_list[["nCount_RNA"]], function(p){
  p + scale_y_continuous(trans = "log10", limits = c(min_umi, 100000)) +
    geom_hline(yintercept = min_umi_postdemux, linetype = 2, color = "red") +
    ggtitle(NULL) +
    theme(axis.text.x = element_blank())
})

max_mito <- 20
vln_list[["percent.mt"]] <- lapply(vln_list[["percent.mt"]], function(p){
  p + scale_y_continuous(limits = c(0, 90)) +
    geom_hline(yintercept = max_mito, linetype = 2, color = "red") +
    ggtitle(NULL)
})

# Plot
vln_list <- unlist(vln_list, recursive = FALSE)
wrap_plots(vln_list, nrow = 3, guides = "collect") & 
  theme(legend.position = "bottom") 

seurat_list <- future_lapply(seurat_list, future.seed = TRUE, function(x){
  x <- subset(x, nCount_RNA > min_umi_postdemux & percent.mt < max_mito)
  x <- subset(x, HTO_classification.global == "Negative", invert = TRUE)
  return(x)
})

# ==============================================================================
# 7. DOUBLET DETECTION
# ==============================================================================

seurat_list <- lapply(seurat_list, detect_doublets, BPPARAM = bpparam)

# Remove doublets annotated by both hashing demultiplexing and scDblFinder
seurat_list <- future_lapply(seurat_list, future.seed = TRUE, function(x){
  Idents(x) <- "HTO_classification.global"
  x <- subset(x, idents = "Doublet", invert = TRUE)
  Idents(x) <- "scDblFinder.class"
  x <- subset(x, idents = "doublet", invert = TRUE)
  return(x)
})

# ---- Post-filter QC violins --------------------------------------------------
vlnfeature_list <- lapply(seurat_list, function(x){
  Idents(x) <- "hash.ID"
  p <- VlnPlot(x, features = "nFeature_RNA", y.max = 8000, pt.size = 0) +
    ggtitle(unique(x@meta.data$orig.ident)) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank()
      ) 
  return(p)
})
wrap_plots(vlnfeature_list, nrow = 2) & theme(text = element_text(size = 10))

# Visualize nCount RNA
vlncount_list <- lapply(seurat_list, function(x){
  Idents(x) <- "hash.ID"
  p <- VlnPlot(x, features = "nCount_RNA", log = FALSE, pt.size = 0) +
    scale_y_continuous(limits = c(min_umi, 70000), trans = "log10") +
    ggtitle(unique(x@meta.data$orig.ident)) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank()
      ) 
  return(p)
})
wrap_plots(vlncount_list, nrow = 2) & theme(text = element_text(size = 10))

# Visualize mitochondrial RNA
vlnmito_list <- lapply(seurat_list, function(x) {
  Idents(x) <- "hash.ID"
  p <- VlnPlot(x, features = "percent.mt", pt.size = 0) +
    ggtitle(unique(x@meta.data$orig.ident)) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank()
    )
  return(p)
})
wrap_plots(vlnmito_list, nrow = 2) & theme(text = element_text(size = 10))

# ==============================================================================
# 8. ASSESS BATCH AND RELATED SURROGATE VARIABLE EFFECTS
# ==============================================================================

# Merge all Seurat objects
seurat_merged <- merge(seurat_list[[1]], y = seurat_list[c(-1)])

# Conserve memory parameter makes a big difference in getting this to run
# Too much memory consumed due to parallelization...must be set to 1 core
plan("sequential")
seurat_merged <- SCTransform(seurat_merged,
                             vst.flavor = "v2",
                             vars.to.regress = "percent.mt", 
                             conserve.memory = TRUE,
                             verbose = FALSE)
seurat_merged <- RunPCA(seurat_merged)
# ElbowPlot(seurat_merged, ndims = 50)
# DimHeatmap(seurat_merged, dims = 28:36, cells = 250, balanced = TRUE)
# Structure starts declining by around PC30
ndims <- 30
seurat_merged <- RunUMAP(seurat_merged, dims = 1:ndims, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:ndims)

# pool effect?
p_pool <- DimPlot(seurat_merged, group.by = "orig.ident", shuffle = TRUE)
# series effect?
p_series <- DimPlot(seurat_merged, group.by = "series", shuffle = TRUE)
# BC HTO effects?
p_hto <- DimPlot(seurat_merged, group.by = "hash.ID", split.by = "donor_timepoint_series", shuffle = TRUE, ncol = 3) & NoLegend()

# ==============================================================================
# 9. DATA INTEGRATION
# ==============================================================================

plan("multisession", workers = ncores)

# ---- Split libraries by sample -----------------------------------------------
# Split libriaries by sample
timepoint_list <- future_lapply(seurat_list, future.seed = TRUE, function(x){
  pool_list <- SplitObject(x, split.by = "donor_timepoint_series")
}) %>% unlist()

# ---- Per-sample SCTransform --------------------------------------------------
# SCTransform
timepoint_list <- future_lapply(timepoint_list, future.seed = TRUE, function(x){
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(x, 
                   vst.flavor = "v2",
                   variable.features.n = 4000,
                   verbose = FALSE,
                   return.only.var.genes = FALSE, 
                   vars.to.regress = "percent.mt")
})

# ---- Cell cycle scoring ------------------------------------------------------
# Add cell cycle scores PER SAMPLE prior to integration
timepoint_list <- future_lapply(timepoint_list, future.seed = TRUE,
                                CellCycleScoring,
                                assay = "SCT",
                                s.features = cc.genes$s.genes,
                                g2m.features = cc.genes$g2m.genes)

# ---- Integration features and anchors (RPCA) ---------------------------------
# Sequential processing for memory requirements
plan("sequential")

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(
  object.list = timepoint_list, nfeatures = 4000,
  assay = rep.int("SCT", times = length(timepoint_list)),
  verbose = FALSE
)

timepoint_list <- PrepSCTIntegration(
  object.list = timepoint_list,
  anchor.features = features,
  assay = "SCT",
  verbose = FALSE
)

# Run PCA on each dataset using these features
timepoint_list <- future_lapply(timepoint_list, future.seed = TRUE, function(x){
  x <- RunPCA(x, verbose = FALSE)
})

# Identify anchors using time 0 hpi samples as the reference
ref_samples <- names(timepoint_list) %>% grep(pattern = "_0_", value = TRUE)

anchors <- FindIntegrationAnchors(
  object.list = timepoint_list,
  reference = which(names(timepoint_list) %in% ref_samples),
  normalization.method = "SCT",
  anchor.features = features,
  dims = 1:30, # Seurat Default
  k.anchor = 10, # Tested several values; results generally similar
  reduction = "rpca",
  verbose = FALSE
)

# Integrate all samples from each library
seurat_integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  dims = 1:30,
  verbose = TRUE
)

# ---- PCA / UMAP on integrated data -------------------------------------------
# Run PCA/UMAP on integrated data
DefaultAssay(seurat_integrated) <- "integrated"
seurat_integrated <- RunPCA(
  seurat_integrated, npcs = 100,
  verbose = FALSE
)
ndims <- 50

plan("multisession", workers = ncores)
seurat_integrated <- RunUMAP(seurat_integrated,
                          reduction = "pca",
                          dims = 1:ndims, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:ndims, verbose = FALSE)

factors_to_plot <- c("series", "library", "donor",
                     "timepoint", "donor_timepoint_series")
p_umap_list <- lapply(factors_to_plot, function(x){
  p <- DimPlot(seurat_integrated, group.by = x,
               shuffle = TRUE, pt.size = 0.00001) +
    ggtitle(x) +
    scale_color_manual(values = cols25()) +
    theme(
      legend.position = "right"
    ) + coord_fixed()
  return(p)
  })
wrap_plots(p_umap_list, ncol = 2)

# ==============================================================================
# 10. CLUSTERING (FIRST ROUND)
# ==============================================================================

seurat_integrated <- FindClusters(
  seurat_integrated,
  resolution = seq(from = 0.2, to = 3, by = 0.2),
  verbose = FALSE
)
clustree(seurat_integrated, prefix = "integrated_snn_res.")

ref_res <- 2.2 # Not that stable but err on higher resolution for initial QC
seurat_integrated <- assign_cluster_resolution(seurat_integrated, res = ref_res)

p_umap_clus <- DimPlot(
  seurat_integrated,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = FALSE
) + 
  NoLegend() +
  ggtitle("First round clustering") +
  coord_fixed()
p_umap_clus

# ---- Cluster representation across samples -----------------------------------
# by donor_timepoint_series
p <- ggplot(
  seurat_integrated@meta.data,
  aes(x = seurat_clusters, fill = donor_timepoint_series)
) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols25()) + 
  theme(
    rect = element_blank(),
    axis.line = element_line(),
    panel.grid.major.y = element_line(color = "grey"),
    panel.grid.minor.y = element_line(color = "grey"),
    legend.position = "bottom"
  ) +
  xlab("Unannotated Clusters") +
  ylab("Fraction of Cluster")
p

# by orig.ident
p <- ggplot(
  seurat_integrated@meta.data,
  aes(x = seurat_clusters, fill = orig.ident)
) +
  geom_bar(position = "fill") +
   scale_fill_manual(values = cols25()) + 
  theme(
    rect = element_blank(),
    axis.line = element_line(),
    panel.grid.major.y = element_line(color = "grey"),
    panel.grid.minor.y = element_line(color = "grey"),
    legend.position = "bottom"
  ) +
  xlab("Unannotated Clusters") +
  ylab("Fraction of Cluster")
p

# ==============================================================================
# 11. MARKER GENES AND CLUSTER QC (SECOND ROUND)
# ==============================================================================

# ---- Canonical marker genes --------------------------------------------------
# We evaluate general cell populations to “get oriented” by expression of established marker genes.
# Ciliated cells: FOXJ1, CAPS
# FOXJ1, along with other transcription factors such as RFX and MYB, activate the transcriptional program for assembly of cilia (DOI:10.1513/AnnalsATS.201802-128AW).
# CAPS, a calcium-binding protein, is an additional marker gene used to identify ciliated cells (DOI:10.1038/s41577-020-00477-9).
# Basal cells: KRT5, TP63
# TP63 and KRT5 are used as marker genes for basal cells (DOI:10.1513/AnnalsATS.201802-128AW, DOI:10.1038/s41577-020-00477-9, DOI:10.1016/j.tcb.2018.04.007).
# TP63 is a transcription factor used to define basal cells which can serve as progenitors of other cell types such as goblet and ciliated cells.
# Secretory cells: SCGB1A1, MUC5B, BPIFA1
# Mucins are highly glycosylated proteins produced by airway epithelial cells and submucosal glands, the most abundant of which are MUC5B and MUC5AC.
# While SCGB1A1 has historically been used as a club cell-specific marker (DOI:10.1513/AnnalsATS.201802-128AW, DOI:10.1038/s41577-020-00477-9), more recent RNA-FISH studies of human airway epithelium have found it to be expressed more broadly in secretory cells, depending on the region of the airway/bronchial tree sampled (DOI: 10.1164/rccm.201804-0734OC).
# Ionocytes: CFTR, ASCL3
# CFTR, a chloride-ion transporter, and ASCL3, a helix-loop-helix transcription factor, were recentlly identified as markers of ionocytes (DOI:10.1038/s41586-018-0394-6, DOI:10.1038/s41586-018-0393-7, DOI:10.1038/s41577-020-00477-9).
# Neuroendocrine: CHGA, CALCA
# Neuroendocrine (NE) cells are defined by ASCL1 expression, a transcription factor that is essential for the formation of pulmonary neuroendocrine cells, also know as PNECs (DOI:10.1038/s41577-020-00477-9).
# PNECs and pulmonary ionocytes are suggested to be developmentally related, as both arise during differentiation from POU2F3+ tuft-like cells (DOI:10.1038/s41467-020-16239-z, DOI:10.1038/s41577-020-00477-9).

# Make vectors of marker genes
ciliated <- c("FOXJ1", "CAPS")
basal <- c("KRT5", "TP63")
secretory <- c("SCGB1A1", "MUC5B", "BPIFA1")
ionocytes <- c("CFTR", "ASCL3")
neuroendocrine <- c("CHGA", "CALCA")
tuftlike <- c("POU2F3")
marker_genes <- c(ciliated, basal, secretory,
                  ionocytes, neuroendocrine, tuftlike)

# Plot Marker Gene Expression
DefaultAssay(seurat_integrated) <- "SCT"
p_umap_markers <- FeaturePlot(
  seurat_integrated,
  features = marker_genes,
  order = FALSE,
  ncol = 3
)  &
  scale_color_viridis(option = "magma") &
  theme(
    text = element_text(size =  10),
    axis.text = element_blank(),
    axis.title = element_blank()
  )
p_umap_clus + p_umap_markers

# ---- Cluster-based QC --------------------------------------------------------

# Based on data exploration
t0_isghidoublet_clusters <- c("27")

# Based on gene expression exploration and related doublet scoring data
potential_doublet_clusters <- c("6", "13")

# UMI counts
p <- VlnPlot(
  seurat_integrated,
  features = "nCount_RNA",
  pt.size = 0,
  y.max = 50000,
  group.by = "seurat_clusters",
  log = TRUE
) +
  xlab("Unannotated Clusters") +
  labs(
    title = "Distribution of UMI Counts by Cluster",
    caption = "Log10 transformed Y Axis"
  ) +
  NoLegend()
p

# % mito
# Percent mito
p <- VlnPlot(
  seurat_integrated,
  features = "percent.mt",
  pt.size = 0,
  group.by = "seurat_clusters",
  log = TRUE
) +
  xlab("Unannotated Clusters") +
  labs(
    caption = "Log10 transformed Y Axis"
  ) +
  NoLegend()
p
hi_mito_clusters <- c("32", "35")

clusters_to_exclude <- c(hi_mito_clusters, t0_isghidoublet_clusters, potential_doublet_clusters)
Idents(seurat_integrated) <- "seurat_clusters"
seurat_integrated <- subset(seurat_integrated,
                            idents = clusters_to_exclude, invert = TRUE)

# ---- Second round clustering -------------------------------------------------

# Run PCA/UMAP on integrated data
DefaultAssay(seurat_integrated) <- "integrated"
seurat_integrated <- RunPCA(
  seurat_integrated, npcs = 100,
  verbose = FALSE
)
ndims <- 100 
# ElbowPlot(seurat_integrated, ndims = 100) + 
# geom_vline(xintercept = ndims, color = "#ED3438")
# DimHeatmap(seurat_integrated, dims = (ndims - 4):(ndims + 4), 
# balanced = TRUE, cells = 250)

plan("multisession", workers = ncores)
seurat_integrated <- RunUMAP(seurat_integrated,
                          reduction = "pca",
                          dims = 1:ndims, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:ndims, verbose = FALSE)

seurat_integrated <- FindClusters(
  seurat_integrated,
  resolution = seq(from = 0.2, to = 3, by = 0.2),
  verbose = FALSE
)
clustree(seurat_integrated, prefix = "integrated_snn_res.")

ref_res <- 1.4
seurat_integrated <- assign_cluster_resolution(seurat_integrated, res = ref_res)

p_umap_clus <- DimPlot(
  seurat_integrated,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = FALSE
) + 
  NoLegend() + 
  coord_fixed() +
  labs(title = "UMAP Clustering without Annotation")
p_umap_clus

# ---- Cell cycle scores by cluster --------------------------------------------
p1 <- VlnPlot(
  seurat_integrated,
  features = "S.Score",
  group.by = "seurat_clusters",
  pt.size = 0
) +
  xlab("Unannotated Clusters") +
  labs(title = "S Phase Geneset Score") +
  NoLegend()

p2 <- VlnPlot(
  seurat_integrated,
  features = "G2M.Score",
  group.by = "seurat_clusters",
  pt.size = 0,
) +
  xlab("Unannotated Clusters") +
  labs(title = "G2M Phase Geneset Score") +
  NoLegend()

p1 + p2 & theme(axis.title.x = element_blank()) & NoLegend()

proliferating_clusters <- c("14", "17", "19", "20", "23")

# ==============================================================================
# 12. MAJOR GROUP ANNOTATION
# ==============================================================================

# General workflow:

# Identify “major cell groups” based on:
# Expression of established cell type marker genes (as presented above)
# Label transfer (Seurat) from previously annotated scRNA-Seq human airway datasets
# Carraro et al. scRNA-Seq HAE datasets (Carraro et al. 2021)

# For each major cell group:
# Subset and determine optimal cluster resolution
# Identify genes that distinguish subclusters and annotate based on literature and additional sources
# Where “biologically appropriate” (e.g. with limited biologically significant differences in gene expression), subclusters will be merged.
# Subclusters will be compared to reference classifications and assessed to address any discrepancies.
# Evaluate and annotate any remaining clusters that may not fit within major cell groups.

# ---- Label transfer from reference datasets ----------------------------------
# downloads and preprocesses list of seurat objects in `carraro_datasets`
source(here::here("R/qc_prep_hae_references.R"))

plan("sequential")

for (x in names(carraro_datasets)) {
  seurat_integrated <- label_transfer_seurat(
  seurat = seurat_integrated, reference = carraro_datasets[[x]],
  query_assay = "integrated", ref_assay = "integrated", norm_method = "SCT",
  labels = "major", name = x)
  }

label_transfers <- c(
  "carraro_ALI_predicted.id",
  "carraro_hBE_predicted.id"
  )

p_umap_labeltransfers_list <- lapply(label_transfers, function(x){
  p <- DimPlot(seurat_integrated, group.by = x, shuffle = TRUE) +
    theme(
      legend.position = "bottom",
    )
  return(p)
})

p_umap_labeltransfers_list[[
  length(p_umap_labeltransfers_list) + 1
  ]] <- p_umap_clus

wrap_plots(p_umap_labeltransfers_list, nrow = 1)

# ---- Assign major groups -----------------------------------------------------
seurat_integrated <- assign_major_group(
  seurat_integrated,
  reference_name = "carraro_ALI_predicted.id",
  clusters = "seurat_clusters")

# Assign intermediate major group
intermediate_clusters <- c("4", "10", "11")
Idents(seurat_integrated) <- "seurat_clusters"
intermediate_clusters_cells <- WhichCells(seurat_integrated,
                                          idents = intermediate_clusters)
seurat_integrated@meta.data[
  colnames(seurat_integrated) %in% intermediate_clusters_cells, "group_major"] <-
  "Intermediate"

# Remove "+" from FOXN4 for naming purposes
seurat_integrated$group_major <- gsub(seurat_integrated$group_major,
                                      pattern = "\\+", replacement = "")

p_group_major <- DimPlot(
  seurat_integrated,
  reduction = "umap",
  group.by = "group_major",
  label = TRUE,
  repel = TRUE
) + 
  scale_color_manual(values = cols25()) + 
  ggtitle("Major cell group assignment") +
  theme(
      legend.position = "bottom",
      axis.text = element_text(size = 10),
      title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) + coord_fixed()
p_group_major

# ---- Major group stability / subcluster resolution ---------------------------
p_tree <- clustree(seurat_integrated,
                   prefix = "integrated_snn_res.",
                   node_colour = "group_major",
                   node_colour_aggr = "as.character",
                   show_axis = TRUE) +
  scale_color_manual(values = cols25())
p_tree

ref_res <- 1.4 # For all major groups except ciliated
seurat_integrated$group_major_clusters <- paste(
  seurat_integrated$group_major,
  as.character(seurat_integrated[[paste0("integrated_snn_res.", ref_res)]][,]),
  sep = "_"
  )
ciliated_res <- 1
Idents(seurat_integrated) <- "group_major"
ciliated_cells <- WhichCells(seurat_integrated, idents = "Ciliated")
seurat_integrated[["group_major_clusters"]][ciliated_cells,] <- paste(
  "Ciliated",
  as.character(seurat_integrated[[paste0("integrated_snn_res.", ciliated_res)]][ciliated_cells,]),
  sep = "_"
  )

# ---- Ionocyte / tuft-like subclustering --------------------------------------
seurat_integrated_ionocyte <- subset(seurat_integrated,
                                     group_major == "Ionocyte")
# Run PCA/UMAP on subset data
DefaultAssay(seurat_integrated_ionocyte) <- "integrated"
seurat_integrated_ionocyte <- RunPCA(
  seurat_integrated_ionocyte, npcs = 100,
  verbose = FALSE
)

ndims <- 4 # Based on data exploration.
plan("multisession", workers = ncores)
seurat_integrated_ionocyte <- RunUMAP(seurat_integrated_ionocyte,
                          reduction = "pca",
                          dims = 1:ndims, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:ndims, verbose = FALSE)

seurat_integrated_ionocyte <- FindClusters(
  seurat_integrated_ionocyte,
  resolution = seq(from = 0.2, to = 3, by = 0.2),
  verbose = FALSE
)
clustree(seurat_integrated_ionocyte, prefix = "integrated_snn_res.")

ref_res <- 1.4 # Very stable
seurat_integrated_ionocyte <- 
  assign_cluster_resolution(seurat_integrated_ionocyte, res = ref_res)

p_umap_clus <- DimPlot(
  seurat_integrated_ionocyte,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + 
  NoLegend() +
  coord_fixed() +
  scale_color_manual(values = cols25())
p_umap_clus

seurat_integrated_ionocyte <- PrepSCTFindMarkers(seurat_integrated_ionocyte)
markers <- FindAllMarkers(seurat_integrated_ionocyte, assay = "SCT")

# Based on differential gene expression, subclusters are annotated as:
# - 0: Tuft-like
# - 1: Ionocyte
# - 2: Ionocyte
# - 3: Ionocyte (precursor?)
# - 4: Tuft-like
# - 5: Artifact (ciliated)
# - 6: Artifact (secretory)

Idents(seurat_integrated_ionocyte) <- "seurat_clusters"
tuftlike_cells <- WhichCells(seurat_integrated_ionocyte, idents = c("0", "4"))
artifact_cells <- WhichCells(seurat_integrated_ionocyte, idents = c("5", "6"))
seurat_integrated[["group_major"]][tuftlike_cells,] <- "Tuft_like"
seurat_integrated[["group_major"]][artifact_cells,] <- "artifact"
seurat_integrated[["group_major_clusters"]][colnames(seurat_integrated_ionocyte),] <- 
    seurat_integrated[["group_major"]][colnames(seurat_integrated_ionocyte),]

# Remove artifact cells from main Seurat object
Idents(seurat_integrated) <- "group_major"
seurat_integrated <- 
  subset(seurat_integrated, idents = "artifact", invert = TRUE)
seurat_integrated$group_major <- factor(seurat_integrated$group_major)
seurat_integrated$group_major_clusters <- 
  factor(seurat_integrated$group_major_clusters)

# ==============================================================================
# 13. CLUSTER-LEVEL ANNOTATION
# ==============================================================================

groups_to_test <- c("Basal", "Intermediate", "Secretory", "Ciliated")

plan("multisession", workers = ncores)
group_major_subcluster_markers <- lapply(groups_to_test, function(group){
  seurat <- subset(seurat_integrated, group_major == group)
  seurat <- PrepSCTFindMarkers(seurat)
  Idents(seurat) <- "group_major_clusters"
  subclusters_to_test <- unique(as.character(seurat$group_major_clusters))
  subcluster_markers <- lapply(subclusters_to_test, function(subcluster){
    markers <- FindConservedMarkers(seurat,
                         ident.1 = subcluster,
                         grouping.var = "library",
                         assay = "SCT",
                         min.cells.group = 1,
                         only.pos = FALSE
                         )
    return(markers)
    })
  names(subcluster_markers) <- subclusters_to_test
  return(subcluster_markers)
  })
names(group_major_subcluster_markers) <- groups_to_test

# ---- Gene set module scores --------------------------------------------------

# Stress gene set from O’Flanagan et al. 2019 
oflanagan_list <- read.csv(file = here::here("supporting-data/Campbell_StressCoreGeneSet.csv"))
oflanagan_stress_geneset <- list(
  as.character(oflanagan_list$gene_symbol[!is.na(oflanagan_list$gene_symbol)])
)
seurat_integrated <- AddModuleScore(
  seurat_integrated,
  features = oflanagan_stress_geneset,
  name = "stress_geneset_score",
  assay = "SCT"
)

# HALLMARK EMT gene set
HALLMARK_EMT_geneset <- read.table(
  file =  here::here("supporting-data/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  skip = 2)[,1]
HALLMARK_EMT_geneset <- list(as.character(HALLMARK_EMT_geneset))

seurat_integrated <- AddModuleScore(
  seurat_integrated,
  features = HALLMARK_EMT_geneset,
  name = "HALLMARK_EMT_geneset_score",
  assay = "SCT"
)

# ---- Basal major group -------------------------------------------------------

clusters <- names(group_major_subcluster_markers[["Basal"]])
output_list <- lapply(clusters, function(x){
  top_markers(
    list_markers = group_major_subcluster_markers[["Basal"]], cluster = x)
  })
names(output_list) <- clusters
print(output_list)

seurat_integrated <- PrepSCTFindMarkers(seurat_integrated)
# Pairwise comparisons of differentiating basal clusters
Idents(seurat_integrated) <- "group_major_clusters"
G2M_degs <- FindConservedMarkers(seurat_integrated,
                                ident.1 = "Basal_20", ident.2 = "Basal_23",
                                grouping.var = "library", assay = "SCT")
S_degs <- FindConservedMarkers(seurat_integrated,
                                ident.1 = "Basal_14", ident.2 = "Basal_17",
                                grouping.var = "library", assay = "SCT")

clusters_to_test <- c("Basal_2", "Basal_5", "Basal_8")

plan("multisession", workers = ncores)
Idents(seurat_integrated) <- "group_major_clusters"
seurat <- subset(seurat_integrated, group_major_clusters %in% clusters_to_test)
seurat <- PrepSCTFindMarkers(seurat)
Idents(seurat) <- "group_major_clusters"
subclusters_to_test <- unique(as.character(seurat$group_major_clusters))
subcluster_markers <- lapply(subclusters_to_test, function(subcluster){
    markers <- FindConservedMarkers(seurat,
                         ident.1 = subcluster,
                         grouping.var = "library",
                         assay = "SCT",
                         min.cells.group = 1,
                         only.pos = TRUE
                         )
    return(markers)
    })
  names(subcluster_markers) <- subclusters_to_test

Idents(seurat_integrated) <- "group_major_clusters"
p1 <- VlnPlot(
  seurat_integrated,
  features = "S.Score",
  group.by = "group_major_clusters",
  idents = clusters,
  pt.size = 0
) +
  labs(title = "S Phase Geneset Score") 

p2 <- VlnPlot(
  seurat_integrated,
  features = "G2M.Score",
  group.by = "group_major_clusters",
  idents = clusters,
  pt.size = 0,
) +
  labs(title = "G2M Phase Geneset Score")

p1 + p2 & vln_theme & NoLegend()

Idents(seurat_integrated) <- "group_major_clusters"
p1 <- VlnPlot(
  seurat_integrated,
  features = "HALLMARK_EMT_geneset_score1",
  group.by = "group_major_clusters",
  idents = clusters,
  pt.size = 0
) +
  labs(title = "HALLMARK EMT geneset score") 

p1 & vln_theme & NoLegend()

# ---- Intermediate major group ------------------------------------------------

clusters <- names(group_major_subcluster_markers[["Intermediate"]])
output_list <- lapply(clusters, function(x){
  top_markers(
    list_markers = group_major_subcluster_markers[["Intermediate"]], cluster = x)
  })
names(output_list) <- clusters
print(output_list)

# ---- Secretory major group ---------------------------------------------------

clusters <- names(group_major_subcluster_markers[["Secretory"]])
output_list <- lapply(clusters, function(x){
  top_markers(
    list_markers = group_major_subcluster_markers[["Secretory"]], cluster = x)
  })
names(output_list) <- clusters
print(output_list)

# Additional secretory markers
features_to_plot <- c("MSMB", "MUC5AC", "MUC5B",
                      "SCGB1A1", "CYP2F1", "ALDH3A1",
                      "SPDEF",  "BPIFA1", "BPIFB1",
                      "IL32", "TGFB2")
p_dot <- DotPlot(seurat_integrated, features = rev(features_to_plot),
                 assay = "SCT", idents = clusters, cluster.idents = TRUE, 
                 scale = FALSE) +
  scale_color_viridis() + 
  ggtitle("Select secretory cell genes") +
  coord_flip() +
  theme(
    plot.title = element_text(size = 10, face = "plain"),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    aspect.ratio = 1
  )
p_dot

Idents(seurat_integrated) <- "group_major_clusters"
p1 <- VlnPlot(
  seurat_integrated,
  features = "stress_geneset_score1",
  group.by = "group_major_clusters",
  idents = clusters,
  pt.size = 0
) +
  labs(title = "O'Flanagan stress geneset score") 

p1 & vln_theme & NoLegend()

# ---- Ciliated major group ----------------------------------------------------

clusters <- names(group_major_subcluster_markers[["Ciliated"]])
output_list <- lapply(clusters, function(x){
  top_markers(
    list_markers = group_major_subcluster_markers[["Ciliated"]], cluster = x)
  })
names(output_list) <- clusters
print(output_list)

Idents(seurat_integrated) <- "group_major_clusters"
p1 <- VlnPlot(
  seurat_integrated,
  features = "stress_geneset_score1",
  group.by = "group_major_clusters",
  idents = clusters,
  pt.size = 0
) +
  labs(title = "O'Flanagan stress geneset score")

p1 & vln_theme & theme(aspect.ratio = 1.5) & NoLegend()

# ---- Additional cell types ---------------------------------------------------
# FOXN4, Ionocyte, Tuft-like, and Neuroendocrine major groups do
# not contain subclusters for additional annotation.

# ==============================================================================
# 14. ASSIGN CLUSTERS TO MINOR GROUP ANNOTATIONS
# ==============================================================================

df_metadata <- seurat_integrated@meta.data
df_metadata <- df_metadata %>%
  mutate(
    group_clusters = case_when(
      # BASAL MAJOR GROUP
      group_major_clusters == "Basal_2" ~ "core",
      group_major_clusters == "Basal_5" ~ "core",
      group_major_clusters == "Basal_7" ~ "differentiating",
      group_major_clusters == "Basal_8" ~ "core",
      group_major_clusters == "Basal_14" ~ "proliferating(S)",
      group_major_clusters == "Basal_17" ~ "proliferating(S/G2M)",
      group_major_clusters == "Basal_20" ~ "proliferating(G2M)",
      group_major_clusters == "Basal_21" ~ "EMT",
      group_major_clusters == "Basal_23" ~ "differentiating proliferating(G2M)",
      # INTERMEDIATE MAJOR GROUP
      group_major_clusters == "Intermediate_4" ~ "secretory_like",
      group_major_clusters == "Intermediate_10" ~ "S100A8_hi",
      group_major_clusters == "Intermediate_11" ~ "basal_like",
      # SECRETORY MAJOR GROUP
      group_major_clusters == "Secretory_1" ~ "MSMB_hi MUC5B_hi MUC5AC_lo",
      group_major_clusters == "Secretory_3" ~ "MSMB_lo MUC5B_lo MUC5AC_lo",
      group_major_clusters == "Secretory_6" ~ "MSMB_hi MUC5B_hi MUC5AC_hi",
      group_major_clusters == "Secretory_9" ~ "MSMB_lo MUC5B_hi MUC5AC_lo",
      group_major_clusters == "Secretory_13" ~ "IGFBP3_hi",
      group_major_clusters == "Secretory_16" ~ "MSMB_hi MUC5B_lo MUC5AC_hi",
      group_major_clusters == "Secretory_18" ~ "stress signature",
      # CILIATED MAJOR GROUP
      group_major_clusters == "Ciliated_0" ~ "core",
      group_major_clusters == "Ciliated_18" ~ "stress signature",
      # FOXN4 MAJOR GROUP
      group_major == "FOXN4" ~ "core",
      # IONOCYTE MAJOR GROUP
      group_major == "Ionocyte" ~ "core",
      # TUFTLIKE MAJOR GROUP
      group_major == "Tuft_like" ~ "core",
      # NE MAJOR GROUP
      group_major == "NE" ~ "core"
    )
  )
df_metadata$cluster_annotation <- 
  paste(df_metadata$group_major, df_metadata$group_clusters, sep = " ") %>%
  gsub(pattern = "core", replacement = "") %>%
  trimws() %>%
  as.factor()

seurat_integrated <- AddMetaData(seurat_integrated,
                                  metadata = df_metadata
                                  )

# ---- Factor level ordering ---------------------------------------------------
# Set level ordering manually for consistency across figures/analyses

# Group major
seurat_integrated$group_major <- factor(seurat_integrated$group_major, levels =
                                        c("Basal",
                                        "Intermediate",
                                        "Secretory",
                                        "Ciliated",
                                        "FOXN4",
                                        "Ionocyte",
                                        "Tuft_like",
                                        "NE"))

# Timepoint
seurat_integrated$timepoint <- factor(seurat_integrated$timepoint,
                                      levels = c("0", "2", "6"))

# ---- Clusters to exclude -----------------------------------------------------

clusters_to_exclude <- c(
  "Basal EMT",
  "Ciliated stress signature",
  "Secretory IGFBP3_hi",
  "Secretory stress signature"
  )
seurat_integrated@meta.data <- seurat_integrated@meta.data %>%
  mutate(
    exclude_cluster = case_when(
      cluster_annotation %in% clusters_to_exclude ~ TRUE,
      TRUE ~ FALSE
      )
  )

# ==============================================================================
# 15. EXPORT
# ==============================================================================
# Remove intermediate cluster columns from metadata
columns_to_remove <-
  grep(
    x = colnames(seurat_integrated@meta.data),
    pattern = "^integrated_snn_res*|prediction.score",
    value = TRUE
    )
seurat_integrated@meta.data[, columns_to_remove] <- NULL

# Output RDS file
saveRDS(seurat_integrated, file = here::here("outputs/haeifn_seurat_integrated.rds"))

# ==============================================================================
# 16. REFERENCES
# ==============================================================================
# 1. Carraro, Gianni, Justin Langerman, Shan Sabri, Zareeb Lorenzana, Arunima Purkayastha, Guangzhu Zhang, Bindu Konda, et al. 2021. “Transcriptional Analysis of Cystic Fibrosis Airways at Single-Cell Resolution Reveals Altered Epithelial Cell States and Composition.” Nat. Med. 27 (5): 806–14.
# 2. McGinnis, Christopher S, David M Patterson, Juliane Winkler, Daniel N Conrad, Marco Y Hein, Vasudha Srivastava, Jennifer L Hu, et al. 2019. “MULTI-Seq: Sample Multiplexing for Single-Cell RNA Sequencing Using Lipid-Tagged Indices.” Nat. Methods 16 (7): 619–26.
# 3. O’Flanagan, Ciara H., Kieran R. Campbell, Allen W. Zhang, et al. 2019. “Dissociation of Solid Tumor Tissues with Cold Active Protease for Single-Cell RNA-Seq Minimizes Conserved Collagenase-Associated Stress Responses.” Genome Biology 20 (1): 210. https://doi.org/10.1186/s13059-019-1830-0.
