# Celltype_IRP_Identification.R
#
# Code base used to identify IFN-responsive peaks (IRPs) in basal, suprabasal, ciliated 
# and secretory cells (Fig. 5A)
#
# Bejjani et al. (2026) "Gene regulatory networks define human airway epithelial
# cell types and their distinct responses to type I interferon"
# Author: Brad Rosenberg, Icahn School of Medicine at Mount Sinai, NY. Modified by 
# Anthony Bejjani, Cincinnati Children's Hospital Medical Center, OH

# REQUIRED DATA: 
# Download 'GSE330156_snATAC_HAEifn_Seurat.rds' from GEO (GSE330156)
# Place the downloaded file in the "inputs" folder of this repository.

# load necessary packages
library(Seurat)
library(Signac)
library(pals)
library(tidyverse)
library(sva)
library(here)

# helper functions
source(here('helper_differential_gene_expression.R'))

# set output directory
outdir <- here('outputs')
dir.create(outdir, recursive=T)

# load seurat object
seurat_dge <- readRDS(here('inputs','GSE330156_snATAC_HAEifn_Seurat.rds'))

# Require a minimum of 10 cells per population per condition (cluster level)
min_cells <- 10
seurat_dge$donor_timepoint <- paste(seurat_dge$Patient, seurat_dge$Time, sep='_')
mat <- table(seurat_dge$minor_group,
             seurat_dge$donor_timepoint) %>% as.matrix()

# Less than 10 cells in at least 2 samples
clusters_to_exclude_count <- rowSums(mat < 10) > 2 
clusters_to_exclude_count <- which(
  clusters_to_exclude_count == TRUE) %>% names()
clusters_to_exclude_count <- setdiff(clusters_to_exclude_count)

# Exclude cell types not relevant to present study
clusters_to_exclude_type <- unique(seurat_dge$minor_group) %>%
  grep(pattern = "proliferating|EMT|stress|IGFBP3|Intermediate", value = TRUE)

# Subset seurat object
Idents(seurat_dge) <- "minor_group"
seurat_dge <- subset(seurat_dge, 
                     idents = unique(c(clusters_to_exclude_count,
                                       clusters_to_exclude_type)),
                     invert = TRUE)

# Generate new column for testing groups ("dge_group")
seurat_dge$dge_group <- seurat_dge$minor_group
seurat_dge$dge_group <- gsub(seurat_dge$dge_group,
                             pattern = "Basal differentiating",
                             replacement = "Basal_differentiating")

# Trim off specific cluster annotation to define dge group
seurat_dge$dge_group <- gsub(seurat_dge$dge_group,
                             pattern = " (.*)",
                             replacement = "")
seurat_dge$dge_group <- factor(seurat_dge$dge_group)

seurat_dge <- RunUMAP(object = seurat_dge, reduction = 'integrated_lsi', dims = 2:30)
seurat_integrated <- RunUMAP(object = seurat_integrated, reduction = 'integrated_lsi', dims = 2:30)

p_umap_prefilter <- DimPlot(seurat_integrated,
                            group.by = "major_group",
                            label = TRUE, repel = TRUE) +
  labs(title = "HAE major groups pre-filter")


p_umap_postfilter <- DimPlot(seurat_dge,
                             group.by = "major_group",
                             label = TRUE, repel = TRUE) +
  labs(title = "HAE groups post-filter for DGE analysis")

(p_umap_prefilter + NoLegend()) + p_umap_postfilter

Idents(seurat_dge) <- 'major_group'
seurat_dge <- subset(seurat_dge, idents=c('Intermediate','NE'), invert=T)

# Set parallel settings to avoid memory crash
plan("sequential")
seurat_dge$major_group <- factor(seurat_dge$major_group)

# Sum counts for pseudobulk dataset
summed <- read.delim(here('inputs','scatac_IFN_bulk_major_group_raw.txt'),sep='\t',header=T, row.names = 1)
summed <- summed[,!grepl('NE',colnames(summed))]
summed <- summed[,!grepl('Intermediate',colnames(summed))]
curr_meta <- read.delim(here('inputs','scatac_IFN_metabulk.txt'),sep='\t',header=T, row.names = 1)
curr_meta <- curr_meta[colnames(summed),]
curr_meta$cell_time <- paste(curr_meta$major_group, curr_meta$Time, sep='_')
summed <- SingleCellExperiment(assays = list(counts = summed), colData=curr_meta)

# Reorder object by celltype then by timepoint
summed <- summed[, order(summed$major_group, summed$Time)]
sample_metadata <- colData(summed) %>%
  as.data.frame() %>%
  select_if(~ !any(is.na(.)))

## Compound design
# Build compound factors
sample_metadata$dge_group_timepoint <- paste(sample_metadata$major_group,
                                             sample_metadata$Time,
                                             sep = "_")
cols <- which(colnames(sample_metadata) == "donor.1"):ncol(sample_metadata)

# Set design matrix
design <- model.matrix(~0 + dge_group_timepoint + Patient, data = sample_metadata)
colnames(design) <- gsub(colnames(design),
                         pattern = "dge_group_timepoint",
                         replacement = "")

# Set factor vectors
## dge_groups/celltypes to include in DGE analysis
celltype <- unique(sample_metadata$major_group)

## timepoints
timepoint <- "0"

## For determining celltype markers, contrast
## each celltype with the AVERAGE of all other cell types
## Celltype 'marker genes' are the baseline contrast (t = 0)
contrasts_celltype_markers <- lapply(celltype, function(test){
  other <- celltype[celltype != test] # get cell types other than self, x
  # Build strings for each contrast and run through makecontrasts()
  contrast_list <- lapply(timepoint, function(x){
    contrast_string <- paste0(
      test, "_", x, 
      " - (", paste(paste0(other, "_", x), collapse = " + "),")/",
      length(other))
    cmd <- paste(
      "contrast <- makeContrasts(", contrast_string, ", levels = design)", 
      sep = '"'
    )
    eval(parse(text = cmd))
  })
  # Combine contrasts into dataframe
  contrast_df <- do.call(cbind, contrast_list)
  return(contrast_df)
})
# Combine all celltype contrasts
mat_contrasts_celltype_markers <- do.call(cbind, contrasts_celltype_markers)

# Set column names
colnames(mat_contrasts_celltype_markers) <- paste0("celltype_",
                                                   gsub(colnames(mat_contrasts_celltype_markers),
                                                        pattern = " - ", replacement = "_V_")
)

# Set significance cutoffs
logFC_threshold <- 1
FDR_threshold <- 0.05

# Run edgeR fit
fit <- run_standard_edgeR_fit(counts = counts(summed),
                              group = sample_metadata$dge_group_timepoint)

# Pairwise tests for IRPs within each cell type
contrasts_celltype_markers <- colnames(mat_contrasts_celltype_markers)

res_list_celltype_markers <- lapply(
  contrasts_celltype_markers, function(x){
    res <- glmQLFTest(fit, contrast = mat_contrasts_celltype_markers[,x])
    df_res <- topTags(res, n = Inf) %>% as.data.frame() %>%
      # Apply significance filters
      filter(FDR < FDR_threshold) %>%
      filter(logFC > logFC_threshold) # Only positive logFC for markers
    df_res <- cbind("peak" = rownames(df_res), df_res)
    return(df_res)
  })
names(res_list_celltype_markers) <- celltype

saveRDS(res_list_celltype_markers, file.path(outdir,'celltype_markers.rds'))

# Set factor vectors
## dge_groups/celltypes to include in clustering analysis
celltype <- unique(sample_metadata$major_group)
## timepoints
timepoint <- unique(sample_metadata$Time)

## For determining IRPs, within a given cell type, IFN timepoint vs baseline
contrasts_time <- lapply(celltype, function(celltype){
  timepoint <- timepoint[timepoint != "0"] # Remove baseline
  contrast_list <- lapply(timepoint, function(x){
    # Build strings for each contrast and run through makecontrasts()
    contrast_string <- paste0(celltype, "_", x, " - ", celltype, "_", "0")
    cmd <- paste(
      "contrast <- makeContrasts(", contrast_string, ", levels = design)",
      sep = '"')
    eval(parse(text = cmd))
  })
  contrast_df <- do.call(cbind, contrast_list)
})
# Combine timepoint contrasts for each cell type
contrasts_time <- do.call(cbind, contrasts_time)
# Set column names
colnames(contrasts_time) <- paste0(
  "time_", gsub(colnames(contrasts_time), pattern = " - ", replacement = "_V_"))

# Set significance cutoffs
logFC_threshold <- 1
FDR_threshold <- 0.05

# Run standard edgeR fit
fit <- run_standard_edgeR_fit(counts = counts(summed),
                              group = sample_metadata$cell_time)

list_celltype_IRPs <- lapply(celltype, function(x){
  celltype_contrast <- grep(colnames(contrasts_time),
                            pattern = paste0("^time_", x, "_[0-9]"),
                            value = TRUE)
  res <- glmQLFTest(fit, contrast = contrasts_time[,celltype_contrast])
  df <- as.data.frame(topTags(res, n = Inf)) %>%
    filter(FDR < FDR_threshold) %>%
    filter(if_any(starts_with("logFC"), ~ abs(.) > logFC_threshold))
  df <- cbind("gene" = rownames(df), df)
  return(df)
})
names(list_celltype_IRPs) <- celltype
saveRDS(list_celltype_IRPs, file.path(outdir,'celltype_irps.rds'))



