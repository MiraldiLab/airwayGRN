# Fig_4AB_IRG_identification.R
#
# Code base used to identify IFN-responsive genes (IRGs) in basal, suprabasal, ciliated, deuterosomal 
# and secretory cells (Fig. 4A-B)
#
# Bejjani et al. (2026) "Gene regulatory networks define human airway epithelial
# cell types and their distinct responses to type I interferon"
# Author: Brad Rosenberg, Icahn School of Medicine at Mount Sinai, NY. Modified by 
# Anthony Bejjani, Cincinnati Children's Hospital Medical Center, OH

# REQUIRED DATA: 
# Download 'GSE330155_scRNA_HAEifn_Seurat.rds' from GEO (GSE330155)
# Place the downloaded file in the "inputs" folder of this repository

# load necessary packages
library(Seurat)
library(tidyverse)
library(sva)
library(dplyr)
library(purrr)
library(edgeR)
library(here)

# helper functions
source(here('helper_differential_gene_expression.R'))

# set and create output directory
outdir <- here('outputs')
dir.create(outdir, recursive=T)

# Load Seurat scRNA-seq object (update directory to scRNA-seq object)
seurat_dge <- readRDS(here("inputs", "GSE330155_scRNA_HAEifn_Seurat.rds"))

# subset seurat object to exclude second run of donor B
series_to_exclude <- "BRRO04_HAEintB"
Idents(seurat_dge) <- "orig.ident"
seurat_dge <- subset(seurat_dge, 
                     idents = series_to_exclude, invert = TRUE)

# Require a minimum of 10 cells per population per condition (cluster level)
min_cells <- 10
mat <- table(seurat_dge$cluster_annotation,
             seurat_dge$donor_timepoint_series) %>% as.matrix()
             
# Less than 10 cells in at least 2 samples
clusters_to_exclude_count <- rowSums(mat < 10) > 2 
clusters_to_exclude_count <- which(
  clusters_to_exclude_count == TRUE) %>% names()
  
# Exclude cell types not relevant to present study
clusters_to_exclude_type <- levels(seurat_dge$cluster_annotation) %>%
  grep(pattern = "proliferating|EMT|stress|IGFBP3|Intermediate", value = TRUE)

# Subset seurat object
Idents(seurat_dge) <- "cluster_annotation"
seurat_dge <- subset(seurat_dge, 
                     idents = unique(c(clusters_to_exclude_count,
                                       clusters_to_exclude_type)),
                     invert = TRUE)

# Generate new column for testing groups ("dge_group")
seurat_dge$dge_group <- seurat_dge$cluster_annotation
seurat_dge$dge_group <- gsub(seurat_dge$dge_group,
                             pattern = "Basal differentiating",
                             replacement = "Suprabasal")
seurat_dge$dge_group <- gsub(seurat_dge$dge_group,
                             pattern = "FOXN4",
                             replacement = "Deuterosomal")

# Trim off specific cluster annotation to define dge group
seurat_dge$dge_group <- gsub(seurat_dge$dge_group,
                             pattern = " (.*)",
                             replacement = "")
seurat_dge$dge_group <- factor(seurat_dge$dge_group)

# Sum counts for pseudobulk dataset
summed <- prep_pseudobulk(seurat_dge,
                          ident = "dge_group",
                          populations = levels(seurat_dge$dge_group),
                          condition = "timepoint", block = "donor",
                          mincells = 0, min_prop = 0)

# Reorder object by celltype then by timepoint
summed <- summed[, order(summed$dge_group, summed$timepoint)]
sample_metadata <- colData(summed) %>%
  as.data.frame() %>%
  select_if(~ !any(is.na(.)))

## Compound design
# Build compound factors
# Build compound factors
sample_metadata$dge_group_timepoint <- paste(sample_metadata$dge_group,
                                             sample_metadata$timepoint,
                                             sep = "_")
cols <- which(colnames(sample_metadata) == "donor.1"):ncol(sample_metadata)

# Set design matrix
design <- model.matrix(~0 + dge_group_timepoint + donor, data = sample_metadata)
colnames(design) <- gsub(colnames(design),
                         pattern = "dge_group_timepoint",
                         replacement = "")

# Set factor vectors
## dge_groups/celltypes to include in clustering analysis
celltype <- unique(sample_metadata$group_major)
## BUILD CONTRASTS FOR EACH CELL TYPE (Timepoint vs Baseline)
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

## Add the 6 vs 2 contrast manually to the matrices
for(i in 1:6){
  contrasts_time[[i]] <- cbind(contrasts_time[[i]], as.data.frame(rep(0,21)))
  colnames(contrasts_time[[i]])[3] <- paste0(celltype[i], "_", '6', " - ", celltype[i], "_", "2")
  contrasts_time[[i]][paste0(celltype[i],'_2'),3] <- -1
  contrasts_time[[i]][paste0(celltype[i],'_6'),3] <- 1
}

## Combine timepoint contrasts for each cell type into one large matrix
contrasts_time <- do.call(cbind, contrasts_time)

## Format column names safely for downstream use
colnames(contrasts_time) <- paste0(
  "time_", gsub(colnames(contrasts_time), pattern = " - ", replacement = "_V_"))

# Set significance cutoffs
logFC_threshold <- 1
FDR_threshold <- 1

# Run the overall fit
fit <- run_standard_edgeR_fit(counts = counts(summed),
                              group = sample_metadata$dge_group_timepoint)

# EXTRACT AND MERGE DIFFERENTIAL EXPRESSION RESULTS
list_celltype_irgs <- lapply(celltype, function(x){
  
  # A. Identify the specific contrasts for this cell type (e.g., 2v0, 6v0, 6v2)
  celltype_contrast <- grep(colnames(contrasts_time),
                            pattern = paste0("^time_", x, "_[0-9]"),
                            value = TRUE)
  
  # B. Iterate through each contrast individually to get accurate FDRs/PValues
  res_list <- lapply(celltype_contrast, function(cont){
    
    # Run the test for the single contrast (drop=FALSE prevents vector coercion)
    res <- glmQLFTest(fit, contrast = contrasts_time[, cont, drop = FALSE])
    
    # Extract results and move rownames to a proper 'gene' column
    df <- as.data.frame(topTags(res, n = Inf)) 
    df <- cbind("gene" = rownames(df), df)
    
    # Rename statistical columns dynamically to match the contrast
    colnames(df)[colnames(df) == "logFC"] <- paste0("logFC_", cont)
    colnames(df)[colnames(df) == "PValue"] <- paste0("PValue_", cont)
    colnames(df)[colnames(df) == "FDR"] <- paste0("FDR_", cont)
    colnames(df)[colnames(df) == "F"] <- paste0("F_", cont)
    
    return(df)
  })
  
  # C. Merge the separate contrast results into one wide data frame
  # Joining by both 'gene' and 'logCPM' prevents duplicate columns
  # Explicitly calling purrr::reduce avoids Bioconductor namespace clashes
  combined_df <- purrr::reduce(res_list, full_join, by = c("gene", "logCPM"))
  
  return(combined_df)
})

# Name the final outer list by the cell types
names(list_celltype_irgs) <- celltype
saveRDS(list_celltype_irgs, file.path(outdir,'celltype_irgs.rds'))

# Combine all matrices into 1 matrix
combined_all <- imap_dfr(list_celltype_irgs, function(df, ct_name) {
  
  # Remove the specific cell type from the column headers
  # e.g., "logFC_time_Basal_2_V_0" becomes "logFC_time_2_V_0"
  df_clean <- df %>%
    rename_with(~ gsub(paste0("_", ct_name), "", .x))
  
  return(df_clean)
  
}, .id = "celltype")
write.table(combined_all, file.path(outdir, 'celltype_IRGs_combined.tsv'), sep='\t', col.names = T, row.names = F, quote=F)

