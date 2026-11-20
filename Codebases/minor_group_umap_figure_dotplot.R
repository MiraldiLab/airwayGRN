# ==============================================================================
# Minor group UMAP and dot plots
# Alec P. Pankow & Brad R. Rosenberg
# Icahn School of Medicine at Mount Sinai, New York, NY USA
# ==============================================================================

# Packages
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pals)
library(Seurat)
library(forcats)

# point to input .rds file here
seurat_integrated <- readRDS(here("outputs/haeifn_seurat_integrated.rds"))

# point to desired output name (infers format from suffix)
umap_output_fh <- here("outputs/minor_group_umap.pdf")
dotplot_output_fh <- here("outputs/minor_group_dotplot.pdf")

### UMAP helper functions ######################################################

#' Custom DimPlot wrapper to generate figure-ready plots. UMAPs with circle labels.
#'
#' @param object A Seurat object
#' @param group.ident A grouping variable to color points by. If colored by cluster identity, object must be processed by FindClusters.
#' @param colors A vector of colors
#' @param pt.size Size of points on UMAP
#' @param legend.col Number of columns of legend
#' @param label Logical to add labels
#' @param ciclelabels Logical to add circle labels to cluster labels. Works effectively for distinct clusters on UMAP space.
#' @param byrow Legend specification. List variables by row (TRUE) or by column (FALSE).
#' @param inset Logical to specify placement of legend. Inset (TRUE) or out by the side.
#' @param legend.position Placement of the plot. Inherits from inset parameter. Default: top_left
#' @param shape Shape of plotted points
#' @param legend.title Title of legend
#' @param legend.text.align Alignment of legend text
#'
#' @return A publication ready UMAP visualization, where points are colored by group.identity.
#'
#' @examples
#' library(Seurat)
#' data("pbmc_small")
#' CustomUMAP(object = pbmc_small, group.ident = "RNA_snn_res.1", inset = F)
#'
#' @export CustomUMAP

# Alternative, more refined function to add circle labels to UMAP visualizations
CustomUMAP <- function(object,
                       group.ident = "integrated_snn_res.1",
                       colors = legocolors$hex[-1],
                       pt.size = 0.01,
                       legend.col = 4,
                       plot.title = "Unsupervised clustering",
                       ciclelabels = TRUE,
                       byrow = TRUE,
                       label = FALSE,
                       #inset = TRUE,
                       shape = 21,
                       legend.title = "Cluster",
                       legend.position = "top left",
                       legend.text.align = 0,
                       pad = 0.1) {
  require(Seurat)
  require(ggplot2)
  #library(legocolors)
  
  cellnumber <- dim(object)[2]
  group.levels <- levels(object@meta.data[, group.ident])
  p2 <-
    DimPlot(
      object = object,
      pt.size = pt.size,
      label = label,
      group.by = group.ident,
    ) +
    scale_color_manual(values = colors, labels = paste0(1:length(group.levels), ". ", group.levels)) +
    scale_fill_manual(values = colors) +
    xlab(expression("UMAP"[1])) +
    ylab(expression("UMAP"[2])) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_text(hjust = 0, vjust = 0),
      axis.title.y = element_text(hjust = 0, vjust = 0),
      legend.title = element_text(face = "bold"),
      legend.text.align = legend.text.align,
      legend.text = element_text(
        color = "black", inherit.blank = T,
        margin = margin(t = 0, r = 0, b = 0, l = -5)
      ),
      legend.key.width = unit(0, "mm"),
      legend.key = element_rect(fill = NA, inherit.blank = T),
      legend.background = element_rect(fill = NA, color = "black", size = 0.25),
      legend.margin = margin(0, 0, 0, 0),
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA)
    ) +
    labs(color = legend.title) +
    ggtitle(paste(plot.title, " (", scales::comma_format(accuracy = 12)(cellnumber), " cells)", sep = "")) +
    guides(colour = guide_legend(
      override.aes = list(
        shape = 21,
        size = 6,
        color = "black",
        fill = colors[1:length(group.levels)]
      ),
      label = T,
      byrow = byrow,
      ncol = legend.col
    ))
  if (ciclelabels) {
    p2 <- custom.LabelClusters(
      plot = p2,
      id = group.ident,
      clusters = levels(object@meta.data[, group.ident]),
      circle.size = 7,
      text.color = "black",
      text.size = 4,
      shape = shape,
      fill = colors[1:length(levels(object@meta.data[, group.ident]))],
      repel = T,
      pad = pad
    )
  }
  return(p2)
}

# Modify Seurat functions to plot geom_point (shape) and geom_text for label clusters, underlayed with a circles.
# This is key, as it identifies where the labels are plotted on the UMAP
GetXYAesthetics <- function(plot, geom = "GeomPoint", plot.first = TRUE) {
  geoms <- sapply(
    X = plot$layers,
    FUN = function(layer) {
      return(class(x = layer$geom)[1])
    }
  )
  geoms <- which(x = geoms == geom)
  if (length(x = geoms) == 0) {
    stop("Cannot find a geom of class ", geom)
  }
  geoms <- min(geoms)
  # mappings (used by newer Seurat DimPlot) to the plain column name
  if (plot.first) {
    x <- rlang::as_label(plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)
    y <- rlang::as_label(plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)
  } else {
    x <- rlang::as_label(plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)
    y <- rlang::as_label(plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)
  }
  return(list("x" = x, "y" = y))
}

# Custom function to add circles labels to UMAP visualizations
custom.LabelClusters <- function(
  plot, # Use DimPlot to generate base ggplot to apply function
  id, # The seurat cluster identifier
  clusters = NULL,
  labels = NULL,
  split.by = NULL,
  repel = F,
  pad = 0.1,
  colors = colors,
  circle.size = circle.size,
  text.size = text.size,
  ...) {
  xynames <- unlist(x = GetXYAesthetics(plot = plot), use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[, id])))
  groups <- levels(data[, id]) # clusters %||% as.character(x = na.omit(object = unique(x = data[, id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in% possible.clusters], collapse = ","))
  }
  labels.loc <- lapply(
    X = groups,
    FUN = function(group) {
      data.use <- data[data[, id] == group, , drop = FALSE]
      data.medians <- if (!is.null(x = split.by)) {
        do.call(
          what = "rbind",
          args = lapply(
            X = unique(x = data.use[, split.by]),
            FUN = function(split) {
              medians <- apply(
                X = data.use[data.use[, split.by] == split, xynames, drop = FALSE],
                MARGIN = 2,
                FUN = median,
                na.rm = TRUE
              )
              medians <- as.data.frame(x = t(x = medians))
              medians[, split.by] <- split
              return(medians)
            }
          )
        )
      } else {
        as.data.frame(x = t(x = apply(
          X = data.use[, xynames, drop = FALSE],
          MARGIN = 2,
          FUN = function(x){ median(x, na.rm = TRUE) + pad}
        )))
      }
      data.medians[, id] <- factor(group, levels = groups)
      return(data.medians)
    }
  )
  labels.loc <- do.call(what = "rbind", args = labels.loc)
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels), ") must be equal to the number of clusters being labeled (", length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  geom.use <- ifelse(test = repel, yes = geom_point, no = geom_label)
  plot <- plot + geom.use(
    data = labels.loc, size = circle.size, shape = 21, stroke = 0.66, col = "gray17",
    mapping = aes_string(x = xynames["x"], y = xynames["y"]),
    ...
  ) + geom_text(
    size = text.size,
    data = labels.loc, col = "gray17",
    mapping = aes_string(x = xynames["x"], y = xynames["y"], label = labels.loc[[id]] %>% as.numeric()),
    ...
  )
  return(plot)
}

###### UMAP color palettes #####################################################

major_cell_type_cols <- c(
  'Basal'="#edae49",
  'Suprabasal'='#66a182',
  'Ciliated'='#d1495b',
  'Secretory'='#00798c',
  'Deuterosomal'='#8d96a3',
  'NE'='#F76A2E',
  'Tuft_like'='#9A9624',
  'Intermediate'='#694698'
)

# be careful about changing order without changing factor... 
manual_cols_23 <- c(
  "Basal" = "#edae49", #"#ffe5a3",
  "Basal differentiating" = "#3a7324", #"#edae49",
  "Basal differentiating proliferating(G2M)" = "#d9ae79", #"#f4bf6b",
  "Basal EMT" = "grey60", 
  "Basal proliferating(G2M)" = "#edae49",
  "Basal proliferating(S)" = "#d18128",
  "Basal proliferating(S/G2M)" = "#b05509",
  
  "Ciliated" = "#d1495b",
  "Ciliated stress signature" = "grey60",
  
  "Deuterosomal" = "#23b852", #"#8d96a3",
  
  "Intermediate basal-like" = "#694698",
  "Intermediate S100A8^hi" = "#9a8fdb",
  "Intermediate secretory-like" = "#ab6ac9",
  
  "Ionocyte" = "#476994", #"#2e4057",
  "NE" = "#F76A2E", 
  
  "Secretory IGFBP3_hi" = "grey60", #"#1891a1",
  "Secretory MSMB_hi MUC5B_hi MUC5AC_hi" = "#b3e5e8",
  "Secretory MSMB_hi MUC5B_hi MUC5AC_lo" = "#80d6dd",
  "Secretory MSMB_hi MUC5B_lo MUC5AC_hi" = "#4dc7d1",
  "Secretory MSMB_lo MUC5B_hi MUC5AC_lo" = "#1ab8c5",
  "Secretory MSMB_lo MUC5B_lo MUC5AC_lo" = "#0099a8",
  "Secretory stress signature" = "grey60",
  
  "Tuft-like" = "#9A9624"
)

cb_palette_23 <- c(
  "#D55E00", "#CC79A7", "#999999", "#66C2A5", "#FC8D62",
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
  "#66A61E", "#40E0D0", "#E6AB02", "#40E0D0", "#A6761D"
)


#### Make UMAP plot  ###########################################################

meta <- FetchData(seurat_integrated, vars = c("cluster_annotation")) %>%
  mutate(
    cluster_pretty = fct_recode(
      cluster_annotation,
      "Intermediate basal-like" = "Intermediate basal_like",
      "Intermediate S100A8^hi" = "Intermediate S100A8_hi",
      "Intermediate secretory-like" = "Intermediate secretory_like",
      "Tuft-like" = "Tuft_like",
      "Deuterosomal" = "FOXN4"
    )
  )

seurat_integrated$cluster_pretty <- meta$cluster_pretty

p2 <- CustomUMAP(
  seurat_integrated,
  group.ident = "cluster_pretty",
  colors = manual_cols_23,
  pt.size = 0.1,
  legend.col = 2,
  plot.title = "Minor HAE group annotations",
  ciclelabels = TRUE,
  byrow = FALSE,
  label = FALSE,
  #inset = TRUE,
  shape = 21,
  legend.title = "Cell annotation",
  legend.position = "top left",
  legend.text.align = 0,
  pad = 0.7
) & coord_fixed()

ggsave(umap_output_fh, plot = p2, height = 6, width = 20)

#### Make dot plot #############################################################

# helper functions
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Scalable ggplot theme for dot plot
dot_theme_fxn <- function(base_size = 14) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      # Figure title
      plot.title = element_text(
        size = rel(1), face = "bold",
        margin = margin(0, 0, 5, 0), hjust = 0
      ),
      # Plotting area
      panel.border = element_rect(size = 0.25, fill = "transparent"),
      panel.grid = element_blank(),
      # Axis labels
      axis.title = element_blank(),
      axis.text = element_text(size = rel(0.60), face = "plain"),
      axis.text.x = element_text(angle = 60, hjust = 0, vjust = 0),
      axis.line = element_blank(), # axis lines from panel border
      # axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(size = 0.25),
      # Legend
      legend.title = element_text(size = rel(0.66), face = "plain"),
      legend.text = element_text(size = rel(0.66), face = "plain"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(5, "mm"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin()
    )
}

cluster_long_df <- function(df, id_col, names_col, values_col) {
  df.wide <- df %>%
    select(all_of(c(id_col, names_col, values_col))) %>%
    pivot_wider(id_cols = id_col, names_from = names_col, values_from = values_col)
  return(df.wide)
}

preprocess_seurat <- function(seurat, features, group_vars = "group_major") {
  DefaultAssay(seurat) <- "SCT"
  
  exp.dat <- FetchData(seurat, vars = c(group_vars, features)) %>%
    pivot_longer(any_of(features), names_to = "feature", values_to = "value")
  
  exp.dat %>%
    group_by(across(all_of(group_vars)), feature) %>%
    summarize(
      percent.expressed = sum(value > 0) / n() * 100,
      mean.expression = mean(value)
    ) %>% group_by(feature) %>%
    mutate(
      scaled.expression =  (mean.expression - mean(mean.expression)) / sd(mean.expression)
    )
}

generate_dotplot <- function(plot.data, group_var, limits = c(-2, 2)) {
  ggplot(
    plot.data,
    aes_string(group_var, "feature", color = "scaled.expression", size = "percent.expressed")
  ) +
    geom_point() +
    geom_tile(aes(fill = group_major, y = length(sel.features) + 1), size = 1, stroke = NA, color = NA) +
    dot_theme_fxn() +
    scale_color_distiller(type = "div", oob = scales::oob_squish, limits = limits, guide = guide_colorbar())
}
  
#### Helper functions end ######################################################

# in this object, data have already been transformed...
DefaultAssay(seurat_integrated) <- "SCT"
  
# select markers to plot
sel.features <- c(
  
  # basal markers
  "KRT5",
  "TP63",
  "KRT13",
  #"KRT14",
  "KRT17",
  
  # proliferation
  "PCNA",
  "CDK1",
  "TOP2A",
  
  # ciliated markers
  "FOXJ1",
  "CAPS",
  
  # deuterosomal markers
  "DEUP1",
  "FOXN4",
  
  # secretory markers
  "SPDEF",
  "SCGB1A1",
  "SCGB3A1",
  "MSMB",
  "MUC5AC",
  "MUC5B",
  
  # ionocyte markers
  "FOXI1",
  "ASCL3",
  "CFTR",
  
  # NE markers
  "ASCL1",
  "CHGA",
  
  # tuft-like marker
  "POU2F3",
  "ASCL2",
  
  # miscellaneous #
  
  # EMT
  "VIM",
  "COL4A2",
  
  # 'stress' 
  "IL32",
  "ATF3",
  #"stress_geneset_score1",
  "TGFB2",
  "S100A8",
  "IGFBP3"
)

# create and add recoded cluster column for display
meta <- FetchData(seurat_integrated, vars = c("cluster_annotation")) %>%
  mutate(
    cluster_pretty = fct_recode(
      cluster_annotation,
      "Intermediate basal-like" = "Intermediate basal_like",
      "Intermediate S100A8^hi" = "Intermediate S100A8_hi",
      "Intermediate secretory-like" = "Intermediate secretory_like",
      "Tuft-like" = "Tuft_like",
      "Deuterosomal" = "FOXN4"
    )
  )
seurat_integrated[["cluster_pretty"]] <- meta$cluster_pretty

# summarize data for dot plot
plot.data <- preprocess_seurat(seurat = seurat_integrated, features = sel.features, group_vars = c("cluster_pretty", "donor")) %>%
  mutate(
    cluster_pretty = factor(cluster_pretty, levels = names(manual_cols_23))
  )
base.size <- 10

# get unique facet levels in the order they appear in the data
facet_levels <- unique(plot.data[["cluster_pretty"]])
numbers <- seq_along(facet_levels)

# mapping: original level -> number (as character so it renders as strip text)
num_map <- setNames(as.character(numbers), facet_levels)

# strip shows only the number:
numeric_labeller <- as_labeller(num_map)

plot.data <- plot.data %>%
  mutate(
    dge_group = paste(donor, cluster_pretty, sep = "|"),
    group_major = factor(cluster_pretty)
    #donor = str_split(donor_timepoint_series, "_") %>% map_chr(1),
    #series = str_split(donor_timepoint_series, "_") %>% map_chr(3)
  ) %>%
  mutate(feature = factor(feature, levels = rev(sel.features)))

# call the plot
p <- generate_dotplot(plot.data, "donor") + 
  dot_theme_fxn(base_size = base.size) + 
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = manual_cols_23) +
  scale_size_area(max_size = 4) +
  facet_wrap(~ group_major, nrow = 1, scales = "free_x", labeller = numeric_labeller) +
  theme(
    panel.spacing = unit(0, "mm"),
    axis.text.x = element_blank(), #element_text(size = base.size - 2),
    axis.title.x.top = element_text(size = base.size),
    axis.ticks.x = element_blank(), #element_text(),
    axis.text.y = element_text(size = base.size, face = "oblique"),
    panel.border = element_rect(color = "grey70", linetype = 1, size = 0.5),
    strip.background = element_blank(),
    strip.text.x = element_text(size = base.size),
    legend.text = element_text(size = base.size),
    legend.title = element_text(size = base.size)
  ) + 
  guides(
    color = guide_colorbar(title = "Scaled expression"),
    fill = "none", #guide_legend(title = "Cell population"), 
    size = guide_legend(title = "Percent expressed")
  ) + 
  labs(x = "Cell annotation (by donor, all timepoints)")

# save figure
ggsave(dotplot_output_fh, plot = p, height = 8, width = 10.5)
