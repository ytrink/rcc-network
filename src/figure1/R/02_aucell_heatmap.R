# figures/figure1/R/01_aucell_heatmap.R
#
# Generates AUCell regulon activity heatmap ordered by broad cell type (Figure 1).
#
# Input:  data/li/Li_integrated.rds
#         network_analysis/outputs/auc_matrix_standardized.csv
#         data/leiden_RB_03_pos.csv
# Output: figures/figure1/R/heatmap_ordered_by_celltype.pdf

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(dplyr)
library(here)

source(file.path(here::here(), "utils", "complex_heatmap_functions.R"))

set.seed(101)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir    <- here::here()
data_dir    <- file.path(root_dir, "data")
network_dir <- file.path(root_dir, "network_analysis", "outputs")
out_dir     <- file.path(root_dir, "figures", "figure1", "R")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Parameters ─────────────────────────────────────────────────────────────────
N_CELLS_PER_TYPE <- 5000
MIN_REGULON_SIZE <- 5

CELL_TYPE_ORDER <- c(
  "B_cell", "T_NK", "Endothelial",
  "Fibroblast", "Macrophage_Monocyte", "Tumor"
)

CELL_TYPE_COLORS <- c(
  "B_cell"              = "darkorange",
  "T_NK"                = "dodgerblue",
  "Endothelial"         = "forestgreen",
  "Fibroblast"          = "purple",
  "Macrophage_Monocyte" = "deeppink",
  "Tumor"               = "black"
)

COMMUNITY_COLORS <- c(
  "Community_0" = "#1f77b4",
  "Community_1" = "#ff7f0e",
  "Community_2" = "#9467bd",
  "Community_3" = "#2ca02c",
  "Community_4" = "#d62728",
  "Community_5" = "#8c564b",
  "Community_6" = "#17becf"
)

# ── Load Seurat object and downsample ──────────────────────────────────────────
rna.li <- readRDS(file.path(data_dir, "li", "Li_integrated.rds"))
colnames(rna.li) <- paste0("LI_", rna.li$patient, "_", colnames(rna.li))

# Broad cell type grouping
rna.li$annotation_broad_1 <- case_when(
  rna.li$annotation_broad %in% c("CD4_T", "CD8_T", "NK") ~ "T_NK",
  rna.li$annotation_broad %in% c("Macrophage", "Monocyte") ~ "Macrophage_Monocyte",
  TRUE ~ rna.li$annotation_broad
)

# Downsample to N_CELLS_PER_TYPE per broad cell type
cells_keep <- rna.li@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  group_by(annotation_broad_1) %>%
  slice_sample(n = N_CELLS_PER_TYPE) %>%
  pull(barcode)

rna.li <- rna.li[, cells_keep]

# ── Load AUC matrix ────────────────────────────────────────────────────────────
heatmap_matrix <- fread(file.path(network_dir, "auc_matrix_standardized.csv")) %>%
  as.data.frame()
rownames(heatmap_matrix) <- heatmap_matrix[[1]]
heatmap_matrix[[1]]      <- NULL

# ── Filter regulons ────────────────────────────────────────────────────────────
regulon_sizes        <- as.numeric(sub(".*\\((\\d+)g\\)", "\\1", rownames(heatmap_matrix)))
heatmap_matrix       <- heatmap_matrix[regulon_sizes >= MIN_REGULON_SIZE, ]
heatmap_matrix       <- heatmap_matrix[, colnames(rna.li), drop = FALSE]

# ── Build annotation column ────────────────────────────────────────────────────
annotationCol <- data.frame(
  annotation_broad_1 = factor(rna.li$annotation_broad_1, levels = CELL_TYPE_ORDER),
  row.names          = colnames(rna.li)
)

stopifnot(all(rownames(annotationCol) == colnames(heatmap_matrix)))

# ── Order cells by broad cell type ────────────────────────────────────────────
ord_cols       <- order(annotationCol$annotation_broad_1)
hm_ordered     <- as.matrix(heatmap_matrix[, ord_cols])
ann_ordered    <- annotationCol[ord_cols, , drop = FALSE]

# ── Leiden community row annotation ───────────────────────────────────────────
leiden_clusters  <- read.csv(file.path(data_dir, "leiden_RB_03_pos.csv"))
all_tfs          <- as.vector(as.matrix(leiden_clusters))
gene_names       <- sub(" \\(.*\\)", "", rownames(heatmap_matrix))

cluster_assignment <- sapply(gene_names, function(g) {
  if (g %in% all_tfs) {
    names(leiden_clusters)[sapply(leiden_clusters, function(col) g %in% col)][1]
  } else {
    "None"
  }
})

cluster_df <- data.frame(
  leiden_cluster = factor(cluster_assignment, levels = names(COMMUNITY_COLORS)),
  row.names      = rownames(heatmap_matrix)
)

row_ha <- rowAnnotation(
  leiden_cluster        = cluster_df$leiden_cluster,
  col                   = list(leiden_cluster = COMMUNITY_COLORS),
  show_annotation_name  = FALSE
)

# ── Top annotation ─────────────────────────────────────────────────────────────
top_annotation <- HeatmapAnnotation(
  df  = ann_ordered,
  col = list(annotation_broad_1 = CELL_TYPE_COLORS)
)

# ── Plot ───────────────────────────────────────────────────────────────────────
draw_complex_heatmap.pdf.nonclustered_cells(
  heatmap_matrix    = hm_ordered,
  filename          = file.path(out_dir, "heatmap_ordered_by_celltype.pdf"),
  top_annotation    = top_annotation,
  right_annotation  = row_ha,
  distance.method   = "pearson",
  clustering.method = "average"
)

cat("Saved -> ", file.path(out_dir, "heatmap_ordered_by_celltype.pdf"), "\n")