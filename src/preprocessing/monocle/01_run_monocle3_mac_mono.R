# preprocessing/monocle/01_run_monocle3_mac_mono.R
#
# Monocle3 cell embedding, UMAP reduction, trajectory analysis, and
# pseudotime ordering for Macrophage/Monocyte cells (Li et al. cohort).
#
# Input:  data/Li_mac_mono.rds
# Output: preprocessing/monocle/mac_mono/cds_mac_mono_umap_sampled30k.csv
#         preprocessing/monocle/mac_mono/cds_mac_mono_sampled15k.csv
#         figures/figure2/umap_macro_groups.pdf

library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(here)

set.seed(101)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir <- here::here()
data_dir <- file.path(root_dir, "data")
out_dir  <- file.path(root_dir, "preprocessing", "monocle", "mac_mono")
fig_dir  <- file.path(root_dir, "figures", "figure2")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load Seurat Object ─────────────────────────────────────────────────────────
cat("Loading Li et al. Macrophage/Monocyte data...\n")
rna.li <- readRDS(file.path(data_dir, "Li_mac_mono.rds"))
DefaultAssay(rna.li) <- "RNA"

# ── Initialize Monocle3 CellDataSet ───────────────────────────────────────────
cat("Constructing Monocle3 CellDataSet...\n")
counts_mat  <- GetAssayData(rna.li, slot = "counts")
cell_meta   <- rna.li@meta.data
gene_meta   <- data.frame(gene_short_name = rownames(counts_mat), row.names = rownames(counts_mat))

cds <- new_cell_data_set(
  expression_data = counts_mat,
  cell_metadata   = cell_meta,
  gene_metadata   = gene_meta
)

# ── Monocle3 Preprocessing, Alignment & Dimensionality Reduction ──────────────
cat("Running preprocess_cds(), align_cds(), and reduce_dimension()...\n")
cds <- preprocess_cds(cds, num_dim = 30)
if ("patient" %in% colnames(colData(cds))) {
  cds <- align_cds(cds, alignment_group = "patient")
}
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# ── Pseudotime Ordering ────────────────────────────────────────────────────────
cat("Ordering cells along pseudotime trajectory...\n")
cds <- order_cells(cds)

# ── Export Downsampled Datasets for downstream analysis & figures ─────────────
cat("Exporting downsampled UMAP coordinates and expression matrices...\n")
set.seed(123)

# 30k sampled UMAP coordinates
cell_ids_30k <- sample(colnames(cds), size = min(30000, ncol(cds)))
cds_30k      <- cds[, cell_ids_30k]
umap_30k     <- as.data.frame(reducedDims(cds_30k)$UMAP)
colnames(umap_30k) <- c("UMAP_1", "UMAP_2")
write.csv(umap_30k, file.path(out_dir, "cds_mac_mono_umap_sampled30k.csv"))

# 15k sampled expression matrix (for perturbation regressor models)
cell_ids_15k <- sample(colnames(cds), size = min(15000, ncol(cds)))
cds_15k      <- cds[, cell_ids_15k]
expr_15k     <- exprs(cds_15k)

fwrite(
  cbind(gene = rownames(expr_15k), as.data.frame(as.matrix(expr_15k))),
  file = file.path(out_dir, "cds_mac_mono_sampled15k.csv")
)

# Save processed Monocle3 object
saveRDS(cds, file.path(out_dir, "cds_mac_mono_ordered.rds"))
cat("Saved Monocle3 object and sampling outputs ->", out_dir, "\n")
