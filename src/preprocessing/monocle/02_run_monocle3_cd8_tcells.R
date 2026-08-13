# preprocessing/monocle/02_run_monocle3_cd8_tcells.R
#
# Monocle3 cell embedding, UMAP reduction, trajectory analysis, and
# pseudotime ordering for CD8+ T-cells (Li et al. cohort).
#
# Input:  data/Li_TNK.rds
# Output: preprocessing/monocle/cd8_tcell/cds_tcells_umap_sampled30k.csv
#         preprocessing/monocle/cd8_tcell/cds_normalized_tcells_sampled30k.csv
#         figures/figure3/umap_cd8_exhaustion.pdf

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
out_dir  <- file.path(root_dir, "preprocessing", "monocle", "cd8_tcell")
fig_dir  <- file.path(root_dir, "figures", "figure3")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load Seurat Object ─────────────────────────────────────────────────────────
cat("Loading Li et al. T/NK cell data...\n")
rna_path <- file.path(data_dir, "Li_TNK.rds")
if (!file.exists(rna_path)) {
  rna_path <- file.path(data_dir, "Li_T.rds")
}
rna.li <- readRDS(rna_path)
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
write.csv(umap_30k, file.path(out_dir, "cds_tcells_umap_sampled30k.csv"))

# 30k sampled expression matrix
expr_30k <- exprs(cds_30k)
fwrite(
  cbind(gene = rownames(expr_30k), as.data.frame(as.matrix(expr_30k))),
  file = file.path(out_dir, "cds_normalized_tcells_sampled30k.csv")
)

# Save processed Monocle3 object
saveRDS(cds, file.path(out_dir, "cds_tcells_ordered.rds"))
cat("Saved Monocle3 object and sampling outputs ->", out_dir, "\n")
