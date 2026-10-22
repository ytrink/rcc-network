# src/preprocessing/monocle/01_run_monocle3_mac_mono.R
#
# Monocle3 embedding of Li et al. macrophages/monocytes (UMAP of Fig. 2A/F) and the cell subsets used by
# the perturbation analysis (Fig. 2E/F).
# Port of stuff_to_run_on_the_side/monocle_scplus/run_monocle_mac_mono.R (object cds_2.rds, verified:
# 44,448 cells, 50 PCA/aligned dims) and plot_regulon_communities_LI/../monocle/mac_mono/Li/run_monocle.R (exports).
#
# Note: Monocle's exprs() returns raw counts, so the exported "sampled" matrices are raw counts
# (these were used to train and run the perturbation regressors).
#
# Input:  data/li/Li_mac_mono.rds, data/li/metadata_merged.csv
# Output: results/monocle/mac_mono/cds_mac_mono.rds
#         data/monocle/cds_mac_mono_umap_sampled30k.csv                       (UMAP for Fig. 2A)
#         data/perturbation/mac_mono/cds_mac_mono_sampled15k.csv              (regressor training)
#         data/perturbation/mac_mono/cds_mac_mono_sampled10k.csv, cds_mac_mono_umap_sampled10k.csv (simulation, arrows)

library(monocle3)
library(Seurat)
library(data.table)
library(here)

source(here::here("src", "utils", "paths.R"))
out_dir  <- results_dir("monocle", "mac_mono")
umap_dir <- data_path("monocle");                  dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
pert_dir <- data_path("perturbation", "mac_mono"); dir.create(pert_dir, recursive = TRUE, showWarnings = FALSE)

# ── Cell data set (raw counts + Li metadata) ───────────────────────────────────
rna.li <- readRDS(data_path("li", "Li_mac_mono.rds"))
colnames(rna.li) <- paste0("LI_", rna.li$patient, "_", colnames(rna.li))
rna.li <- JoinLayers(rna.li)
expression.matrix <- GetAssayData(rna.li, assay = "RNA", layer = "counts")

sample_metadata <- read.csv(data_path("li", "metadata_merged.csv"))
rownames(sample_metadata) <- sample_metadata$barcode
cell_meta <- sample_metadata[colnames(expression.matrix), ]
gene_meta <- data.frame(gene_short_name = rownames(expression.matrix), row.names = rownames(expression.matrix))

cds <- new_cell_data_set(expression.matrix, cell_metadata = cell_meta, gene_metadata = gene_meta)

# ── Embedding (parameters of the published object cds_2.rds) ──────────────────
cds <- preprocess_cds(cds, num_dim = 50, method = "PCA")
cds <- align_cds(cds, alignment_group = "patient")
cds <- reduce_dimension(cds, umap.min_dist = 0.3, umap.n_neighbors = 50L)
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE)
saveRDS(cds, file.path(out_dir, "cds_mac_mono.rds"))

# ── Exports (each with set.seed(123), as in the original) ─────────────────────
export_sample <- function(n, umap_file = NULL, expr_file = NULL) {
  set.seed(123)
  cds_s <- cds[, sample(colnames(cds), size = n)]
  if (!is.null(umap_file)) write.csv(as.data.frame(reducedDims(cds_s)$UMAP), umap_file)
  if (!is.null(expr_file)) {
    m <- exprs(cds_s)
    fwrite(cbind(gene = rownames(m), as.data.frame(as.matrix(m))), file = expr_file)
  }
}
export_sample(30000, umap_file = file.path(umap_dir, "cds_mac_mono_umap_sampled30k.csv"))
export_sample(10000, umap_file = file.path(pert_dir, "cds_mac_mono_umap_sampled10k.csv"),
                     expr_file = file.path(pert_dir, "cds_mac_mono_sampled10k.csv"))
export_sample(15000, expr_file = file.path(pert_dir, "cds_mac_mono_sampled15k.csv"))

cat("Saved ->", out_dir, umap_dir, pert_dir, "\n")
