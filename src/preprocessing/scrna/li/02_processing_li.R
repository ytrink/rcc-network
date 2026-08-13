# preprocessing/scrna/li/02_processing_li.R
#
# Loads Li et al. raw count matrix, creates Seurat object, runs SCTransform
# per patient, and integrates with Harmony.
#
# Input:  data/li/matrix/adata_X.mtx, barcodes.csv, genes.csv, metadata.csv
# Output: data/li/Li_integrated.rds

library(Seurat)
library(harmony)
library(Matrix)
library(dplyr)
library(future)
library(here)

options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB
set.seed(101)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir   <- here::here()
matrix_dir <- file.path(root_dir, "data", "li", "matrix")
out_dir    <- file.path(root_dir, "data", "li")

# ── Parameters ─────────────────────────────────────────────────────────────────
PCA_DIMS    <- 1:30
RESOLUTIONS <- c(0.2, 0.3, 0.5, 0.8)

# ── Load matrix ────────────────────────────────────────────────────────────────
expr_matrix           <- readMM(file.path(matrix_dir, "adata_X.mtx"))
rownames(expr_matrix) <- readLines(file.path(matrix_dir, "barcodes.csv"))
colnames(expr_matrix) <- readLines(file.path(matrix_dir, "genes.csv"))
metadata              <- read.csv(file.path(matrix_dir, "metadata.csv"))

# ── Create Seurat object ───────────────────────────────────────────────────────
seurat_obj             <- CreateSeuratObject(counts = t(expr_matrix), meta.data = metadata)
seurat_obj@meta.data   <- metadata
rm(expr_matrix, metadata)

# ── SCTransform per patient ────────────────────────────────────────────────────
seurat_list <- SplitObject(seurat_obj, split.by = "patient")
seurat_list <- lapply(seurat_list, SCTransform, verbose = FALSE)
rm(seurat_obj)

# ── Merge ──────────────────────────────────────────────────────────────────────
rna.li <- merge(
  seurat_list[[1]],
  y            = seurat_list[-1],
  add.cell.ids = names(seurat_list)
)
rm(seurat_list)

DefaultAssay(rna.li) <- "SCT"

# ── PCA ────────────────────────────────────────────────────────────────────────
rna.li <- rna.li %>%
  FindVariableFeatures(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

# ── Harmony integration ────────────────────────────────────────────────────────
rna.li <- RunHarmony(rna.li, group.by.vars = "patient", plot_convergence = FALSE)

# ── Clustering ─────────────────────────────────────────────────────────────────
rna.li <- rna.li %>%
  RunUMAP(reduction = "harmony", dims = PCA_DIMS, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = PCA_DIMS, verbose = FALSE)

for (res in RESOLUTIONS) {
  rna.li <- FindClusters(rna.li, resolution = res, verbose = FALSE)
  rna.li[[paste0("clusters_res_", res)]] <- Idents(rna.li)
}

# ── Save ───────────────────────────────────────────────────────────────────────
saveRDS(rna.li, file.path(out_dir, "Li_integrated.rds"))
cat("Saved -> ", file.path(out_dir, "Li_integrated.rds"), "\n")