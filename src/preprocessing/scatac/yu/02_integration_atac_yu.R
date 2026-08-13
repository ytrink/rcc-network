# preprocessing/scatac/yu/02_integration_atac_yu.R
#
# TF-IDF normalization, LSI dimensionality reduction, Harmony batch correction,
# and clustering of merged Yu et al. scATAC-seq data.
#
# Input:  data/yu/scatac/atac_merged_filtered.rds
# Output: data/yu/scatac/atac_integrated.rds

library(Seurat)
library(Signac)
library(harmony)
library(ggplot2)
library(here)

set.seed(101)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir <- here::here()
data_dir <- file.path(root_dir, "data", "yu", "scatac")
qc_dir   <- file.path(root_dir, "preprocessing", "scatac", "yu", "qc_plots")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# ── Parameters ─────────────────────────────────────────────────────────────────
LSI_DIMS    <- 2:30
RESOLUTIONS <- seq(0.2, 1.2, by = 0.1)

# ── Load ───────────────────────────────────────────────────────────────────────
atac <- readRDS(file.path(data_dir, "atac_merged_filtered.rds"))
DefaultAssay(atac) <- "ATAC"

# ── Normalization and dimensionality reduction ─────────────────────────────────
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q1")
atac <- RunSVD(atac)

# ── Pre-integration QC plots ───────────────────────────────────────────────────
pdf(file.path(qc_dir, "pre_integration_atac.pdf"))
print(ElbowPlot(atac, ndims = 30, reduction = "lsi"))
print(DepthCor(atac))
print(DimPlot(atac, group.by = "batch", pt.size = 0.1,
              reduction = "umap") + ggtitle("Unintegrated"))
dev.off()

# ── Harmony integration ────────────────────────────────────────────────────────
# Note: LSI reduction is copied to 'pca' slot as workaround for Harmony
atac@reductions$pca <- atac@reductions$lsi

atac <- RunHarmony(
  atac,
  group.by.vars = "batch",
  assay.use     = "ATAC",
  project.dim   = FALSE
)

atac <- RunUMAP(atac, dims = LSI_DIMS, reduction = "harmony", verbose = FALSE)

# ── Post-integration plot ──────────────────────────────────────────────────────
pdf(file.path(qc_dir, "post_integration_atac.pdf"))
print(DimPlot(atac, group.by = "batch", pt.size = 0.1) + ggtitle("Harmony integration"))
dev.off()

# ── Clustering ─────────────────────────────────────────────────────────────────
atac <- FindNeighbors(atac, reduction = "harmony", dims = LSI_DIMS, verbose = FALSE)

for (res in RESOLUTIONS) {
  atac <- FindClusters(atac, resolution = res, algorithm = 3, verbose = FALSE)
  atac[[paste0("clusters_res_", res)]] <- Idents(atac)
}

# ── Save ───────────────────────────────────────────────────────────────────────
saveRDS(atac, file.path(data_dir, "atac_integrated.rds"))
cat("Saved ->", file.path(data_dir, "atac_integrated.rds"), "\n")