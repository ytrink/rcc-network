# preprocessing/scrna/yu/02_integration_yu.R
#
# Harmony integration and clustering of merged Yu et al. scRNA-seq data.
#
# Input:  data/yu/scrna/rna_merged.rds
# Output: data/yu/scrna/rna_integrated.rds

library(Seurat)
library(harmony)
library(dplyr)

set.seed(101)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir <- here::here()
data_dir <- file.path(root_dir, "data", "yu", "scrna")

# ── Parameters ─────────────────────────────────────────────────────────────────
PCA_DIMS    <- 1:30
RESOLUTIONS <- c(0.2, 0.3, 0.5, 0.8)

# ── Load ───────────────────────────────────────────────────────────────────────
ccRCC.merged <- readRDS(file.path(data_dir, "rna_merged.rds"))
DefaultAssay(ccRCC.merged) <- "SCT"

# ── PCA ────────────────────────────────────────────────────────────────────────
ccRCC.merged <- ccRCC.merged %>%
  FindVariableFeatures() %>%
  RunPCA(verbose = FALSE)

# ── Harmony integration ────────────────────────────────────────────────────────
ccRCC.rna <- RunHarmony(
  ccRCC.merged,
  group.by.vars    = "batch",
  plot_convergence = FALSE
)

# ── Clustering at multiple resolutions ────────────────────────────────────────
ccRCC.rna <- ccRCC.rna %>%
  RunUMAP(reduction = "harmony", dims = PCA_DIMS, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = PCA_DIMS, verbose = FALSE)

for (res in RESOLUTIONS) {
  ccRCC.rna <- FindClusters(ccRCC.rna, resolution = res, verbose = FALSE)
  ccRCC.rna[[paste0("clusters_res_", res)]] <- Idents(ccRCC.rna)
}

# ── Save ───────────────────────────────────────────────────────────────────────
saveRDS(ccRCC.rna, file.path(data_dir, "rna_integrated.rds"))
cat("Saved integrated object ->", file.path(data_dir, "rna_integrated.rds"), "\n")