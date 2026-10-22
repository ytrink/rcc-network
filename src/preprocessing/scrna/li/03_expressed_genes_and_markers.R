# src/preprocessing/scrna/li/03_expressed_genes_and_markers.R
#
# Small gene lists derived from the Li et al. data that are used downstream:
#   1. {mac_mono,CD8}_expressed_genes_10perc.csv - genes detected (counts > 0) in >= 10% of cells of the
#      cell type (centrality subnetworks, Figs. S6/S10; CD8 perturbation). Method as in
#      D:/SC_datasets/ccRCC_li/get_expressed_genes_per_celltype.R (re-computation matches the published lists:
#      mac/mono identical, CD8 Jaccard 0.988).
#   2. degs_per_cell_type.csv - marker genes per broad cell type (B_cell, Endothelial, Fibroblast, Mac_mono,
#      T_NK, Tumor), used to train the perturbation regressors (Figs. 2-3). The original script was not saved;
#      FindAllMarkers parameters were inferred from the published file (no negative log2FC, min log2FC = 1.0,
#      min pct.1 = 0.2). Please use the published file (Zenodo) to reproduce the figures exactly.
#
# Input:  data/li/Li_mac_mono.rds, data/li/Li_CD8T_normalized.rds, data/li/Li_integrated.rds
# Output: data/li/mac_mono_expressed_genes_10perc.csv, data/li/CD8_expressed_genes_10perc.csv,
#         data/li/degs_per_cell_type.csv (only written if not already present)

library(Seurat)
library(dplyr)
library(here)

source(here::here("src", "utils", "paths.R"))

expressed_10perc <- function(obj) {
  counts <- GetAssayData(JoinLayers(obj, assay = "RNA"), assay = "RNA", layer = "counts")
  names(which(Matrix::rowMeans(counts > 0) >= 0.10))
}

for (ct in c("mac_mono", "CD8")) {
  f <- data_path("li", paste0(ct, "_expressed_genes_10perc.csv"))
  if (file.exists(f)) next
  obj <- readRDS(data_path("li", if (ct == "mac_mono") "Li_mac_mono.rds" else "Li_CD8T_normalized.rds"))
  write.csv(expressed_10perc(obj), f)
}

markers_file <- data_path("li", "degs_per_cell_type.csv")
if (!file.exists(markers_file)) {
  li <- readRDS(data_path("li", "Li_integrated.rds"))
  li$annotation_broad_1 <- case_when(
    li$annotation_broad %in% c("CD4_T", "CD8_T", "NK")   ~ "T_NK",
    li$annotation_broad %in% c("Macrophage", "Monocyte") ~ "Mac_mono",
    TRUE ~ li$annotation_broad
  )
  li <- subset(li, subset = annotation_broad_1 %in% c("B_cell", "Endothelial", "Fibroblast", "Mac_mono", "T_NK", "Tumor"))
  DefaultAssay(li) <- "RNA"
  li <- NormalizeData(JoinLayers(li))
  Idents(li) <- li$annotation_broad_1
  degs <- FindAllMarkers(li, only.pos = TRUE, logfc.threshold = 1, min.pct = 0.2)
  write.csv(degs, markers_file, row.names = FALSE)
}
cat("Done\n")
