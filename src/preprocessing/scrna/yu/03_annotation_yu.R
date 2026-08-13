# preprocessing/scrna/yu/03_annotation_yu.R
#
# Marker gene feature plots for manual cell type annotation across
# multiple clustering resolutions (Yu et al. scRNA-seq).
#
# Input:  data/yu/scrna/rna_integrated.rds
# Output: preprocessing/scrna/yu/annotation_plots/<resolution>/<cell_type>.pdf

library(Seurat)
library(here)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir <- here::here()
data_dir <- file.path(root_dir, "data", "yu", "scrna")
plot_dir <- file.path(root_dir, "preprocessing", "scrna", "yu", "annotation_plots")

# ── Marker genes by cell type ──────────────────────────────────────────────────
MARKER_GENES <- list(
  immune            = c("PTPRC"),
  t_cell            = c("CD3D", "CD3E"),
  cd8_t_cell        = c("CD8A", "CD8B"),
  cd4_t_cell        = c("CD4", "CD40LG", "IL7R"),
  treg              = c("FOXP3", "IL2RA"),
  nk_cell           = c("GNLY", "KLRD1"),
  proliferating_t   = c("CD8A", "MKI67", "TOP2A", "STMN1"),
  b_cell            = c("MS4A1", "CD79A", "MZB1", "IGKC", "JCHAIN"),
  tumor             = c("CA9", "NDUFA4L2"),
  tumor_other       = c("BAP1", "PAX8", "PAX2", "MKI67", "HIF1A", "VHL"),
  macrophage        = c("CSF1R", "CD68", "CD163"),
  monocyte          = c("S100A8", "S100A9", "S100A12", "FCGR3A", "LST1", "LILRB2"),
  dendritic         = c("CD1C", "CLEC10A", "CLEC9A", "IDO1"),
  neutrophil        = c("CXCL8", "FCGR3B"),
  mast              = c("TPSAB1", "KIT"),
  endothelial       = c("PTPRB", "PECAM1", "KDR", "VCAM1"),
  mesangial         = c("PDGFRB", "ACTA2"),
  fibroblast        = c("COL1A1", "COL1A2", "FBLN1", "FBLN2")
)

RESOLUTIONS <- c(0.2, 0.3, 0.5, 0.8)

# ── Load ───────────────────────────────────────────────────────────────────────
data.combined <- readRDS(file.path(data_dir, "rna_integrated.rds"))

# ── Plot marker genes per resolution ──────────────────────────────────────────
for (res in RESOLUTIONS) {
  res_label <- gsub("\\.", "", as.character(res))  # e.g. 0.2 -> "02"
  res_dir   <- file.path(plot_dir, paste0("res_", res_label))
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  
  Idents(data.combined) <- data.combined[[paste0("clusters_res_", res)]]
  
  # UMAP overview
  pdf(file.path(res_dir, "basic.pdf"), width = 20, height = 12)
  print(DimPlot(data.combined, group.by = "ident",  label = TRUE, raster = TRUE))
  print(DimPlot(data.combined, group.by = "batch",  label = TRUE, raster = TRUE))
  dev.off()
  
  # Marker gene feature plots
  for (cell_type in names(MARKER_GENES)) {
    pdf(file.path(res_dir, paste0(cell_type, ".pdf")), width = 20, height = 12)
    print(FeaturePlot(data.combined, features = MARKER_GENES[[cell_type]], raster = TRUE))
    dev.off()
  }
  
  cat("Done resolution", res, "\n")
}