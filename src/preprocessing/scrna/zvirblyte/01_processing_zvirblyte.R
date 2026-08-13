# preprocessing/scrna/zvirblyte/01_processing_zvirblyte.R
#
# Loads Zvirblyte et al. h5ad, adds broad cell type annotations,
# and processes endothelial and stromal subsets with Harmony integration.
#
# Input:  data/zvirblyte/GSE242299_all_cells_50236_33538.h5ad
# Output: data/zvirblyte/Zvir_all.rds
#         data/zvirblyte/Zvir_endo.rds
#         data/zvirblyte/Zvir_stromal.rds

library(Seurat)
library(SeuratDisk)
library(harmony)
library(here)

set.seed(101)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir <- here::here()
data_dir <- file.path(root_dir, "data", "zvirblyte")

# ── Parameters ─────────────────────────────────────────────────────────────────
PCA_DIMS   <- 1:30
RESOLUTION <- 0.5

# ── Broad cell type mapping ────────────────────────────────────────────────────
BROAD_CELLTYPE <- c(
  "TAM 1"                            = "Macrophage/Monocyte",
  "TAM 2"                            = "Macrophage/Monocyte",
  "TAM 3"                            = "Macrophage/Monocyte",
  "TAM 4"                            = "Macrophage/Monocyte",
  "Classical monocytes"              = "Macrophage/Monocyte",
  "Non-classical monocytes"          = "Macrophage/Monocyte",
  "CD8 T cells"                      = "T_cells",
  "Resting/memory T cells"           = "T_cells",
  "Cytotoxic T cells"                = "T_cells",
  "Regulatory T cells"               = "T_cells",
  "NK cells"                         = "NK_cells",
  "B cells"                          = "B_cells",
  "Plasma cells"                     = "B_cells",
  "IGHG-high plasma cells"           = "B_cells",
  "Mast cells"                       = "Other_immune",
  "Tumor cells 1"                    = "Tumor",
  "Tumor cells 2"                    = "Tumor",
  "Tumor cells 3"                    = "Tumor",
  "Tumor vasculature 1"              = "Endothelial",
  "Tumor vasculature 2"              = "Endothelial",
  "Tumor vasculature 3"              = "Endothelial",
  "Tumor vasculature 4"              = "Endothelial",
  "Tumor AVR-like vasculature"       = "Endothelial",
  "Glomerular endothelium"           = "Endothelial",
  "AVR"                              = "Endothelial",
  "DVR"                              = "Endothelial",
  "vSMCs"                            = "Stromal",
  "Mesangial/vSMCs"                  = "Stromal",
  "Myofibroblasts"                   = "Stromal",
  "Proximal tubule"                  = "Epithelial",
  "Epithelial progenitor-like cells" = "Epithelial",
  "Podocytes"                        = "Epithelial",
  "Principal cells"                  = "Epithelial",
  "Type A-ICs"                       = "Epithelial",
  "OM Type A-ICs"                    = "Epithelial",
  "Type B-IC"                        = "Epithelial",
  "TAL of LOH"                       = "Epithelial",
  "tAL of LOH"                       = "Epithelial",
  "DCT/CNT"                          = "Epithelial",
  "Cycling"                          = "Other"
)

# ── Helper: SCTransform + Harmony + cluster ────────────────────────────────────
process_subset <- function(obj, batch_var = "sample",
                           pca_dims = PCA_DIMS, resolution = RESOLUTION) {
  DefaultAssay(obj) <- "RNA"
  obj %>%
    SCTransform(assay = "RNA", verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunHarmony(group.by.vars = batch_var, verbose = FALSE) %>%
    RunUMAP(reduction = "harmony", dims = pca_dims, verbose = FALSE) %>%
    FindNeighbors(reduction = "harmony", dims = pca_dims, verbose = FALSE) %>%
    FindClusters(resolution = resolution, verbose = FALSE)
}

# ── Load and convert ───────────────────────────────────────────────────────────
h5ad_path    <- file.path(data_dir, "GSE242299_all_cells_50236_33538.h5ad")
h5seurat_path <- sub(".h5ad", ".h5Seurat", h5ad_path)

if (!file.exists(h5seurat_path)) {
  Convert(h5ad_path, dest = "h5Seurat")
}

seurat_obj <- LoadH5Seurat(h5seurat_path)

# ── Prefix barcodes with patient/library info ──────────────────────────────────
colnames(seurat_obj) <- paste0(
  seurat_obj$library2, "_",
  seurat_obj$library,  "_",
  colnames(seurat_obj)
)

# ── Remove pre-computed reductions ────────────────────────────────────────────
for (red in c("draw_graph_fa", "pca_harmony", "umap", "pca")) {
  if (red %in% names(seurat_obj@reductions)) {
    seurat_obj@reductions[[red]] <- NULL
  }
}

# ── Add broad cell type annotations ───────────────────────────────────────────
unmapped <- setdiff(unique(seurat_obj$cell_type), names(BROAD_CELLTYPE))
if (length(unmapped)) {
  message("Unmapped cell type labels: ", paste(unmapped, collapse = ", "))
}

seurat_obj$broad_cell_type_2 <- factor(
  unname(BROAD_CELLTYPE[seurat_obj$cell_type]),
  levels = c("B_cells", "T_cells", "NK_cells", "Macrophage/Monocyte",
             "Endothelial", "Stromal", "Epithelial", "Tumor",
             "Other_immune", "Other")
)

seurat_obj <- NormalizeData(seurat_obj)

# ── Save full object and metadata ─────────────────────────────────────────────
saveRDS(seurat_obj, file.path(data_dir, "Zvir_all.rds"))
write.csv(seurat_obj@meta.data, file.path(data_dir, "metadata_all_cells.csv"))
cat("Saved full object\n")

# ── Process endothelial subset ────────────────────────────────────────────────
endo <- process_subset(subset(seurat_obj, subset = cell_group == "Endothelium"))
saveRDS(endo, file.path(data_dir, "Zvir_endo.rds"))
cat("Saved endothelial subset\n")

# ── Process stromal subset ────────────────────────────────────────────────────
stromal <- process_subset(subset(seurat_obj, subset = cell_group == "Stromal"))
saveRDS(stromal, file.path(data_dir, "Zvir_stromal.rds"))
cat("Saved stromal subset\n")