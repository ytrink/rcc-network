# preprocessing/scatac/yu/03_annotation_transfer_yu.R
#
# Transfers cell type labels from annotated scRNA-seq reference to
# scATAC-seq data using Seurat label transfer (CCA anchors).
#
# Input:  data/yu/scatac/atac_integrated.rds
#         data/yu/scrna/rna_annotated.h5Seurat
# Output: data/yu/scatac/atac_annotated.rds

library(Seurat)
library(Signac)
library(SeuratDisk)
library(ggplot2)
library(here)

options(future.globals.maxSize = 8 * 1024^3)  # 8 GB
set.seed(123)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir  <- here::here()
atac_dir  <- file.path(root_dir, "data", "yu", "scatac")
rna_dir   <- file.path(root_dir, "data", "yu", "scrna")
plot_dir  <- file.path(root_dir, "preprocessing", "scatac", "yu", "qc_plots")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load ───────────────────────────────────────────────────────────────────────
cat("Loading ATAC object...\n")
data.atac <- readRDS(file.path(atac_dir, "atac_integrated.rds"))

cat("Loading RNA reference...\n")
data.rna  <- LoadH5Seurat(file.path(rna_dir, "rna_annotated.h5Seurat"))
DefaultAssay(data.rna) <- "SCT"

# ── Find transfer anchors ──────────────────────────────────────────────────────
cat("Finding transfer anchors...\n")
transfer.anchors <- FindTransferAnchors(
  reference      = data.rna,
  query          = data.atac,
  features       = VariableFeatures(data.rna),
  reference.assay = "RNA",
  query.assay    = "ACTIVITY",
  reduction      = "rpca",
  verbose        = TRUE
)

saveRDS(transfer.anchors, file.path(atac_dir, "transfer_anchors.rds"))

# ── Transfer labels ────────────────────────────────────────────────────────────
cat("Transferring cell type labels...\n")
celltype.predictions <- TransferData(
  anchorset        = transfer.anchors,
  refdata          = data.rna$cell_type.2,
  weight.reduction = data.atac[["lsi"]],
  dims             = 2:30
)

data.atac <- AddMetaData(data.atac, metadata = celltype.predictions)
Idents(data.atac) <- data.atac$predicted.id

# ── Plot ───────────────────────────────────────────────────────────────────────
pdf(file.path(plot_dir, "atac_predicted_celltypes.pdf"))
print(DimPlot(data.atac, group.by = "predicted.id", label = TRUE) +
        NoLegend() + ggtitle("Predicted cell type annotation"))
dev.off()

# ── Save ───────────────────────────────────────────────────────────────────────
saveRDS(data.atac, file.path(atac_dir, "atac_annotated.rds"))
cat("Saved ->", file.path(atac_dir, "atac_annotated.rds"), "\n")