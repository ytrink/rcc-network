# src/figure4/R/01_umap.R
#
# Endothelial cells (Zvirblyte et al.): healthy EC vs tumor EC (TEC) embedding and UMAP panels.
# Port of the original analysis (plot_regulon_communities_LI/endo/renin_score_v3.R, lines 1-35 and 250-492)
# and plot_umap_panels_consistent.R.
#
#   1. Keep healthy ECs (AVR, DVR, glomerular endothelium) and TECs (tumor vasculature 1-4, tumor AVR-like)
#   2. SCTransform -> PCA -> Harmony (by patient) -> UMAP on Harmony 1:10
#   3. TEC score = AddModuleScore(PLVAP, VWF, SPARC, INSR, ANGPT2) on SCT
#   4. Export matrices used by the perturbation scripts (figure4/python) and draw the UMAP panels
#      of Fig. 4A, S12A-C and S15A-B with identical axes (Seurat DimPlot/FeaturePlot style, as Figs. 2A/3A)
#
# Input:  data/zvirblyte/zvir_endo.rds
#         optional: data/zvirblyte/X_umap_published.csv  (published coordinates; if present they are used
#                   for plotting so the panels match the paper exactly - UMAP can differ slightly between
#                   package versions)
# Output: results/figure4/ zvir_endo_norm.rds, X_umap.csv, X_counts.csv, metadata.csv, umap_panels/*.png|pdf

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(here)

source(here::here("src", "utils", "paths.R"))
out_dir   <- results_dir("figure4")
panel_dir <- results_dir("figure4", "umap_panels")

# ── Parameters ─────────────────────────────────────────────────────────────────
HEALTHY_EC_TYPES <- c("AVR", "DVR", "Glomerular endothelium")
TUMOR_EC_TYPES   <- c(
  "Tumor vasculature 1", "Tumor vasculature 2",
  "Tumor vasculature 3", "Tumor vasculature 4",
  "Tumor AVR-like vasculature"
)
TEC_MODULE <- c("PLVAP", "VWF", "SPARC", "INSR", "ANGPT2")

# ── Subset ─────────────────────────────────────────────────────────────────────
endo <- readRDS(data_path("zvirblyte", "zvir_endo.rds"))
DefaultAssay(endo) <- "RNA"
endo$pb_group <- case_when(
  endo$cell_type %in% HEALTHY_EC_TYPES ~ "Healthy_EC",
  endo$cell_type %in% TUMOR_EC_TYPES   ~ "Tumor_EC",
  TRUE ~ NA_character_
)
endo.pb <- subset(endo, cells = rownames(endo@meta.data)[!is.na(endo$pb_group)])

# ── Re-embedding (SCTransform -> PCA -> Harmony -> UMAP) ───────────────────────
endo.pb <- SCTransform(endo.pb, assay = "RNA", verbose = FALSE)
endo.pb <- RunPCA(endo.pb, assay = "SCT", verbose = FALSE)
endo.pb <- RunHarmony(endo.pb, group.by.vars = "patient")
endo.pb <- RunUMAP(endo.pb, reduction = "harmony", dims = 1:10)   # Seurat default seed.use = 42
endo.pb <- FindNeighbors(endo.pb, reduction = "harmony", dims = 1:15)

# ── TEC score ──────────────────────────────────────────────────────────────────
endo.pb <- AddModuleScore(endo.pb, features = list(TEC_MODULE), name = "zvir_module", assay = "SCT")
endo.pb$tec_score <- endo.pb$zvir_module1

# ── Exports for the perturbation step (figure4/python) ────────────────────────
saveRDS(endo.pb, file.path(out_dir, "zvir_endo_norm.rds"))
write.csv(endo.pb@meta.data, file.path(out_dir, "metadata.csv"))
write.csv(as.matrix(Embeddings(endo.pb, "umap")), file.path(out_dir, "X_umap.csv"))

endo.rna <- JoinLayers(endo.pb, assay = "RNA")
X_counts <- GetAssayData(endo.rna, assay = "RNA", layer = "counts")
write.csv(as.data.frame(as.matrix(X_counts)), file.path(out_dir, "X_counts.csv"))   # regressors were trained on raw counts

# ── UMAP panels (Fig. 4A, S12A-C, S15A-B) ─────────────────────────────────────
published <- data_path("zvirblyte", "X_umap_published.csv")
if (file.exists(published)) {
  umap <- read.csv(published, row.names = 1)
  stopifnot(all(colnames(endo.pb) %in% rownames(umap)))
  endo.pb[["umap"]] <- CreateDimReducObject(as.matrix(umap[colnames(endo.pb), ]), key = "UMAP_", assay = "SCT")
  cat("Using published UMAP coordinates for the panels\n")
}

endo.pb$broad <- factor(
  ifelse(endo.pb$pb_group == "Tumor_EC", "Tumor endothelial cells\n(TECs)", "Endothelial cells (ECs)"),
  levels = c("Endothelial cells (ECs)", "Tumor endothelial cells\n(TECs)")
)
endo.pb$cell_type2 <- factor(endo.pb$cell_type, levels = c(HEALTHY_EC_TYPES, "Tumor AVR-like vasculature",
                                                           paste("Tumor vasculature", 1:4)))

# Same axis limits for every panel (so ticks match), 1:1 aspect; titles are added in the figure layout
emb  <- Embeddings(endo.pb, "umap")
XLIM <- range(emb[, 1]) + c(-0.4, 0.4)
YLIM <- range(emb[, 2]) + c(-0.4, 0.4)

finish <- function(p) {
  p + ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2") +
    scale_x_continuous(limits = XLIM) + scale_y_continuous(limits = YLIM) +
    coord_fixed(ratio = 1)
}
discrete <- list(guides(colour = guide_legend(override.aes = list(size = 5))),
                 theme(legend.text = element_text(size = 14)))
save_panel <- function(p, name, width = 6.5, height = 5.5) {
  ggsave(file.path(panel_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(panel_dir, paste0(name, ".pdf")), p, width = width, height = height, bg = "white")
}

save_panel(finish(DimPlot(endo.pb, reduction = "umap", group.by = "broad")) + discrete,
           "Fig4A_S12B_S15A_broad_annotation", width = 8.5)
save_panel(finish(FeaturePlot(endo.pb, features = "tec_score", reduction = "umap")),
           "Fig4A_S15B_TEC_score")
save_panel(finish(DimPlot(endo.pb, reduction = "umap", group.by = "cell_type2")) + discrete,
           "S12A_zvirblyte_annotations", width = 9)
for (gene in TEC_MODULE) {
  save_panel(finish(FeaturePlot(endo.pb, features = gene, reduction = "umap")), paste0("S12C_", gene))
}

cat("Saved ->", out_dir, "\n")
