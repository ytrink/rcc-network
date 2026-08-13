# figures/figure4/R/01_umap.R
#
# Generates UMAP plots for normal-to-tumor endothelial cell transition (Figure 4).
# Plots: (1) tumor vasculature module score, (2) Healthy/Tumor EC group assignment
#
# Input:  data/zvirblyte/Zvir_endo.rds
# Output: figures/figure4/R/umap_tumor_ec_module.pdf
#         figures/figure4/R/umap_ec_groups.pdf

library(Seurat)
library(dplyr)
library(ggplot2)
library(here)

set.seed(101)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir <- here::here()
data_dir <- file.path(root_dir, "data")
out_dir  <- file.path(root_dir, "figures", "figure4", "R")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Parameters ─────────────────────────────────────────────────────────────────
HEALTHY_EC_TYPES <- c("AVR", "DVR", "Glomerular endothelium")
TUMOR_EC_TYPES   <- c(
  "Tumor vasculature 1", "Tumor vasculature 2",
  "Tumor vasculature 3", "Tumor vasculature 4",
  "Tumor AVR-like vasculature"
)

TUMOR_EC_MODULE <- c("PLVAP", "VWF", "SPARC", "INSR", "ANGPT2")

# ── Load and subset ────────────────────────────────────────────────────────────
endo <- readRDS(file.path(data_dir, "zvirblyte", "Zvir_endo.rds"))
DefaultAssay(endo) <- "RNA"
endo <- NormalizeData(endo, verbose = FALSE)
endo$pb_group <- case_when(
  endo$cell_type %in% HEALTHY_EC_TYPES ~ "Healthy_EC",
  endo$cell_type %in% TUMOR_EC_TYPES   ~ "Tumor_EC",
  TRUE ~ NA_character_
)

endo.pb <- subset(endo, cells = rownames(endo@meta.data)[!is.na(endo$pb_group)])

# ── Tumor EC module score ──────────────────────────────────────────────────────
endo.pb <- AddModuleScore(
  endo.pb,
  features = list(TUMOR_EC_MODULE),
  name     = "TumorEC_module",
  assay    = "SCT"
)

# ── Plot 1: Module score ───────────────────────────────────────────────────────
p1 <- FeaturePlot(endo.pb, features = "TumorEC_module1") +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("")

ggsave(file.path(out_dir, "umap_tumor_ec_module.pdf"), p1, width = 6, height = 5)

# ── Plot 2: EC group assignment ────────────────────────────────────────────────
p2 <- DimPlot(endo.pb, group.by = "pb_group") +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("") +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave(file.path(out_dir, "umap_ec_groups.pdf"), p2, width = 6, height = 5)

cat("Saved UMAP plots ->", out_dir, "\n")