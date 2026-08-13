# figures/figure3/R/01_umap.R
#
# Generates UMAP plots for CD8+ T cell exhaustion analysis (Figure 3).
# Plots: (1) exhaustion score, (2) Ex/Non-ex group assignment
#
# Input:  data/li/Li_CD8T_normalized.rds
#         data/monocle/cd8/cds_tcells_umap.csv
# Output: figures/figure3/R/umap_exhaustion_score.pdf
#         figures/figure3/R/umap_ex_groups.pdf

library(Seurat)
library(dplyr)
library(ggplot2)
library(here)

set.seed(101)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir <- here::here()
data_dir <- file.path(root_dir, "data")
out_dir  <- file.path(root_dir, "figures", "figure3", "R")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Parameters ─────────────────────────────────────────────────────────────────
EXH_GENES <- c("TIGIT", "LAG3")
MIN_CELLS <- 50

# ── Load data ──────────────────────────────────────────────────────────────────
rna.li <- readRDS(file.path(data_dir, "li", "Li_CD8T_normalized.rds"))
colnames(rna.li) <- paste0("LI_", rna.li$patient, "_", colnames(rna.li))
DefaultAssay(rna.li) <- "RNA"
rna.li <- NormalizeData(rna.li, verbose = FALSE)
rna.li <- JoinLayers(rna.li)

# ── Exhaustion score and group assignment ──────────────────────────────────────
exh_genes <- intersect(EXH_GENES, rownames(rna.li))
rna.li    <- AddModuleScore(rna.li, features = list(exh_genes), name = "ExhScore")
score     <- rna.li$ExhScore1
med_val   <- median(score)

rna.li$ex_group <- ifelse(score >= med_val, "Ex", "Non_ex")

# ── Filter to paired patients with sufficient cells ────────────────────────────
tab             <- table(rna.li$patient, rna.li$ex_group)
paired_patients <- rownames(tab)[rowSums(tab > 0) == 2]
rna.li          <- subset(rna.li, subset = patient %in% paired_patients)

md           <- rna.li@meta.data
md$sample_id <- paste(md$patient, md$ex_group, sep = "__")
keep_cells   <- rownames(md)[md$sample_id %in% names(which(table(md$sample_id) >= MIN_CELLS))]
rna.li       <- subset(rna.li, cells = keep_cells)

# ── Load Monocle UMAP and embed ────────────────────────────────────────────────
X_umap      <- read.csv(file.path(data_dir, "monocle", "cd8", "cds_tcells_umap.csv"), row.names = 1)
common_cells <- intersect(rownames(X_umap), colnames(rna.li))
rna.li.umap  <- rna.li[, common_cells]
X_umap       <- X_umap[colnames(rna.li.umap), ]

stopifnot(all(rownames(X_umap) == colnames(rna.li.umap)))

rna.li.umap[["umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(X_umap),
  key        = "UMAP_",
  assay      = DefaultAssay(rna.li.umap)
)

# ── Plot 1: Exhaustion score ───────────────────────────────────────────────────
p1 <- FeaturePlot(rna.li.umap, features = "ExhScore1") +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("")

ggsave(file.path(out_dir, "umap_exhaustion_score.pdf"), p1, width = 6, height = 5)

# ── Plot 2: Ex / Non-ex groups ────────────────────────────────────────────────
p2 <- DimPlot(rna.li.umap, group.by = "ex_group") +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("")

ggsave(file.path(out_dir, "umap_ex_groups.pdf"), p2, width = 6, height = 5)

cat("Saved UMAP plots ->", out_dir, "\n")