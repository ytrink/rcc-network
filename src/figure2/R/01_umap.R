# Figure 2 - Panel A: Macrophage/Monocyte UMAP plots
# Generates: (1) macro group UMAP, (2) TREM2 expression UMAP

library(Seurat)
library(dplyr)
library(ggplot2)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir  <- "https://github.com/ytrink/rcc-network"
data_dir  <- file.path(root_dir, "data")
umap_dir  <- file.path(root_dir, "preprocessing", "monocle", "mac_mono")
out_dir   <- file.path(root_dir, "figures", "figure2")

# ── Load data ──────────────────────────────────────────────────────────────────
rna.li <- readRDS(file.path(data_dir, "Li_mac_mono.rds"))
colnames(rna.li) <- paste0("LI_", rna.li$patient, "_", colnames(rna.li))
DefaultAssay(rna.li) <- "RNA"

X_umap <- read.csv(file.path(umap_dir, "cds_mac_mono_umap_sampled30k.csv"), row.names = 1)

# ── Subset to common cells ─────────────────────────────────────────────────────
common_cells  <- intersect(rownames(X_umap), colnames(rna.li))
rna.li.umap   <- rna.li[, common_cells]
X_umap        <- X_umap[colnames(rna.li.umap), ]

stopifnot(all(rownames(X_umap) == colnames(rna.li.umap)))

# ── Embed Monocle UMAP into Seurat object ──────────────────────────────────────
rna.li.umap[["umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(X_umap),
  key        = "UMAP_",
  assay      = DefaultAssay(rna.li.umap)
)
rna.li.umap <- NormalizeData(rna.li.umap)

# ── Cell type grouping ─────────────────────────────────────────────────────────
rna.li.umap$macro_group <- case_when(
  grepl("Mono",   rna.li.umap$annotation) ~ "Monocytes",
  grepl("TR Mac", rna.li.umap$annotation) ~ "Tissue Resident Macrophages",
  grepl("TAM",    rna.li.umap$annotation) ~ "TAMs",
  TRUE ~ NA_character_
)

rna.li.umap$macro_group <- factor(
  rna.li.umap$macro_group,
  levels = c("Monocytes", "Tissue Resident Macrophages", "TAMs")
)

# ── Plot 1: Macro group UMAP ───────────────────────────────────────────────────
p1 <- DimPlot(rna.li.umap, group.by = "macro_group") +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("") +
  theme(
    legend.position = c(0.62, 0.85),
    legend.text     = element_text(size = 14)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave(file.path(out_dir, "umap_macro_groups.pdf"), p1, width = 6, height = 5)

# ── Plot 2: TREM2 expression UMAP ─────────────────────────────────────────────
p2 <- FeaturePlot(rna.li.umap, features = "TREM2") +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("")

ggsave(file.path(out_dir, "umap_TREM2.pdf"), p2, width = 6, height = 5)