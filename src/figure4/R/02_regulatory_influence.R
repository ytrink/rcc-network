# results/figure4/02_regulatory_influence.R
#
# Computes and plots TF regulatory influence scores for the
# normal-to-tumor endothelial cell transition (Healthy EC → Tumor EC).
#
# Input:  data/zvirblyte/zvir_endo.rds
#         data/scplus/scplus_direct_extended.graphml
# Output: results/figure4/tf_regulatory_influence_scores.csv
#         results/figure4/regulatory_influence_barplot.pdf

library(Seurat)
library(DESeq2)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(here)

set.seed(101)

# ── Paths ──────────────────────────────────────────────────────────────────────
source(here::here("src", "utils", "paths.R"))
root_dir    <- ROOT_DIR
data_dir    <- file.path(root_dir, "data")
network_dir <- file.path(root_dir, "data", "scplus")
out_dir     <- results_dir("figure4")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Parameters ─────────────────────────────────────────────────────────────────
HEALTHY_EC_TYPES <- c("AVR", "DVR", "Glomerular endothelium")
TUMOR_EC_TYPES   <- c(
  "Tumor vasculature 1", "Tumor vasculature 2",
  "Tumor vasculature 3", "Tumor vasculature 4",
  "Tumor AVR-like vasculature"
)

PADJ_CUTOFF  <- 0.05
LFC_CUTOFF   <- 0.5   # same rule as the CD8 comparison; gives identical scores to the original run (padj only)
N_LABEL_TOP  <- 5
N_LABEL_BOT  <- 5

# ── Load data ──────────────────────────────────────────────────────────────────
endo <- readRDS(file.path(data_dir, "zvirblyte", "zvir_endo.rds"))
DefaultAssay(endo) <- "RNA"

# ── Group assignment ───────────────────────────────────────────────────────────
endo$pb_group <- case_when(
  endo$cell_type %in% HEALTHY_EC_TYPES ~ "Healthy_EC",
  endo$cell_type %in% TUMOR_EC_TYPES   ~ "Tumor_EC",
  TRUE ~ NA_character_
)

endo.pb <- subset(endo, cells = rownames(endo@meta.data)[!is.na(endo$pb_group)])

# ── Pseudobulk DE ─────────────────────────────────────────────────────────────
pb <- AggregateExpression(
  endo.pb,
  assays        = "RNA",
  slot          = "counts",
  group.by      = c("patient", "pb_group"),
  return.seurat = TRUE
)

# Note: AggregateExpression() replaces "_" with "-" in group names (Tumor_EC -> Tumor-EC)
Idents(pb) <- pb$pb_group

bulk_ec_de <- FindMarkers(
  object   = pb,
  ident.1  = "Tumor-EC",
  ident.2  = "Healthy-EC",
  test.use = "DESeq2"
)

bulk_ec_de$gene      <- rownames(bulk_ec_de)
bulk_ec_de$direction <- ifelse(
  bulk_ec_de$avg_log2FC > 0, "Tumor_EC_up", "Healthy_EC_up"
)

# DEG lists used to train the perturbation regressors (figure4/python/01_get_regressor.py).
# As in the original run these use padj < 0.05 only (78 TEC-up / 81 EC-up genes).
degs_padj <- subset(bulk_ec_de, p_val_adj < PADJ_CUTOFF)
write.csv(degs_padj$gene[degs_padj$direction == "Tumor_EC_up"],   file.path(out_dir, "degs_tec.csv"))
write.csv(degs_padj$gene[degs_padj$direction == "Healthy_EC_up"], file.path(out_dir, "degs_hec.csv"))

# DEGs for the influence score
sig_degs <- subset(bulk_ec_de, p_val_adj < PADJ_CUTOFF & abs(avg_log2FC) > LFC_CUTOFF)
deg_dir  <- sig_degs[, c("gene", "direction")] %>%
  mutate(dir_num = ifelse(direction == "Tumor_EC_up", 1, -1))

# ── Load GRN and extract TF edges ─────────────────────────────────────────────
g         <- read_graph(file.path(network_dir, "scplus_direct_extended.graphml"), format = "graphml")
node_type <- setNames(V(g)$type, V(g)$name)

edges <- igraph::as_data_frame(g, what = "edges") %>%
  mutate(
    from_type = node_type[from],
    to_type   = node_type[to]
  ) %>%
  filter(from_type == "TF") %>%
  left_join(deg_dir, by = c("to" = "gene")) %>%
  filter(!is.na(direction))

# ── TF average expression ──────────────────────────────────────────────────────
DefaultAssay(pb) <- "RNA"
pb <- NormalizeData(pb, verbose = FALSE)

avg_expr <- AverageExpression(
  pb,
  assays   = "RNA",
  slot     = "data",
  group.by = "pb_group"
)$RNA %>%
  as.data.frame() %>%
  mutate(gene = rownames(.))

tf_expr <- avg_expr %>%
  filter(gene %in% edges$from) %>%
  mutate(tf_avg_expr = (`Healthy-EC` + `Tumor-EC`) / 2) %>%
  select(gene, tf_avg_expr)

# ── Regulatory influence scores ────────────────────────────────────────────────
tf_scores <- edges %>%
  left_join(tf_expr, by = c("from" = "gene")) %>%
  mutate(weighted_edge_score = rho_TF2G * dir_num * tf_avg_expr) %>%
  group_by(from) %>%
  summarise(
    final_score = sum(weighted_edge_score, na.rm = TRUE),
    n_targets   = n(),
    .groups     = "drop"
  ) %>%
  arrange(final_score)

write.csv(tf_scores, file.path(out_dir, "tf_regulatory_influence_scores.csv"), row.names = FALSE)

# ── Fig. 4C input: direct MEF2C targets up-regulated in TECs (ToppGene) ────────
# All MEF2C out-neighbours in the GRN that are TEC-up DEGs (padj < 0.05, log2FC > 0); 33 genes, all activating.
# The GO enrichment itself was run on the ToppGene web server (ToppFun, accessed 31 Aug 2025; Table S2).
mef2c_targets <- names(neighbors(g, "MEF2C", mode = "out"))
mef2c_tec     <- intersect(mef2c_targets, degs_padj$gene[degs_padj$direction == "Tumor_EC_up"])
write.csv(mef2c_tec, file.path(out_dir, "mef2c_tec_degs.csv"))

# ── Plot ───────────────────────────────────────────────────────────────────────
tf_scores_plot <- tf_scores %>%
  mutate(
    rank = row_number(),
    side = ifelse(final_score >= 0, "pos", "neg")
  )

labels_df <- bind_rows(
  slice_max(tf_scores_plot, final_score, n = N_LABEL_TOP, with_ties = FALSE),
  slice_min(tf_scores_plot, final_score, n = N_LABEL_BOT, with_ties = FALSE)
)

p <- ggplot(tf_scores_plot, aes(x = rank, y = final_score, fill = final_score)) +
  geom_col(width = 1) +
  scale_fill_gradient2(
    low      = "#5D3FD3",
    mid      = "grey90",
    high     = "#E67E22",
    midpoint = 0
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_text_repel(
    data          = labels_df,
    aes(label     = from),
    size          = 3.5,
    fontface      = "bold.italic",
    box.padding   = 0.8,
    point.padding = 0.5,
    min.segment.length = 0,
    segment.color = "grey50",
    segment.size  = 0.3,
    nudge_x       = ifelse(labels_df$side == "pos", -20, 20),
    direction     = "y"
  ) +
  labs(
    x = "Transcription Factors Ranked by Influence",
    y = "Regulatory Influence Score"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.title      = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "regulatory_influence_barplot.pdf"), p, width = 6, height = 5)
cat("Saved ->", file.path(out_dir, "regulatory_influence_barplot.pdf"), "\n")