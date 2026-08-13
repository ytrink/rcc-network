# figures/figure4/R/02_regulatory_influence.R
#
# Computes and plots TF regulatory influence scores for the
# normal-to-tumor endothelial cell transition (Healthy EC → Tumor EC).
#
# Input:  data/zvirblyte/Zvir_endo.rds
#         data/scplus/scplus_direct_extended.graphml
# Output: figures/figure4/R/tf_regulatory_influence_scores.csv
#         figures/figure4/R/regulatory_influence_barplot.pdf

library(Seurat)
library(DESeq2)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(here)

set.seed(101)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir    <- here::here()
data_dir    <- file.path(root_dir, "data")
network_dir <- file.path(root_dir, "data", "scplus")
out_dir     <- file.path(root_dir, "figures", "figure4", "R")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Parameters ─────────────────────────────────────────────────────────────────
HEALTHY_EC_TYPES <- c("AVR", "DVR", "Glomerular endothelium")
TUMOR_EC_TYPES   <- c(
  "Tumor vasculature 1", "Tumor vasculature 2",
  "Tumor vasculature 3", "Tumor vasculature 4",
  "Tumor AVR-like vasculature"
)

PADJ_CUTOFF  <- 0.05
N_LABEL_TOP  <- 5
N_LABEL_BOT  <- 5

# ── Load data ──────────────────────────────────────────────────────────────────
endo <- readRDS(file.path(data_dir, "zvirblyte", "Zvir_endo.rds"))
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

Idents(pb) <- pb$pb_group

bulk_ec_de <- FindMarkers(
  object   = pb,
  ident.1  = "Tumor_EC",
  ident.2  = "Healthy_EC",
  test.use = "DESeq2"
)

bulk_ec_de$gene      <- rownames(bulk_ec_de)
bulk_ec_de$direction <- ifelse(
  bulk_ec_de$avg_log2FC > 0, "Tumor_EC_up", "Healthy_EC_up"
)

sig_degs <- subset(bulk_ec_de, p_val_adj < PADJ_CUTOFF)
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
  mutate(tf_avg_expr = (Healthy_EC + Tumor_EC) / 2) %>%
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