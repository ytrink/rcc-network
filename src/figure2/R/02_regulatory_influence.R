# Figure 2 - Panel B: Regulatory Influence Analysis (Monocyte → TAM transition)
# Generates: regulatory influence barplot

library(Seurat)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir     <- "https://github.com/ytrink/rcc-network"
data_dir     <- file.path(root_dir, "data")
network_dir  <- file.path(root_dir, "network_analysis")
out_dir      <- file.path(root_dir, "figures", "figure2")

# ── Load data ──────────────────────────────────────────────────────────────────
rna.li <- readRDS(file.path(data_dir, "Li_mac_mono.rds"))
colnames(rna.li) <- paste0("LI_", rna.li$patient, "_", colnames(rna.li))
DefaultAssay(rna.li) <- "RNA"

# ── Cell type grouping ─────────────────────────────────────────────────────────
mono_labels <- c(
  "Classical Mono.1", "Classical Mono.2", "Classical Mono.3",
  "Classical Mono.4", "Non-classical Mono"
)
tam_labels <- c(
  "RGS+ TAM", "Pro-infla. TAM", "MHC-II TAM",
  "GPNMB+ TAM", "SPP1+ TAM", "FN1+ TAM"
)

rna.li$mono_tam_group <- case_when(
  rna.li$annotation %in% mono_labels ~ "Monocyte",
  rna.li$annotation %in% tam_labels  ~ "TAM",
  TRUE ~ "TR Macrophage"
)

rna.li.subset <- subset(rna.li, subset = mono_tam_group %in% c("Monocyte", "TAM"))

# ── Pseudobulk differential expression ────────────────────────────────────────
rna.li.pb <- AggregateExpression(
  rna.li.subset,
  assays     = "RNA",
  slot       = "counts",
  group.by   = c("patient", "mono_tam_group"),
  return.seurat = TRUE
)

cts     <- GetAssayData(rna.li.pb, slot = "counts")
coldata <- rna.li.pb@meta.data
coldata$condition <- factor(coldata$mono_tam_group, levels = c("Monocyte", "TAM"))

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(cts)),
  colData   = coldata,
  design    = ~ patient + condition
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "TAM", "Monocyte"))

# ── DE gene filtering ──────────────────────────────────────────────────────────
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df$direction <- case_when(
  res_df$padj < 0.1 & res_df$log2FoldChange >  0.5 ~ "TAM_up",
  res_df$padj < 0.1 & res_df$log2FoldChange < -0.5 ~ "Mono_up",
  TRUE ~ "NS"
)

deg_dir <- res_df %>%
  filter(direction != "NS") %>%
  select(gene, direction) %>%
  mutate(dir_num = ifelse(direction == "TAM_up", 1, -1))

# ── Load GRN and extract TF edges ─────────────────────────────────────────────
g <- read_graph(file.path(network_dir, "scplus_direct_extended.graphml"), format = "graphml")

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
DefaultAssay(rna.li.pb) <- "RNA"
rna.li.pb <- NormalizeData(rna.li.pb)

avg_expr <- AverageExpression(
  rna.li.pb,
  assays   = "RNA",
  slot     = "data",
  group.by = "mono_tam_group"
)$RNA %>%
  as.data.frame() %>%
  mutate(gene = rownames(.))

tf_expr <- avg_expr %>%
  filter(gene %in% edges$from) %>%
  mutate(tf_avg_expr = (TAM + Monocyte) / 2) %>%
  select(gene, tf_avg_expr)

# ── Regulatory influence scores ────────────────────────────────────────────────
tf_scores <- edges %>%
  left_join(tf_expr, by = c("from" = "gene")) %>%
  mutate(weighted_edge_score = rho_TF2G * dir_num * tf_avg_expr) %>%
  group_by(from) %>%
  summarise(
    final_score = sum(weighted_edge_score, na.rm = TRUE),
    n_targets   = n()
  ) %>%
  arrange(final_score)

write.csv(tf_scores, file.path(out_dir, "tf_regulatory_influence_scores.csv"), row.names = FALSE)

# ── Plot: Regulatory influence barplot ────────────────────────────────────────
tf_scores_plot <- tf_scores %>%
  mutate(
    rank = row_number(),
    side = ifelse(final_score >= 0, "pos", "neg")
  )

labels_df <- bind_rows(
  slice_max(tf_scores_plot, final_score, n = 5, with_ties = FALSE),
  slice_min(tf_scores_plot, final_score, n = 5, with_ties = FALSE)
)

p <- ggplot(tf_scores_plot, aes(x = rank, y = final_score, fill = final_score)) +
  geom_col(width = 1) +
  scale_fill_gradient2(
    low     = "green",
    mid     = "grey90",
    high    = "orange",
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