# figures/figure3/R/02_regulatory_influence.R
#
# Computes and plots TF regulatory influence scores for the
# CD8+ T cell exhaustion transition (Non-exhausted → Exhausted).
#
# Input:  data/li/Li_CD8T_normalized.rds
#         data/scplus/scplus_direct_extended.graphml
# Output: figures/figure3/R/tf_regulatory_influence_scores.csv
#         figures/figure3/R/regulatory_influence_barplot.pdf

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
out_dir     <- file.path(root_dir, "figures", "figure3", "R")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Parameters ─────────────────────────────────────────────────────────────────
EXH_GENES     <- c("TIGIT", "LAG3")
MEDIAN_SPLIT  <- TRUE       # split on median exhaustion score
MIN_CELLS     <- 50         # minimum cells per patient x group
LFC_CUTOFF    <- 0.5
PADJ_CUTOFF   <- 0.05
N_LABEL_TOP   <- 6
N_LABEL_BOT   <- 5

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

# ── Filter to paired patients with sufficient cells ───────────────────────────
rna.li.sub <- subset(rna.li, cells = colnames(rna.li)[!is.na(rna.li$ex_group)])

# Keep only patients with both groups
tab             <- table(rna.li.sub$patient, rna.li.sub$ex_group)
paired_patients <- rownames(tab)[rowSums(tab > 0) == 2]
rna.li.sub      <- subset(rna.li.sub, subset = patient %in% paired_patients)

# Minimum cells per patient x group
md              <- rna.li.sub@meta.data
md$sample_id    <- paste(md$patient, md$ex_group, sep = "__")
keep_samples    <- names(which(table(md$sample_id) >= MIN_CELLS))
rna.li.sub      <- subset(rna.li.sub, cells = rownames(md)[md$sample_id %in% keep_samples])

# ── Pseudobulk DE ─────────────────────────────────────────────────────────────
rna.li.pb <- AggregateExpression(
  rna.li.sub,
  assays        = "RNA",
  slot          = "counts",
  group.by      = c("patient", "ex_group"),
  return.seurat = TRUE
)

cts     <- GetAssayData(rna.li.pb, assay = "RNA", slot = "counts")
coldata <- rna.li.pb@meta.data
coldata$condition <- factor(coldata$ex_group, levels = c("Non_ex", "Ex"))

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(cts)),
  colData   = coldata,
  design    = ~ patient + condition
)
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Ex", "Non_ex"))

# ── DE gene filtering ──────────────────────────────────────────────────────────
res_df <- as.data.frame(res) %>%
  filter(!is.na(padj)) %>%
  mutate(
    gene      = rownames(.),
    direction = case_when(
      padj < PADJ_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "Ex_up",
      padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "Non_ex_up",
      TRUE ~ "NS"
    )
  )

deg_dir <- res_df %>%
  filter(direction != "NS") %>%
  select(gene, direction) %>%
  mutate(dir_num = ifelse(direction == "Ex_up", 1, -1))

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
DefaultAssay(rna.li.pb) <- "RNA"
rna.li.pb <- NormalizeData(rna.li.pb, verbose = FALSE)

avg_expr <- AverageExpression(
  rna.li.pb,
  assays   = "RNA",
  slot     = "data",
  group.by = "ex_group"
)$RNA %>%
  as.data.frame() %>%
  mutate(gene = rownames(.))

tf_expr <- avg_expr %>%
  filter(gene %in% edges$from) %>%
  mutate(tf_avg_expr = (Ex + Non_ex) / 2) %>%
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
    low      = "orange",
    mid      = "grey90",
    high     = "purple",
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
cat("Saved -> ", file.path(out_dir, "regulatory_influence_barplot.pdf"), "\n")