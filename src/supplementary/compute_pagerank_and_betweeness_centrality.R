# supplementary/compute_pagerank_and_betweeness_centrality.R
#
# Computes Reverse PageRank, Undirected Betweenness Centrality, and 
# plots Regulatory Influence vs. Centrality metrics (Figure S6).
#
# Generates:
#   - Panel A: Regulatory Influence vs. log10(Undirected Betweenness)
#   - Panel B: Regulatory Influence vs. log10(Reverse PageRank)
#
# Input:  data/eRegulon_direct.tsv
#         data/eRegulons_extended.tsv
#         figures/figure2/tf_regulatory_influence_scores.csv (or data/tf_scores.csv)
# Output: supplementary/centrality_results/tf_centrality_metrics.csv
#         supplementary/centrality_results/figure_S6A_influence_vs_betweenness.pdf
#         supplementary/centrality_results/figure_S6B_influence_vs_pagerank.pdf

library(igraph)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(here)

set.seed(101)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir  <- here::here()
data_dir  <- file.path(root_dir, "data")
fig2_dir  <- file.path(root_dir, "figures", "figure2")
out_dir   <- file.path(root_dir, "supplementary", "centrality_results")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load eRegulons & Construct Graph ───────────────────────────────────────────
cat("Loading SCENIC+ eRegulon tables...\n")
direct_path   <- file.path(data_dir, "eRegulon_direct.tsv")
extended_path <- file.path(data_dir, "eRegulons_extended.tsv")

if (!file.exists(extended_path)) {
  extended_path <- file.path(data_dir, "eRegulon_extended.tsv")
}

direct_df   <- read.delim(direct_path, stringsAsFactors = FALSE)
extended_df <- read.delim(extended_path, stringsAsFactors = FALSE)

# Combine unique TF -> Gene edge pairs with importance_x_abs_rho weights
edges_df <- bind_rows(
  direct_df   %>% select(TF, Gene, importance_x_abs_rho),
  extended_df %>% select(TF, Gene, importance_x_abs_rho)
) %>%
  distinct() %>%
  filter(!is.na(importance_x_abs_rho))

# List of all unique TFs
tfs_all <- unique(edges_df$TF)

# Construct directed igraph graph
g <- graph_from_data_frame(edges_df, directed = TRUE)

cat("GRN constructed with", vcount(g), "nodes and", ecount(g), "edges.\n")

# ── 1. Reverse PageRank Calculation ────────────────────────────────────────────
cat("Computing Reverse PageRank...\n")

g_rev <- reverse_edges(g)
w_pr  <- E(g_rev)$importance_x_abs_rho

pr_res <- page_rank(
  g_rev,
  directed = TRUE,
  weights  = w_pr
)

pr_df <- data.frame(
  gene     = names(pr_res$vector),
  pagerank = as.numeric(pr_res$vector),
  stringsAsFactors = FALSE
)

# ── 2. Undirected Betweenness Centrality Calculation ───────────────────────────
cat("Computing Undirected Betweenness Centrality...\n")

eps   <- 1e-12
w_btw <- E(g)$importance_x_abs_rho
E(g)$dist <- 1 / (w_btw + eps)

btw_res <- betweenness(
  g,
  directed   = FALSE,
  weights    = E(g)$dist,
  normalized = TRUE
)

btw_df <- data.frame(
  gene        = names(btw_res),
  betweenness = as.numeric(btw_res),
  stringsAsFactors = FALSE
)

# ── Load Regulatory Influence Scores ───────────────────────────────────────────
infl_file <- file.path(fig2_dir, "tf_regulatory_influence_scores.csv")
if (!file.exists(infl_file)) {
  infl_file <- file.path(data_dir, "tf_regulatory_influence_scores.csv")
}

cat("Loading Regulatory Influence scores from:", infl_file, "\n")
infl_df <- read.csv(infl_file, stringsAsFactors = FALSE)
if ("from" %in% colnames(infl_df)) {
  infl_df <- infl_df %>% rename(gene = from)
}

# ── Merge All Metrics ──────────────────────────────────────────────────────────
tf_metrics <- infl_df %>%
  left_join(pr_df,  by = "gene") %>%
  left_join(btw_df, by = "gene") %>%
  filter(!is.na(pagerank), !is.na(betweenness))

out_csv <- file.path(out_dir, "tf_centrality_metrics.csv")
write.csv(tf_metrics, out_csv, row.names = FALSE)
cat("Saved centrality metrics table ->", out_csv, "\n")

# ── Panel A: Regulatory Influence vs. Undirected Betweenness ──────────────────
cat("Generating Panel A scatterplot (Betweenness)...\n")

cut_inf <- median(abs(tf_metrics$final_score), na.rm = TRUE)
cut_btw <- median(abs(tf_metrics$betweenness), na.rm = TRUE)

tf_plot_btw <- tf_metrics %>%
  mutate(high_both = (final_score >= cut_inf) & (betweenness >= cut_btw))

p_btw <- ggplot(tf_plot_btw, aes(x = final_score, y = log10(betweenness + 1e-6))) +
  geom_point(aes(alpha = high_both), size = 2, color = "grey35") +
  geom_point(data = subset(tf_plot_btw, high_both), size = 2.3, color = "#7B4CCB") + # Purple highlight
  geom_text_repel(
    data          = subset(tf_plot_btw, high_both),
    aes(label     = gene),
    size          = 2.5,
    max.overlaps  = Inf,
    box.padding   = 0.2,
    point.padding = 0.1,
    segment.size  = 0.2
  ) +
  geom_vline(xintercept = cut_inf, linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = log10(cut_btw + 1e-6), linetype = "dashed", linewidth = 0.3) +
  scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 0.9), guide = "none") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title = "Betweenness (Undirected)",
    x     = "Regulatory Influence",
    y     = "log10(Undirected Betweenness)"
  )

ggsave(file.path(out_dir, "figure_S6A_influence_vs_betweenness.pdf"), p_btw, width = 6, height = 5)
ggsave(file.path(out_dir, "figure_S6A_influence_vs_betweenness.png"), p_btw, width = 6, height = 5, dpi = 300)

# ── Panel B: Regulatory Influence vs. Reverse PageRank ────────────────────────
cat("Generating Panel B scatterplot (Reverse PageRank)...\n")

cut_pr <- median(abs(tf_metrics$pagerank), na.rm = TRUE)

tf_plot_pr <- tf_metrics %>%
  mutate(high_both = (final_score >= cut_inf) & (pagerank >= cut_pr))

p_pr <- ggplot(tf_plot_pr, aes(x = final_score, y = log10(pagerank + 1e-6))) +
  geom_point(aes(alpha = high_both), size = 2, color = "grey35") +
  geom_point(data = subset(tf_plot_pr, high_both), size = 2.3, color = "blue") + # Blue highlight
  geom_text_repel(
    data          = subset(tf_plot_pr, high_both),
    aes(label     = gene),
    size          = 2.5,
    max.overlaps  = Inf,
    box.padding   = 0.2,
    point.padding = 0.1,
    segment.size  = 0.2
  ) +
  geom_vline(xintercept = cut_inf, linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = log10(cut_pr + 1e-6), linetype = "dashed", linewidth = 0.3) +
  scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 0.9), guide = "none") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title = "Reverse Pagerank",
    x     = "Regulatory Influence",
    y     = "log10(Reverse Pagerank)"
  )

ggsave(file.path(out_dir, "figure_S6B_influence_vs_pagerank.pdf"), p_pr, width = 6, height = 5)
ggsave(file.path(out_dir, "figure_S6B_influence_vs_pagerank.png"), p_pr, width = 6, height = 5, dpi = 300)

cat("Figure S6 reproduction complete!\n")
