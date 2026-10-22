# src/supplementary/compute_pagerank_and_betweeness_centrality.R
#
# TF centrality in cell-type-specific regulatory subnetworks (Figures S6, S10, S13; Table S3).
# Port of plot_regulon_communities_LI/{mac_mono,CD8,endo}/compute_other_metrics.R
#
# For each cell type:
#   1. TFs = hand-picked TFs whose regulons are active in that cell type in the AUCell heatmap (Fig. 1C;
#      Table S3). These largely overlap with the Leiden communities of Table S1.
#   2. Subnetwork = these TFs + their direct targets in the GRN, restricted to genes (TFs and TGs) detected
#      (counts > 0) in >= 10% of cells of that cell type; induced subgraph (so TFs that are targets of the
#      listed TFs are included and also scored).
#   3. Edge weight = SCENIC+ importance_x_abs_rho. The GRN graphml has one edge per TF->TG pair (built with a
#      networkx DiGraph, so the weight is that of one supporting region-to-gene link).
#   4. Reverse PageRank (edges reversed, weighted) and undirected betweenness (distance = 1 / weight).
#   5. Scatter of regulatory influence vs centrality; TFs above the median on both axes are labelled.
#
# Input:  data/scplus/scplus_direct_extended.graphml
#         expression: data/li/Li_mac_mono.rds, data/li/Li_CD8T_normalized.rds, results/figure4/zvir_endo_norm.rds
#           (or precomputed gene lists data/li/{mac_mono,CD8}_expressed_genes_10perc.csv if present)
#         influence scores: results/figure{2,3,4}/tf_regulatory_influence_scores.csv
# Output: results/supplementary/centrality/<cell_type>_tf_centrality_metrics.csv and scatter plots

library(Seurat)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(here)

source(here::here("src", "utils", "paths.R"))
out_dir <- results_dir("supplementary", "centrality")

# ── Table S3: TFs defining each subnetwork ─────────────────────────────────────
TF_LISTS <- list(
  mac_mono = c("TCF7L2", "CEBPD", "MITF", "ATF3", "FOS", "CREB5", "HIF1A", "MEF2A", "ZFHX3", "ETV6",
               "USF2", "JDP2", "YBX1", "KLF4", "IRF8", "IRF5", "SPI1", "CEBPB", "CEBPA", "MAFB",
               "POU2F2", "REL", "HMGA1", "NR4A3", "MAF"),
  CD8      = c("RAD21", "ID2", "IRF2", "ETS1", "PRDM1", "ELF1", "GATA3", "STAT4", "IKZF2", "TBX21",
               "IKZF3", "EOMES", "BCL11B", "BATF", "CREM", "RUNX3", "IKZF1", "RUNX1", "MAZ", "HMGXB3",
               "BHLHE40"),
  endo     = c("TBX3", "EBF1", "PDLIM5", "ZBTB20", "RARB", "KLF9", "FOXF1", "ETS2", "EBF3", "TCF7L1",
               "ERG", "SOX7", "GATA2", "SOX17", "ELK3", "TAL1", "PBX1", "ZEB1", "SOX13", "HLX", "MEF2C",
               "FLI1", "LMO2", "HMBOX1", "TCF4", "NFAT5", "EBF2", "PRRX1", "HEYL", "NFATC4", "NFIC",
               "HES1", "HOXD8", "NFIA", "NFIB", "SOX4", "FOXC1", "SMAD9", "NR2F2", "NR2F1", "NFIX",
               "TEAD1", "TEAD4", "HOXB7")
)

EXPRESSION <- list(
  mac_mono = data_path("li", "Li_mac_mono.rds"),
  CD8      = data_path("li", "Li_CD8T_normalized.rds"),
  endo     = file.path(RESULTS_DIR, "figure4", "zvir_endo_norm.rds")     # figure4/R/01_umap.R
)
INFLUENCE <- list(
  # Note: the published Fig. S6 used an earlier (non-pseudobulk) version of the mac/mono influence scores;
  # with the pseudobulk scores of Fig. 2B, SPI1 and MAFB are highlighted in both panels as well.
  mac_mono = file.path(RESULTS_DIR, "figure2", "tf_regulatory_influence_scores.csv"),
  CD8      = file.path(RESULTS_DIR, "figure3", "tf_regulatory_influence_scores.csv"),
  endo     = file.path(RESULTS_DIR, "figure4", "tf_regulatory_influence_scores.csv")
)
FIGURE <- c(mac_mono = "S6", CD8 = "S10", endo = "S13")

# Genes detected (counts > 0) in >= 10% of cells
expressed_genes <- function(cell_type) {
  precomputed <- data_path("li", paste0(cell_type, "_expressed_genes_10perc.csv"))
  if (file.exists(precomputed)) return(read.csv(precomputed, row.names = 1)$x)
  obj    <- readRDS(EXPRESSION[[cell_type]])
  counts <- GetAssayData(JoinLayers(obj, assay = "RNA"), assay = "RNA", layer = "counts")
  names(which(Matrix::rowMeans(counts > 0) >= 0.10))
}

scatter <- function(df, metric, label, colour) {
  cut_inf    <- median(abs(df$final_score), na.rm = TRUE)
  cut_metric <- median(abs(df[[metric]]),   na.rm = TRUE)
  df$high_both <- df$final_score >= cut_inf & df[[metric]] >= cut_metric
  ggplot(df, aes(x = final_score, y = log10(.data[[metric]] + 1e-6))) +
    geom_point(aes(alpha = high_both), size = 2, color = "grey35") +
    geom_point(data = subset(df, high_both), size = 2.3, color = colour) +
    geom_text_repel(data = subset(df, high_both), aes(label = gene), size = 2.5, max.overlaps = Inf,
                    box.padding = 0.2, point.padding = 0.1, segment.size = 0.2) +
    geom_vline(xintercept = cut_inf, linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = log10(cut_metric + 1e-6), linetype = "dashed", linewidth = 0.3) +
    scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 0.9), guide = "none") +
    theme_minimal(base_size = 12) +
    labs(x = "Regulatory Influence", y = paste0("log10(", label, ")"), title = "")
}

g <- read_graph(data_path("scplus", "scplus_direct_extended.graphml"), format = "graphml")
tfs_all <- V(g)$name[V(g)$type == "TF"]

for (cell_type in names(TF_LISTS)) {
  message("== ", cell_type)
  tfs     <- TF_LISTS[[cell_type]]
  targets <- setdiff(V(g)$name[unique(unlist(neighborhood(g, order = 1, nodes = tfs, mode = "out")))], tfs)
  keep    <- intersect(unique(c(tfs, targets)), expressed_genes(cell_type))
  sub_g   <- induced_subgraph(g, vids = V(g)[name %in% keep])
  message("subnetwork: ", vcount(sub_g), " nodes, ", ecount(sub_g), " edges; listed TFs kept: ",
          sum(tfs %in% keep), "/", length(tfs))

  # Reverse PageRank
  pr <- page_rank(reverse_edges(sub_g), directed = TRUE, weights = E(sub_g)$importance_x_abs_rho)$vector
  # Undirected betweenness, distance = 1 / weight
  btw <- betweenness(sub_g, directed = FALSE, weights = 1 / (E(sub_g)$importance_x_abs_rho + 1e-12),
                     normalized = TRUE)

  metrics <- data.frame(gene = names(pr), pagerank = as.numeric(pr), betweenness = as.numeric(btw[names(pr)])) %>%
    filter(gene %in% tfs_all)

  tf_metrics <- read.csv(INFLUENCE[[cell_type]]) %>%
    rename(gene = from) %>%
    inner_join(metrics, by = "gene")
  write.csv(tf_metrics, file.path(out_dir, paste0(cell_type, "_tf_centrality_metrics.csv")), row.names = FALSE)

  fig <- FIGURE[[cell_type]]
  for (ext in c("pdf", "png")) {
    ggsave(file.path(out_dir, paste0("Fig", fig, "_", cell_type, "_influence_vs_betweenness.", ext)),
           scatter(tf_metrics, "betweenness", "Undirected Betweenness", "#7B4CCB"), width = 7, height = 5, dpi = 300, bg = "white")
    ggsave(file.path(out_dir, paste0("Fig", fig, "_", cell_type, "_influence_vs_pagerank.", ext)),
           scatter(tf_metrics, "pagerank", "Reverse PageRank", "blue"), width = 7, height = 5, dpi = 300, bg = "white")
  }
}
cat("Saved ->", out_dir, "\n")
