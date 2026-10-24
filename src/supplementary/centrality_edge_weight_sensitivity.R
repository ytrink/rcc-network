# src/supplementary/centrality_edge_weight_sensitivity.R
#
# Robustness of the TF centrality results (Figures S6, S10, S13) to the choice of edge weight.
#
# The GRN graphml has one edge per TF->TG pair (networkx DiGraph); when several regulatory elements (REs) link a TF to
# a TG, the edge carries the importance_x_abs_rho of one of them. This script recomputes reverse PageRank and
# undirected betweenness in the same cell-type subnetworks as compute_pagerank_and_betweeness_centrality.R with three
# weightings per TF->TG edge:
#   w_last : as in the graphml (published)
#   w_max  : maximum importance_x_abs_rho over the supporting REs
#   w_sum  : sum over the supporting REs
# and reports, for each, the TFs above the median on both regulatory influence and centrality.
#
# Result (Methods): the main TFs stay highlighted under all three weightings (SPI1; EOMES, IKZF3; MEF2C, ETS2, ELK3,
# ERG). MAFB stays high in betweenness, but its PageRank passes the median only with the published weighting.
#
# Input:  data/scplus/{scplus_direct_extended.graphml, eRegulon_direct.tsv, eRegulons_extended.tsv}
#         genes expressed in >=10% of cells and influence scores, as in compute_pagerank_and_betweeness_centrality.R
# Output: results/supplementary/centrality/centrality_edge_weight_sensitivity.csv

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(Seurat)
  library(here)
})

source(here::here("src", "utils", "paths.R"))
out_dir <- results_dir("supplementary", "centrality")

g <- read_graph(data_path("scplus", "scplus_direct_extended.graphml"), format = "graphml")

# ── Region-level weights per TF->TG pair ───────────────────────────────────────
cols <- c("Region", "Gene", "TF", "importance_x_abs_rho")
region_links <- distinct(bind_rows(read.delim(data_path("scplus", "eRegulon_direct.tsv"))[, cols],
                                   read.delim(data_path("scplus", "eRegulons_extended.tsv"))[, cols]))
agg <- region_links %>%
  group_by(TF, Gene) %>%
  summarise(w_max = max(importance_x_abs_rho), w_sum = sum(importance_x_abs_rho), n_regions = n(), .groups = "drop")
cat("TF-TG pairs:", nrow(agg), "| pairs with >1 RE:", sum(agg$n_regions > 1), "\n")

edges <- igraph::as_data_frame(g, what = "edges")
m <- match(paste(edges$from, edges$to), paste(agg$TF, agg$Gene))
stopifnot(!any(is.na(m)))
E(g)$w_last <- E(g)$importance_x_abs_rho
E(g)$w_max  <- agg$w_max[m]
E(g)$w_sum  <- agg$w_sum[m]
cat("edges whose graphml weight equals the max RE weight:",
    round(100 * mean(abs(E(g)$w_last - E(g)$w_max) < 1e-12), 1), "%\n")

# ── Subnetworks, as in compute_pagerank_and_betweeness_centrality.R (Table S3) ──
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
  endo     = file.path(RESULTS_DIR, "figure4", "zvir_endo_norm.rds")
)
INFLUENCE <- list(
  mac_mono = file.path(RESULTS_DIR, "figure2", "tf_regulatory_influence_scores.csv"),
  CD8      = file.path(RESULTS_DIR, "figure3", "tf_regulatory_influence_scores.csv"),
  endo     = file.path(RESULTS_DIR, "figure4", "tf_regulatory_influence_scores.csv")
)

expressed_genes <- function(cell_type) {
  precomputed <- data_path("li", paste0(cell_type, "_expressed_genes_10perc.csv"))
  if (file.exists(precomputed)) return(read.csv(precomputed, row.names = 1)$x)
  obj    <- readRDS(EXPRESSION[[cell_type]])
  counts <- GetAssayData(JoinLayers(obj, assay = "RNA"), assay = "RNA", layer = "counts")
  names(which(Matrix::rowMeans(counts > 0) >= 0.10))
}

results <- list()
for (cell_type in names(TF_LISTS)) {
  tfs     <- TF_LISTS[[cell_type]]
  targets <- setdiff(V(g)$name[unique(unlist(neighborhood(g, order = 1, nodes = tfs, mode = "out")))], tfs)
  keep    <- intersect(unique(c(tfs, targets)), expressed_genes(cell_type))
  sub_g   <- induced_subgraph(g, vids = V(g)[name %in% keep])
  infl    <- read.csv(INFLUENCE[[cell_type]]) %>% rename(gene = from)
  cat("\n==", cell_type, ":", vcount(sub_g), "nodes,", ecount(sub_g), "edges\n")

  for (w in c("w_last", "w_max", "w_sum")) {
    weights <- edge_attr(sub_g, w)
    pr  <- page_rank(reverse_edges(sub_g), directed = TRUE, weights = weights)$vector
    btw <- betweenness(sub_g, directed = FALSE, weights = 1 / (weights + 1e-12), normalized = TRUE)
    tm  <- infl %>%
      inner_join(data.frame(gene = names(pr), pagerank = as.numeric(pr), betweenness = as.numeric(btw[names(pr)])),
                 by = "gene")
    cut_inf <- median(abs(tm$final_score))
    tm$high_betweenness <- tm$final_score >= cut_inf & tm$betweenness >= median(abs(tm$betweenness))
    tm$high_pagerank    <- tm$final_score >= cut_inf & tm$pagerank >= median(abs(tm$pagerank))
    cat(sprintf("  [%s] high influence + betweenness: %s\n           high influence + PageRank:    %s\n", w,
                paste(sort(tm$gene[tm$high_betweenness]), collapse = ", "),
                paste(sort(tm$gene[tm$high_pagerank]), collapse = ", ")))
    tm$cell_type <- cell_type
    tm$weight    <- w
    results[[paste(cell_type, w)]] <- tm
  }
}
write.csv(bind_rows(results), file.path(out_dir, "centrality_edge_weight_sensitivity.csv"), row.names = FALSE)
cat("Saved ->", out_dir, "\n")
