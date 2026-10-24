# src/supplementary/tnk_subnetwork_for_cytoscape.R
#
# T/NK cell gene regulatory subnetwork shown in Figure S3 (drawn in Cytoscape).
#
# TFs: the CD8+ T cell TF list of Table S3 (the same list as the centrality analysis, Figure S10).
# Network: edges from these TFs (those expressed in >=10% of CD8+ T cells) to their targets expressed in >=10% of
# CD8+ T cells. Other TFs can appear as targets only. Node attribute "listed_tf" marks the Table S3 TFs; edge colour in
# the figure follows the sign of rho_TF2G (activation / repression).
#
# Input:  data/scplus/scplus_direct_extended.graphml
#         data/li/CD8_expressed_genes_10perc.csv (or computed from data/li/Li_CD8T_normalized.rds)
# Output: results/supplementary/tnk_subnetwork/tnk_subnetwork.graphml (import into Cytoscape)

suppressPackageStartupMessages({
  library(igraph)
  library(Seurat)
  library(here)
})

source(here::here("src", "utils", "paths.R"))
out_dir <- results_dir("supplementary", "tnk_subnetwork")

TFS <- c("RAD21", "ID2", "IRF2", "ETS1", "PRDM1", "ELF1", "GATA3", "STAT4", "IKZF2", "TBX21", "IKZF3", "EOMES",
         "BCL11B", "BATF", "CREM", "RUNX3", "IKZF1", "RUNX1", "MAZ", "HMGXB3", "BHLHE40")

precomputed <- data_path("li", "CD8_expressed_genes_10perc.csv")
if (file.exists(precomputed)) {
  expressed <- read.csv(precomputed, row.names = 1)$x
} else {
  obj       <- readRDS(data_path("li", "Li_CD8T_normalized.rds"))
  counts    <- GetAssayData(JoinLayers(obj, assay = "RNA"), assay = "RNA", layer = "counts")
  expressed <- names(which(Matrix::rowMeans(counts > 0) >= 0.10))
}

g <- read_graph(data_path("scplus", "scplus_direct_extended.graphml"), format = "graphml")
tfs_expressed <- intersect(TFS, expressed)
cat("TFs expressed in >=10% of CD8+ T cells:", length(tfs_expressed), "/", length(TFS),
    "| not expressed:", paste(setdiff(TFS, tfs_expressed), collapse = ", "), "\n")

keep_edges <- E(g)[.from(V(g)[name %in% tfs_expressed]) & .to(V(g)[name %in% expressed])]
sub_g <- subgraph_from_edges(g, eids = keep_edges, delete.vertices = TRUE)
V(sub_g)$listed_tf <- V(sub_g)$name %in% TFS

cat("nodes:", vcount(sub_g), "| edges:", ecount(sub_g),
    "| activating:", sum(E(sub_g)$rho_TF2G > 0), "| repressive:", sum(E(sub_g)$rho_TF2G <= 0), "\n")
write_graph(sub_g, file.path(out_dir, "tnk_subnetwork.graphml"), format = "graphml")
cat("Saved ->", out_dir, "\n")
