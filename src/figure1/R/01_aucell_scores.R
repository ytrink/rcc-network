# src/figure1/R/01_aucell_scores.R
#
# AUCell regulon activity scores for the Li et al. scRNA-seq data (Figure 1C).
# Port of stuff_to_run_on_the_side/AUC_SCPLUS/scplus_fragments_downsampled/auc_li_all/run_auc_all_ccRCC.R
#
#   - cells: Li et al. cells from the 10 patients with ccRCC histology (metadata_merged.csv, Histology == "ccRCC")
#   - expression: RNA, LogNormalize (scale factor 1e4)
#   - gene sets: activating regulons (positive rho_TF2G), i.e. each TF's targets in the positive-link GRN
#     (originally written to scplus_direct_extended_POS.xlsx in analysis.ipynb; built here from the graphml)
#   - AUCell scores are row-standardized (z-score per regulon) for the heatmap
#
# Input:  data/li/Li_integrated.rds              (preprocessing/scrna/li/02_processing_li.R)
#         data/li/metadata_merged.csv
#         data/scplus/scplus_direct_extended.graphml
# Output: results/figure1/aucell/auc_matrix.csv, auc_matrix_stand.csv, cell_annotations.csv

library(Seurat)
library(AUCell)
library(GSEABase)
library(igraph)
library(data.table)
library(here)

source(here::here("src", "utils", "paths.R"))
out_dir <- results_dir("figure1", "aucell")

# ── Gene sets: activating regulons ─────────────────────────────────────────────
g     <- read_graph(data_path("scplus", "scplus_direct_extended.graphml"), format = "graphml")
edges <- igraph::as_data_frame(g, what = "edges")
edges <- edges[edges$regulation == 1, ]                      # = rho_TF2G > 0
regulons  <- split(edges$to, edges$from)
gene_sets <- GeneSetCollection(lapply(names(regulons), function(tf) GeneSet(unique(regulons[[tf]]), setName = tf)))

# ── Cells: ccRCC patients only ─────────────────────────────────────────────────
data.li <- readRDS(data_path("li", "Li_integrated.rds"))
colnames(data.li) <- paste0("LI_", data.li$patient, "_", colnames(data.li))

sample_metadata <- read.csv(data_path("li", "metadata_merged.csv"))
rownames(sample_metadata) <- sample_metadata$barcode
stopifnot(all(colnames(data.li) %in% rownames(sample_metadata)))
sample_metadata <- sample_metadata[colnames(data.li), ]

data.li <- data.li[, sample_metadata$Histology == "ccRCC"]
cat("ccRCC cells:", ncol(data.li), "from", length(unique(data.li$patient)), "patients\n")

DefaultAssay(data.li) <- "RNA"
data.li <- NormalizeData(data.li, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
expr_matrix <- GetAssayData(JoinLayers(data.li, overwrite = TRUE), assay = "RNA", layer = "data")

# ── AUCell ─────────────────────────────────────────────────────────────────────
gene_sets <- subsetGeneSets(gene_sets, rownames(expr_matrix))
gene_sets <- setGeneSetNames(gene_sets, newNames = paste0(names(gene_sets), " (", nGenes(gene_sets), "g)"))

cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE, nCores = 24)
cells_AUC      <- AUCell_calcAUC(gene_sets, cells_rankings)

auc_matrix <- as.data.frame(cells_AUC@assays@data@listData$AUC)
fwrite(data.table(Gene = rownames(auc_matrix), auc_matrix), file.path(out_dir, "auc_matrix.csv"))
write.csv(data.li$annotation, file.path(out_dir, "cell_annotations.csv"))

# ── Row-standardize (z-score per regulon) ──────────────────────────────────────
auc_matrix_std <- t(scale(t(as.matrix(auc_matrix))))
write.csv(auc_matrix_std, file.path(out_dir, "auc_matrix_stand.csv"))
cat("Saved ->", out_dir, "\n")
