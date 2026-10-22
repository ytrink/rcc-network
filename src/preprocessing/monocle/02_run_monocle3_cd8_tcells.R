# src/preprocessing/monocle/02_run_monocle3_cd8_tcells.R
#
# Monocle3 embedding of Li et al. CD8+ T cells (UMAP of Fig. 3A/F) and the cell subsets used by the
# perturbation analysis (Fig. 3E/F).
# Port of stuff_to_run_on_the_side/monocle_scplus/run_monocle.R (object cds_1.rds, verified: 84,849 cells,
# 100 PCA/aligned dims) and monocle/T_NK/Li/run_monocle.R (exports).
#
# Cells: 11 CD8+ T-cell subtypes; cycling CD8+ T cells (and gamma-delta/MAIT cells) are excluded.
# Note: the original script aligns on a "Sample" column that is not part of the Li metadata; the batch
#       variable used here is "patient" (as for macrophages/monocytes).
# Note: Monocle's exprs() returns raw counts, so the exported "normalized" matrices are raw counts.
#
# Input:  data/li/Li_T.rds (Li et al. T cells; any Li object containing the CD8+ subtypes below works),
#         data/li/metadata_merged.csv
# Output: results/monocle/cd8/cds_tcells.rds
#         data/monocle/cds_tcells_umap.csv                                     (UMAP for Fig. 3A)
#         data/perturbation/CD8/cds_normalized_tcells_sampled_7500.csv, cds_tcells_umap_sampled_7500.csv

library(monocle3)
library(Seurat)
library(data.table)
library(here)

source(here::here("src", "utils", "paths.R"))
out_dir  <- results_dir("monocle", "cd8")
umap_dir <- data_path("monocle");             dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
pert_dir <- data_path("perturbation", "CD8"); dir.create(pert_dir, recursive = TRUE, showWarnings = FALSE)

CD8_TYPES <- c("CD8+T_Naïve/CM", "CD8+T_EM", "CD8+T_EMRA", "CD8+T_EFF.1", "CD8+T_EFF.2", "CD8+T_Act-XCL1/2",
               "CD8+T_RM-CXCL13", "CD8+T_preEX-CXCR4", "CD8+T_preEX-PMCH", "CD8+T_EX-IL10", "CD8+T_EX-CCL4L2")

# ── Cell data set (raw counts + Li metadata) ───────────────────────────────────
rna.li <- readRDS(data_path("li", "Li_T.rds"))
rna.li <- subset(rna.li, subset = annotation %in% CD8_TYPES)
colnames(rna.li) <- paste0("LI_", rna.li$patient, "_", colnames(rna.li))
rna.li <- JoinLayers(rna.li)
expression.matrix <- GetAssayData(rna.li, assay = "RNA", layer = "counts")

sample_metadata <- read.csv(data_path("li", "metadata_merged.csv"))
rownames(sample_metadata) <- sample_metadata$barcode
cell_meta <- sample_metadata[colnames(expression.matrix), ]
gene_meta <- data.frame(gene_short_name = rownames(expression.matrix), row.names = rownames(expression.matrix))

cds <- new_cell_data_set(expression.matrix, cell_metadata = cell_meta, gene_metadata = gene_meta)

# ── Embedding (parameters of the published object cds_1.rds) ──────────────────
cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "patient")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
saveRDS(cds, file.path(out_dir, "cds_tcells.rds"))

# ── Exports ────────────────────────────────────────────────────────────────────
write.csv(as.data.frame(reducedDims(cds)$UMAP), file.path(umap_dir, "cds_tcells_umap.csv"))

set.seed(123)
cds_s <- cds[, sample(colnames(cds), size = 7500)]
write.csv(as.data.frame(reducedDims(cds_s)$UMAP), file.path(pert_dir, "cds_tcells_umap_sampled_7500.csv"))
m <- exprs(cds_s)
fwrite(cbind(gene = rownames(m), as.data.frame(as.matrix(m))),
       file = file.path(pert_dir, "cds_normalized_tcells_sampled_7500.csv"))

cat("Saved ->", out_dir, umap_dir, pert_dir, "\n")
