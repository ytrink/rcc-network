# network_analysis/01_aucell_scores.R
#
# Computes AUCell regulon activity scores for Li et al. scRNA-seq data
# using SCENIC+ eRegulon gene sets.
#
# Input:  data/li/Li_integrated.rds
#         data/scplus_direct_extended_POS.xlsx
# Output: network_analysis/outputs/auc_matrix.csv
#         network_analysis/outputs/auc_matrix_standardized.csv
#         network_analysis/outputs/cell_annotations.csv

library(Seurat)
library(AUCell)
library(GSEABase)
library(readxl)
library(data.table)
library(here)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir    <- here::here()
data_dir    <- file.path(root_dir, "data")
out_dir     <- file.path(root_dir, "network_analysis", "outputs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Helper: build gene sets from Excel ────────────────────────────────────────
build_gene_sets <- function(filepath) {
  sheet_names <- excel_sheets(filepath)
  gene_sets   <- lapply(sheet_names, function(sheet) {
    targets <- read_excel(filepath, sheet = sheet)$Targets
    GeneSet(targets, setName = sheet)
  })
  names(gene_sets) <- sheet_names
  GeneSetCollection(gene_sets)
}

# ── Helper: run AUCell ────────────────────────────────────────────────────────
run_aucell <- function(gene_sets_filepath, expr_matrix, out_dir, n_cores = 24) {
  
  # Build gene sets
  gene_set_collection <- build_gene_sets(gene_sets_filepath)
  gene_sets <- subsetGeneSets(gene_set_collection, rownames(expr_matrix))
  gene_sets <- setGeneSetNames(
    gene_sets,
    newNames = paste0(names(gene_sets), " (", nGenes(gene_sets), "g)")
  )
  
  # Rankings
  message("Building cell rankings...")
  cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE, nCores = n_cores)
  
  # AUC scores
  message("Running AUC...")
  cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings)
  
  # Save AUC matrix
  auc_matrix <- as.data.frame(cells_AUC@assays@data@listData$AUC)
  fwrite(
    data.table::data.table(Regulon = rownames(auc_matrix), auc_matrix),
    file = file.path(out_dir, "auc_matrix.csv")
  )
  
  message("Done.")
  return(cells_AUC)
}

# ── Load data ──────────────────────────────────────────────────────────────────
data.li <- readRDS(file.path(data_dir, "li", "Li_integrated.rds"))
colnames(data.li) <- paste0("LI_", data.li$patient, "_", colnames(data.li))

data.li <- NormalizeData(data.li, assay = "RNA", normalization.method = "LogNormalize",
                         scale.factor = 10000, verbose = FALSE)

expr_matrix <- GetAssayData(JoinLayers(data.li), assay = "RNA", slot = "data")

# ── Run AUCell ─────────────────────────────────────────────────────────────────
run_aucell(
  gene_sets_filepath = file.path(data_dir, "scplus_direct_extended_POS.xlsx"),
  expr_matrix        = expr_matrix,
  out_dir            = out_dir
)

# ── Save cell annotations ──────────────────────────────────────────────────────
write.csv(data.li$annotation, file.path(out_dir, "cell_annotations.csv"))

# ── Standardize AUC matrix ────────────────────────────────────────────────────
auc_matrix <- fread(file.path(out_dir, "auc_matrix.csv")) |>
  as.data.frame()
rownames(auc_matrix) <- auc_matrix[[1]]
auc_matrix[[1]]      <- NULL
auc_matrix_std       <- t(scale(t(auc_matrix)))

write.csv(auc_matrix_std, file.path(out_dir, "auc_matrix_standardized.csv"))
cat("Saved standardized AUC matrix\n")