# preprocessing/scrna/yu/04_infercnv_yu.R
#
# Large-scale chromosomal copy number variation (CNV) estimation
# using inferCNV (v1.20.0) to confirm malignancy of putative tumor cells,
# using NK cells as the reference population (Yu et al. cohort).
#
# Input:  data/yu/scrna/rna_annotated.rds
#         data/annotations/hg38_gencode_v27.txt
# Output: preprocessing/scrna/yu/infercnv_results/infercnv_tumor_nk.rds

library(infercnv)
library(Seurat)
library(dplyr)
library(here)

set.seed(101)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir <- here::here()
data_dir <- file.path(root_dir, "data", "yu", "scrna")
ref_dir  <- file.path(root_dir, "data", "annotations")
out_dir  <- file.path(root_dir, "preprocessing", "scrna", "yu", "infercnv_results")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load Annotated Seurat Object ───────────────────────────────────────────────
cat("Loading annotated Seurat object...\n")
rna_path <- file.path(data_dir, "rna_annotated.rds")
if (!file.exists(rna_path)) {
  rna_path <- file.path(data_dir, "rna_integrated.rds")
}
rcc <- readRDS(rna_path)
DefaultAssay(rcc) <- "RNA"

# Extract raw counts matrix
counts_raw <- GetAssayData(rcc, assay = "RNA", slot = "counts")

# ── Select NK Reference and Tumor Cells ────────────────────────────────────────
cat("Selecting NK reference and tumor cells...\n")
cell_labels <- data.frame(
  celltype = rcc$cell_type.2,
  barcode  = colnames(rcc),
  batch    = rcc$batch,
  stringsAsFactors = FALSE
)

idx_ref   <- which(cell_labels$celltype == "NK")
idx_tumor <- which(cell_labels$celltype == "tumor")
idx_take  <- c(idx_ref, idx_tumor)

cell_subset <- cell_labels[idx_take, ]

# Prefix tumor cells with sample batch ID to assess inter-patient heterogeneity
cell_subset <- cell_subset %>%
  mutate(celltype = ifelse(celltype == "tumor", paste0("tumor_", batch), celltype)) %>%
  select(barcode, celltype)

# Write annotation file for inferCNV
ann_file <- file.path(out_dir, "annotations_tumor_nk.txt")
write.table(cell_subset, ann_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ── Gene Ordering File ─────────────────────────────────────────────────────────
gene_order_file <- file.path(ref_dir, "hg38_gencode_v27.txt")
if (!file.exists(gene_order_file)) {
  gene_order_file <- file.path(data_dir, "hg38_gencode_v27.txt")
}

counts_subset <- counts_raw[, cell_subset$barcode]

# ── Create & Run inferCNV ──────────────────────────────────────────────────────
cat("Creating inferCNV object...\n")
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix   = counts_subset,
  annotations_file    = ann_file,
  delim               = "\t",
  gene_order_file     = gene_order_file,
  ref_group_names     = c("NK"),
  max_cells_per_group = 3000
)

cat("Running inferCNV analysis (denoise = TRUE, HMM = FALSE)...\n")
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff            = 0.1,
  out_dir           = out_dir,
  cluster_by_groups = TRUE,
  denoise           = TRUE,
  HMM               = FALSE
)

# ── Save Result ────────────────────────────────────────────────────────────────
out_rds <- file.path(out_dir, "infercnv_tumor_nk.rds")
saveRDS(infercnv_obj, out_rds)
cat("Saved inferCNV object ->", out_rds, "\n")
