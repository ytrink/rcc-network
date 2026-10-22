# src/figure1/R/02_aucell_heatmap.R
#
# Figure 1C: AUCell regulon activity heatmap (Li et al., ccRCC patients), cells and regulons hierarchically
# clustered (Pearson correlation distance, average linkage); rows coloured by Leiden community.
# Port of stuff_to_run_on_the_side/plot_auc_heatmaps_all/plot_auc_heatmap_li_ALL_POS_ccRCC.R
#
#   - cells: up to 5,000 cells per broad cell type (9 types: B_cell, CD4_T, CD8_T, NK, Monocyte, Macrophage,
#     Endothelial, Fibroblast, Tumor). The published sample is data/li/Li_sampled_broad2_5000_barcodes.csv;
#     if that file is absent a new random sample is drawn (set.seed(123)), which will differ slightly.
#   - regulons with >= 5 target genes detected in the dataset
#
# Input:  results/figure1/aucell/auc_matrix_stand.csv   (01_aucell_scores.R)
#         data/li/Li_integrated.rds, data/li/metadata_merged.csv
#         data/scplus/leiden_RB_03_pos.csv               (Table S1)
# Output: results/figure1/heatmap_pearson.pdf (heatmap body rasterized)

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(data.table)
library(here)

source(here::here("src", "utils", "paths.R"))
source(here::here("src", "utils", "complex_heatmap_functions.R"))
out_dir <- results_dir("figure1")

N_CELLS_PER_TYPE <- 5000
MIN_REGULON_SIZE <- 5

# ── AUCell matrix ──────────────────────────────────────────────────────────────
heatmap_matrix <- as.data.frame(fread(file.path(RESULTS_DIR, "figure1", "aucell", "auc_matrix_stand.csv")))
rownames(heatmap_matrix) <- heatmap_matrix[[1]]
heatmap_matrix[[1]] <- NULL

regulon_size   <- as.numeric(sub(".*\\((\\d+)g\\)", "\\1", rownames(heatmap_matrix)))
heatmap_matrix <- heatmap_matrix[regulon_size >= MIN_REGULON_SIZE, ]

# ── Cells: 5,000 per broad cell type ───────────────────────────────────────────
rna.li <- readRDS(data_path("li", "Li_integrated.rds"))
colnames(rna.li) <- paste0("LI_", rna.li$patient, "_", colnames(rna.li))
rna.li <- rna.li[, intersect(colnames(rna.li), colnames(heatmap_matrix))]   # ccRCC cells scored by AUCell

barcodes_file <- data_path("li", "Li_sampled_broad2_5000_barcodes.csv")
if (file.exists(barcodes_file)) {
  cells_keep <- intersect(read.csv(barcodes_file)[[1]], colnames(rna.li))
} else {
  set.seed(123)
  cells_keep <- rna.li@meta.data %>%
    tibble::rownames_to_column("barcode") %>%
    group_by(annotation_broad) %>%
    slice_sample(n = N_CELLS_PER_TYPE) %>%
    pull(barcode)
}
rna.li         <- rna.li[, cells_keep]
heatmap_matrix <- heatmap_matrix[, colnames(rna.li)]

# ── Column annotations (as in Fig. 1C) ─────────────────────────────────────────
sample_metadata <- read.csv(data_path("li", "metadata_merged.csv"))
rownames(sample_metadata) <- sample_metadata$barcode
ann <- sample_metadata[colnames(heatmap_matrix), ]
ann[] <- lapply(ann, function(x) { if (is.character(x)) x[x == ""] <- "Unknown"; x })

ann_df <- data.frame(
  patient            = ann$patient,
  summaryDescription = ann$summaryDescription,   # tissue of origin
  annotation_broad_1 = case_when(
    rna.li$annotation_broad %in% c("CD4_T", "CD8_T", "NK")   ~ "T_NK",
    rna.li$annotation_broad %in% c("Macrophage", "Monocyte") ~ "Macrophage_Monocyte",
    TRUE ~ rna.li$annotation_broad
  ),
  stage              = as.character(ann$stage),
  grade              = as.character(ann$grade),
  region_ii          = recode(ann$region,
    "a" = "tumor core", "c" = "tumor core", "d" = "tumor core", "e" = "tumor core",
    "g" = "tumor normal interface", "n" = "normal kidney", "b" = "peripheral blood",
    "i" = "adrenal metastasis", "h" = "normal adrenal gland", "t" = "thrombus", "f" = "perinephric fat"),
  row.names = colnames(heatmap_matrix)
)
ann_df[] <- lapply(ann_df, as.factor)

ann_colors <- list(
  patient = c(PD43824 = "#1f77b4", PD43948 = "#ff7f0e", PD44966 = "#2ca02c", PD44967 = "#d62728",
              PD45814 = "#9467bd", PD45815 = "#8c564b", PD45816 = "#e377c2", PD47171 = "#7f7f7f",
              PD47465 = "#bcbd22", PD47512 = "#17becf"),
  summaryDescription = c(Blood = "red", Metastasis = "gold", Tumour = "black", "Normal kidney" = "blue",
                         "Tumour-normal" = "green", "Normal adrenal" = "purple", Fat = "orange", Thrombus = "cyan"),
  annotation_broad_1 = c(B_cell = "darkorange", T_NK = "dodgerblue", Endothelial = "forestgreen",
                         Fibroblast = "purple", Macrophage_Monocyte = "deeppink", Tumor = "black"),
  stage = c("3a" = "red", "1a" = "green", "1b" = "blue", Unknown = "black"),
  grade = c("4" = "red", "2" = "green", "3" = "blue"),
  region_ii = c("tumor core" = "black", "tumor normal interface" = "green", "normal kidney" = "blue",
                "peripheral blood" = "red", "adrenal metastasis" = "gold", "normal adrenal gland" = "purple",
                "perinephric fat" = "orange", thrombus = "cyan")
)
top_annotation <- HeatmapAnnotation(df = ann_df, col = ann_colors)

# ── Row annotation: Leiden community (Table S1) ────────────────────────────────
leiden_clusters <- read.csv(data_path("scplus", "leiden_RB_03_pos.csv"))
gene_names      <- sub(" \\(.*\\)", "", rownames(heatmap_matrix))
community <- sapply(gene_names, function(g) {
  hit <- names(leiden_clusters)[sapply(leiden_clusters, function(col) g %in% col)]
  if (length(hit)) hit[1] else NA
})
community_colors <- c(Community_0 = "#1f77b4", Community_1 = "#ff7f0e", Community_2 = "#9467bd",
                      Community_3 = "#2ca02c", Community_4 = "#d62728", Community_5 = "#8c564b",
                      Community_6 = "#17becf")
row_ha <- rowAnnotation(leiden_cluster = factor(community, levels = names(community_colors)),
                        col = list(leiden_cluster = community_colors), show_annotation_name = FALSE)

# ── Plot (clustered cells and regulons) ───────────────────────────────────────
draw_complex_heatmap.pdf(
  heatmap_matrix    = make_heatmap_matrix_from_df(heatmap_matrix, standardize_rows = FALSE),
  filename          = file.path(out_dir, "heatmap_pearson.pdf"),
  top_annotation    = top_annotation,
  right_annotation  = row_ha,
  distance.method   = "pearson",
  clustering.method = "average"
)
cat("Saved ->", file.path(out_dir, "heatmap_pearson.pdf"), "\n")
