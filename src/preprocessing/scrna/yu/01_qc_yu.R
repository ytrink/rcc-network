# preprocessing/scrna/01_qc_yu.R
#
# Quality control, per-sample normalization, doublet removal,
# and merging of 19 ccRCC scRNA-seq samples (Yu et al., GSE207493).
#
# Input:  Raw 10X feature-barcode matrices per sample (data/yu/scrna/<SAMPLE>/)
# Output: data/yu/scrna/rna_merged.rds

library(Seurat)
library(scDblFinder)
library(ggpubr)
library(dplyr)

set.seed(101)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir  <- rprojroot::find_rstudio_root_file()
data_dir  <- file.path(root_dir, "data", "yu", "scrna")
out_dir   <- file.path(root_dir, "data", "yu", "scrna", "per_sample")
qc_dir    <- file.path(root_dir, "preprocessing", "scrna", "qc_plots")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir,  recursive = TRUE, showWarnings = FALSE)

# ── Parameters ─────────────────────────────────────────────────────────────────
SAMPLES <- c(
  "RCC81",  "RCC84",  "RCC86",  "RCC87",  "RCC94",  "RCC96",  "RCC99",
  "RCC100", "RCC101", "RCC103", "RCC104", "RCC106", "RCC112", "RCC113",
  "RCC114", "RCC115", "RCC116", "RCC119", "RCC120"
)

QC_FILTERS <- list(
  min_features  = 200,
  max_features  = 6000,
  min_counts    = 1000,
  max_mt        = 10
)

PCA_DIMS      <- 1:30
RESOLUTIONS   <- seq(0.5, 2, by = 0.1)

# ── Helper: process one sample ─────────────────────────────────────────────────
process_sample <- function(sample_id) {
  
  # Load
  counts <- Read10X(data.dir = file.path(data_dir, sample_id))
  obj    <- CreateSeuratObject(
    counts       = counts,
    project      = sample_id,
    min.cells    = 3,
    min.features = QC_FILTERS$min_features
  )
  
  # Mitochondrial fraction
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # QC plots
  p1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p3 <- ggdensity(obj@meta.data, x = "nCount_RNA",   title = sample_id)
  p4 <- ggdensity(obj@meta.data, x = "nFeature_RNA", title = sample_id)
  p5 <- ggdensity(obj@meta.data, x = "percent.mt",   title = sample_id)
  
  # Filter
  obj <- subset(
    obj,
    subset = nFeature_RNA > QC_FILTERS$min_features &
      nFeature_RNA < QC_FILTERS$max_features &
      nCount_RNA   > QC_FILTERS$min_counts   &
      percent.mt   < QC_FILTERS$max_mt
  )
  
  # Normalize + cluster (for doublet detection)
  obj <- SCTransform(obj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = PCA_DIMS, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  obj <- RunUMAP(obj, dims = PCA_DIMS, verbose = FALSE)
  
  # Elbow plot
  p6 <- ElbowPlot(obj, ndims = 50)
  
  # Resolution sweep plots
  res_plots <- lapply(RESOLUTIONS, function(res) {
    DimPlot(obj, reduction = "umap", label = TRUE,
            group.by = paste0("SCT_snn_res.", res))
  })
  
  # Doublet detection
  set.seed(123)
  sce          <- scDblFinder(GetAssayData(obj), clusters = TRUE)
  obj$doublet  <- sce$scDblFinder.class
  obj          <- subset(obj, subset = doublet == "singlet")
  obj$batch    <- sample_id
  
  # Clean metadata
  obj@meta.data <- obj@meta.data[, c(
    "orig.ident", "batch", "nCount_RNA", "nFeature_RNA", "percent.mt", "doublet"
  )]
  
  list(
    obj   = obj,
    plots = list(p1, p2, p3, p4, p5, p6, res_plots)
  )
}

# ── Process all samples ────────────────────────────────────────────────────────
results <- lapply(SAMPLES, process_sample)
names(results) <- SAMPLES

# ── Save QC plots ──────────────────────────────────────────────────────────────
pdf(file.path(qc_dir, "qc_plots_yu.pdf"))
for (s in SAMPLES) {
  plots <- results[[s]]$plots
  print(plots[[1]])  # VlnPlot
  print(plots[[2]])  # FeatureScatter
  print(plots[[3]])  # density nCount
  print(plots[[4]])  # density nFeature
  print(plots[[5]])  # density percent.mt
  print(plots[[6]])  # ElbowPlot
  for (p in plots[[7]]) print(p)  # resolution sweep
}
dev.off()

# ── Merge all samples ──────────────────────────────────────────────────────────
objects    <- lapply(SAMPLES, function(s) results[[s]]$obj)
rcc.merged <- merge(
  objects[[1]],
  y           = objects[-1],
  add.cell.ids = SAMPLES
)

DefaultAssay(rcc.merged) <- "RNA"

# ── Save ───────────────────────────────────────────────────────────────────────
saveRDS(rcc.merged, file.path(data_dir, "rna_merged.rds"))
cat("Saved merged object ->", file.path(data_dir, "rna_merged.rds"), "\n")