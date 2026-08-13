# preprocessing/scatac/yu/01_qc_atac_yu.R
#
# Creates consensus peak set, builds chromatin accessibility matrices,
# performs QC and filtering for Yu et al. scATAC-seq data.
#
# Input:  data/yu/scatac/<SAMPLE>/peaks.bed
#         data/yu/scatac/<SAMPLE>/singlecell.csv
#         data/yu/scatac/<SAMPLE>/fragments.tsv.gz
#         data/annotations_hsapiensv86.rds
# Output: data/yu/scatac/atac_merged_filtered.rds

library(Seurat)
library(Signac)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(here)

set.seed(101)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir  <- here::here()
atac_dir  <- file.path(root_dir, "data", "yu", "scatac")
out_dir   <- file.path(root_dir, "data", "yu", "scatac")
qc_dir    <- file.path(root_dir, "preprocessing", "scatac", "yu", "qc_plots")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# ── Parameters ─────────────────────────────────────────────────────────────────
SAMPLES <- c(
  "RCC81",  "RCC84",  "RCC86",  "RCC87",  "RCC94",  "RCC96",  "RCC99",
  "RCC100", "RCC101", "RCC103", "RCC104", "RCC106", "RCC112", "RCC113",
  "RCC114", "RCC115", "RCC116", "RCC119", "RCC120"
)

QC_FILTERS <- list(
  min_fragments     = 1000,
  max_fragments     = 20000,
  min_pct_in_peaks  = 15,
  max_blacklist     = 0.05,
  max_nucleosome    = 4,
  min_tss           = 3
)

MIN_PASSED_FILTERS <- 500

# ── Helper: load one sample ────────────────────────────────────────────────────
load_sample <- function(sample_id, combined_peaks) {
  
  sample_dir <- file.path(atac_dir, sample_id)
  
  # Load and filter metadata
  md <- read.table(
    file.path(sample_dir, "singlecell.csv"),
    stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1
  )[-1, ]
  md <- md[md$passed_filters > MIN_PASSED_FILTERS, ]
  
  # Fragment object
  frags <- CreateFragmentObject(
    path  = file.path(sample_dir, "fragments.tsv.gz"),
    cells = rownames(md)
  )
  
  # Count matrix
  counts <- FeatureMatrix(
    fragments = frags,
    features  = combined_peaks,
    cells     = rownames(md)
  )
  
  # Seurat object
  obj        <- CreateSeuratObject(
    CreateChromatinAssay(counts, fragments = frags),
    assay     = "ATAC",
    meta.data = md
  )
  obj$batch  <- sample_id
  obj
}

# ── Build consensus peak set ───────────────────────────────────────────────────
cat("Building consensus peak set...\n")

peak_list <- lapply(SAMPLES, function(s) {
  peaks <- read.table(
    file.path(atac_dir, s, "peaks.bed"),
    col.names = c("chr", "start", "end")
  )
  makeGRangesFromDataFrame(peaks)
})

combined_peaks <- reduce(do.call(c, peak_list))
combined_peaks <- combined_peaks[width(combined_peaks) > 20 &
                                   width(combined_peaks) < 10000]

cat("Consensus peaks:", length(combined_peaks), "\n")

# ── Load all samples ───────────────────────────────────────────────────────────
cat("Loading samples...\n")
objects <- lapply(SAMPLES, load_sample, combined_peaks = combined_peaks)
names(objects) <- SAMPLES

# ── Merge ──────────────────────────────────────────────────────────────────────
combined <- merge(
  objects[[1]],
  y            = objects[-1],
  add.cell.ids = SAMPLES
)

# ── Add genome annotations ─────────────────────────────────────────────────────
Annotation(combined) <- readRDS(file.path(root_dir, "data", "annotations_hsapiensv86.rds"))

# ── QC metrics ────────────────────────────────────────────────────────────────
cat("Computing QC metrics...\n")
combined <- NucleosomeSignal(combined)
combined <- TSSEnrichment(combined, fast = FALSE, assay = "ATAC", verbose = FALSE)

combined$pct_reads_in_peaks <- combined$peak_region_fragments /
  combined$passed_filters * 100
combined$blacklist_ratio    <- combined$blacklist_region_fragments /
  combined$peak_region_fragments
combined$high.tss           <- ifelse(combined$TSS.enrichment > 2, "High", "Low")
combined$nucleosome_group   <- ifelse(combined$nucleosome_signal > 4, "NS > 4", "NS < 4")

# ── QC plots ───────────────────────────────────────────────────────────────────
pdf(file.path(qc_dir, "atac_qc_plots.pdf"))
print(TSSPlot(combined, group.by = "high.tss") + NoLegend())
print(FragmentHistogram(combined, group.by = "nucleosome_group"))
print(VlnPlot(
  combined,
  features = c("pct_reads_in_peaks", "peak_region_fragments",
               "TSS.enrichment", "blacklist_ratio", "nucleosome_signal"),
  pt.size = 0, ncol = 5
))
dev.off()

# ── Filter cells ───────────────────────────────────────────────────────────────
combined.pro <- subset(
  combined,
  subset = peak_region_fragments > QC_FILTERS$min_fragments  &
    peak_region_fragments < QC_FILTERS$max_fragments  &
    pct_reads_in_peaks    > QC_FILTERS$min_pct_in_peaks &
    blacklist_ratio       < QC_FILTERS$max_blacklist  &
    nucleosome_signal     < QC_FILTERS$max_nucleosome &
    TSS.enrichment        > QC_FILTERS$min_tss
)

combined.pro@meta.data$orig.ident <- gsub("_.*", "", rownames(combined.pro@meta.data))
cat("Cells after filtering:", ncol(combined.pro), "\n")

# ── Save ───────────────────────────────────────────────────────────────────────
saveRDS(combined.pro, file.path(out_dir, "atac_merged_filtered.rds"))
cat("Saved ->", file.path(out_dir, "atac_merged_filtered.rds"), "\n")