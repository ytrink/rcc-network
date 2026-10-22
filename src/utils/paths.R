# src/utils/paths.R
#
# Shared path helpers. The repository root is found via the `.here` file, so scripts can be run
# from any working directory. All inputs live under data/ (see README for the expected layout and
# the Zenodo download); all outputs are written under results/.

ROOT_DIR    <- here::here()
DATA_DIR    <- file.path(ROOT_DIR, "data")
RESULTS_DIR <- file.path(ROOT_DIR, "results")

# Path to an input file, e.g. data_path("scplus", "scplus_direct_extended.graphml")
data_path <- function(...) file.path(DATA_DIR, ...)

# Output directory (created if needed), e.g. results_dir("figure2")
results_dir <- function(...) {
  d <- file.path(RESULTS_DIR, ...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
