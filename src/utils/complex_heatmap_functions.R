"""
utils/complex_heatmap_functions.R

Utility functions for building and drawing ComplexHeatmap objects.
Called by figures/figure1/R/01_aucell_heatmap.R and other figure scripts.
"""

library(ComplexHeatmap)
library(circlize)
library(grid)

# ── Matrix preparation ─────────────────────────────────────────────────────────

#' Convert a data frame to a numeric heatmap matrix
#'
#' @param df             data.frame — regulons x cells AUC matrix
#' @param standardize_rows logical — if TRUE, z-score rows (default FALSE)
make_heatmap_matrix_from_df <- function(df, standardize_rows = FALSE) {
  mat <- matrix(
    as.numeric(as.matrix(df)),
    nrow = nrow(df), ncol = ncol(df),
    dimnames = list(rownames(df), colnames(df))
  )
  if (standardize_rows) mat <- t(scale(t(mat)))
  mat
}


#' Load a CSV and convert to numeric heatmap matrix
#'
#' @param counts_filename path to CSV (row names in first column)
#' @param standardize_rows logical — if TRUE, z-score rows (default FALSE)
make_heatmap_matrix <- function(counts_filename, standardize_rows = FALSE) {
  counts <- read.csv(counts_filename, row.names = 1)
  make_heatmap_matrix_from_df(counts, standardize_rows)
}

# ── Core draw function ─────────────────────────────────────────────────────────

#' Draw a ComplexHeatmap to PDF with clustered rows and columns
#'
#' @param heatmap_matrix  numeric matrix — regulons x cells
#' @param filename        output PDF path
#' @param top_annotation  HeatmapAnnotation object for columns (default NULL)
#' @param right_annotation HeatmapAnnotation object for rows (default NULL)
#' @param distance.method clustering distance method (default "pearson")
#' @param clustering.method clustering linkage method (default "average")
#' @param font.size       row label font size (default 6)
draw_complex_heatmap.pdf <- function(heatmap_matrix,
                                     filename,
                                     top_annotation    = NULL,
                                     right_annotation  = NULL,
                                     distance.method   = "pearson",
                                     clustering.method = "average",
                                     font.size         = 6) {
  pdf(filename, width = 25, height = 15)
  
  ht <- Heatmap(
    heatmap_matrix,
    show_row_names            = TRUE,
    show_column_names         = FALSE,
    clustering_distance_rows  = distance.method,
    clustering_distance_columns = distance.method,
    clustering_method_rows    = clustering.method,
    clustering_method_columns = clustering.method,
    top_annotation            = top_annotation,
    right_annotation          = right_annotation,
    na_col                    = "black",
    row_names_gp              = gpar(fontsize = font.size),
    heatmap_legend_param      = list(legend_gp = gpar(fontsize = 8))
  )
  
  ht_drawn <- draw(ht)
  dev.off()
  invisible(ht_drawn)
}


#' Draw a ComplexHeatmap to PDF with clustered rows but pre-ordered columns
#'
#' @inheritParams draw_complex_heatmap.pdf
draw_complex_heatmap.pdf.nonclustered_cells <- function(heatmap_matrix,
                                                        filename,
                                                        top_annotation    = NULL,
                                                        right_annotation  = NULL,
                                                        distance.method   = "pearson",
                                                        clustering.method = "average",
                                                        font.size         = 6) {
  pdf(filename, width = 25, height = 15)
  
  ht <- Heatmap(
    heatmap_matrix,
    show_row_names            = TRUE,
    show_column_names         = FALSE,
    cluster_columns           = FALSE,
    clustering_distance_rows  = distance.method,
    clustering_method_rows    = clustering.method,
    top_annotation            = top_annotation,
    right_annotation          = right_annotation,
    na_col                    = "black",
    row_names_gp              = gpar(fontsize = font.size),
    heatmap_legend_param      = list(legend_gp = gpar(fontsize = 8))
  )
  
  ht_drawn <- draw(ht)
  dev.off()
  invisible(ht_drawn)
}


#' Draw a ComplexHeatmap ordered by a metadata column
#'
#' @param heatmap_matrix  numeric matrix — regulons x cells
#' @param filename        output PDF path
#' @param annotationCol   data.frame — cell metadata (rows = cells)
#' @param metadata_column column name in annotationCol to sort cells by
#' @param colors          named list of colors for annotation
#' @param right_annotation HeatmapAnnotation object for rows (default NULL)
#' @param distance.method row clustering distance method (default "pearson")
#' @param clustering.method row clustering linkage method (default "average")
draw_complex_heatmap.pdf.ordered_by_metadata <- function(heatmap_matrix,
                                                         filename,
                                                         annotationCol,
                                                         metadata_column,
                                                         colors,
                                                         right_annotation  = NULL,
                                                         distance.method   = "pearson",
                                                         clustering.method = "average") {
  stopifnot(all(colnames(heatmap_matrix) %in% rownames(annotationCol)))
  
  col_order        <- order(annotationCol[[metadata_column]])
  ann_ordered      <- annotationCol[col_order, , drop = FALSE]
  mat_ordered      <- heatmap_matrix[, col_order, drop = FALSE]
  top_annotation   <- HeatmapAnnotation(df = ann_ordered, col = colors)
  
  pdf(filename, width = 25, height = 15)
  
  ht <- Heatmap(
    mat_ordered,
    cluster_rows              = TRUE,
    cluster_columns           = FALSE,
    clustering_distance_rows  = distance.method,
    clustering_method_rows    = clustering.method,
    top_annotation            = top_annotation,
    right_annotation          = right_annotation,
    column_title              = paste("Sorted by", metadata_column),
    show_column_names         = FALSE,
    row_names_gp              = gpar(fontsize = 7),
    na_col                    = "black"
  )
  
  draw(ht)
  dev.off()
  message("Saved -> ", filename)
}