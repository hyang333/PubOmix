#' Visualize HOMER Known Motif Results as a ComplexHeatmap Bar Plot
#'
#' Creates a publication-quality visualization of HOMER known motif discovery
#' results using \pkg{ComplexHeatmap}. The plot contains five columns:
#' \enumerate{
#'   \item Rank of the motifs
#'   \item Transcription factor names
#'   \item Motif logos (PNG images)
#'   \item Horizontal bar plot of -log(P-value)
#'   \item Horizontal bar plot of percent of targets
#' }
#'
#' @param homer_data A \code{\link{HomerData}} object.
#' @param top_n Integer or NULL. Number of top motifs to display. If NULL,
#'   all motifs are shown. Default is 20.
#' @param pvalue_col Character. Color for the -log(P-value) bar plot.
#'   Default is \code{"#2166AC"} (blue).
#' @param target_col Character. Color for the target percentage bar plot.
#'   Default is \code{"#B2182B"} (red).
#' @param title Character. Title for the plot. Default is
#'   \code{"HOMER Known Motif Enrichment"}.
#' @param logo_width \code{grid::unit} object. Width of the motif logo column.
#'   Default is \code{unit(2, "cm")}.
#' @param bar_width \code{grid::unit} object. Width of each bar plot column.
#'   Default is \code{unit(4, "cm")}.
#' @param row_height \code{grid::unit} object. Height of each row.
#'   Default is \code{unit(8, "mm")}.
#' @param ... Additional arguments passed to \code{ComplexHeatmap::draw()}.
#'
#' @return A \code{ComplexHeatmap::Heatmap} object (invisibly). The plot is
#'   drawn as a side effect.
#'
#' @examples
#' \dontrun{
#' homer_data <- readHomerDir("/path/to/homer/output")
#' vizHomerBar(homer_data, top_n = 15)
#' }
#'
#' @import ComplexHeatmap
#' @import grid
#' @importFrom png readPNG
#' @importFrom circlize colorRamp2
#' @export
vizHomerBar <- function(homer_data,
                        top_n = 20,
                        pvalue_col = "#2166AC",
                        target_col = "#B2182B",
                        title = "HOMER Known Motif Enrichment",
                        logo_width = grid::unit(2, "cm"),
                        bar_width = grid::unit(4, "cm"),
                        row_height = grid::unit(8, "mm"),
                        ...) {

  # --- Validate input ---
  if (!is(homer_data, "HomerData")) {
    stop("homer_data must be a HomerData object. ",
         "Use readHomerDir() to create one.")
  }

  md <- homer_data@motif_data
  logos <- homer_data@logo_files

  # Subset to top_n
  if (!is.null(top_n) && top_n < nrow(md)) {
    md <- md[seq_len(top_n), ]
    logos <- logos[seq_len(top_n)]
  }

  n_motifs <- nrow(md)

  # --- Prepare data ---
  ranks <- md$rank
  tf_names <- md$motif_name
  neg_log_pval <- md$neg_log_p_value
  target_pct <- md$target_pct

  # --- Column 1: Rank annotation (text) ---
  anno_rank <- anno_text(
    x = as.character(ranks),
    which = "row",
    gp = gpar(fontsize = 9, fontface = "bold"),
    width = unit(1, "cm")
  )

  # --- Column 2: TF name annotation (text) ---
  # Calculate width based on max string length
  max_name_len <- max(nchar(tf_names))
  name_width <- unit(min(max_name_len * 2.2, 80), "mm")

  anno_tf <- anno_text(
    x = tf_names,
    which = "row",
    gp = gpar(fontsize = 8),
    width = name_width
  )

  # --- Column 3: Motif logos ---
  # Check which logo files exist
  logo_exists <- file.exists(logos) & !is.na(logos)

  if (any(logo_exists)) {
    # Use anno_image for logos that exist; NA for missing ones
    logo_paths <- ifelse(logo_exists, logos, NA_character_)
    anno_logo <- anno_image(
      logo_paths,
      which = "row",
      width = logo_width,
      space = unit(1, "mm")
    )
  } else {
    warning("No motif logo files found. Logo column will be empty.")
    anno_logo <- anno_empty(which = "row", width = logo_width)
  }

  # --- Column 4: -log(P-value) bar plot ---
  anno_pval <- anno_barplot(
    x = neg_log_pval,
    which = "row",
    bar_width = 0.8,
    width = bar_width,
    gp = gpar(fill = pvalue_col, col = NA),
    axis_param = list(
      side = "top",
      gp = gpar(fontsize = 7)
    )
  )

  # --- Column 5: Percent of targets bar plot ---
  anno_target <- anno_barplot(
    x = target_pct,
    which = "row",
    bar_width = 0.8,
    width = bar_width,
    gp = gpar(fill = target_col, col = NA),
    axis_param = list(
      side = "top",
      gp = gpar(fontsize = 7)
    )
  )

  # --- Assemble row annotations ---
  ha <- rowAnnotation(
    Rank       = anno_rank,
    TF         = anno_tf,
    Motif      = anno_logo,
    `-log(Pvalue)` = anno_pval,
    `% Targets`    = anno_target,
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    gap = unit(2, "mm")
  )

  # --- Create a minimal dummy heatmap matrix ---
  # ComplexHeatmap requires a matrix; we use a 1-column invisible matrix
  mat <- matrix(NA, nrow = n_motifs, ncol = 1)
  rownames(mat) <- tf_names

  # --- Draw the heatmap ---
  ht <- Heatmap(
    mat,
    name = "dummy",
    show_heatmap_legend = FALSE,
    width = unit(0, "mm"),
    right_annotation = ha,
    row_names_side = NULL,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title = NULL,
    column_title = title,
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    height = row_height * n_motifs
  )

  # Draw the heatmap
  draw(ht, ...)

  invisible(ht)
}
