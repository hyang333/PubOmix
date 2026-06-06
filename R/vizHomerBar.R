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
#' @param logo_width Numeric. Width of the motif logo column in cm.
#'   Default is 3.
#' @param bar_width Numeric. Width of each bar plot column in cm.
#'   Default is 3.
#' @param row_height Numeric. Height of each row in mm.
#'   Default is 6.
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
vizHomerBar <- function(
  homer_data,
  top_n = NULL,
  pvalue_col = "#2166AC",
  target_col = "#B2182B",
  title = "HOMER Known Motif Enrichment",
  logo_width = 3,
  bar_width = 3,
  row_height = 6,
  ...
) {
  # Convert numeric parameters to grid units
  logo_width <- grid::unit(logo_width, "cm")
  bar_width <- grid::unit(bar_width, "cm")
  row_height <- grid::unit(row_height, "mm")

  # --- Validate input ---
  if (!is(homer_data, "HomerData")) {
    stop(
      "homer_data must be a HomerData object. ",
      "Use readHomerDir() to create one."
    )
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
    gp = gpar(fontsize = 12),
    width = name_width
  )

  # --- Column 3: Motif logos ---
  # Check which logo files exist
  logo_exists <- file.exists(logos) & !is.na(logos)

  if (any(logo_exists)) {
    # Use anno_image for logos that exist; NA for missing ones
    logo_paths <- ifelse(logo_exists, logos, NA_character_)
    # Convert SVG files to temporary PNGs to avoid grImport2 warnings
    logo_paths <- .convert_svg_to_png(logo_paths)
    anno_logo <- anno_image(
      logo_paths,
      which = "row",
      border = FALSE,
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
    Rank = anno_rank,
    TF = anno_tf,
    Motif = anno_logo,
    `-log(Pvalue)` = anno_pval,
    `% Targets` = anno_target,
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


#' Convert SVG logo files to temporary PNG files
#'
#' Avoids grImport2 \code{checkValidSVG} warnings by converting SVG files
#' to raster PNG using the \pkg{rsvg} package. If \pkg{rsvg} is not
#' installed, PNG files in the HOMER directory are used as fallback.
#'
#' @param logo_paths Character vector of logo file paths (may include
#'   PNG, SVG, or NA values).
#' @return Character vector of file paths with SVGs replaced by temp PNGs.
#' @keywords internal
.convert_svg_to_png <- function(logo_paths) {
  svg_idx <- which(
    grepl("\\.svg$", logo_paths, ignore.case = TRUE) &
      !is.na(logo_paths)
  )
  if (length(svg_idx) == 0) {
    return(logo_paths)
  }

  # Try to find matching .logo.png files first (HOMER often generates both)
  for (i in svg_idx) {
    png_alt <- sub("\\.svg$", ".png", logo_paths[i], ignore.case = TRUE)
    if (file.exists(png_alt)) {
      logo_paths[i] <- png_alt
    }
  }

  # Re-check which are still SVGs
  svg_idx <- which(
    grepl("\\.svg$", logo_paths, ignore.case = TRUE) &
      !is.na(logo_paths)
  )
  if (length(svg_idx) == 0) {
    return(logo_paths)
  }

  # Convert remaining SVGs to temp PNGs using rsvg
  if (requireNamespace("rsvg", quietly = TRUE)) {
    for (i in svg_idx) {
      tmp_png <- tempfile(fileext = ".png")
      tryCatch(
        {
          rsvg::rsvg_png(logo_paths[i], file = tmp_png, width = 400)
          logo_paths[i] <- tmp_png
        },
        error = function(e) {
          warning(sprintf("Failed to convert SVG to PNG: %s", logo_paths[i]))
        }
      )
    }
  } else {
    message(
      "Tip: install the 'rsvg' package to avoid SVG rendering warnings: ",
      "install.packages('rsvg')"
    )
  }

  return(logo_paths)
}
