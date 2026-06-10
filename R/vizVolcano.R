#' Create a Volcano Plot from DESeq2 Results
#'
#' Generates a publication-quality volcano plot from DESeq2 differential
#' expression results. Genes are colored by significance and fold-change
#' magnitude:
#' \itemize{
#'   \item \strong{Red}: significant (adjusted p < \code{padj_cutoff}) and
#'     \code{|log2FoldChange| > lfc_threshold} (large change)
#'   \item \strong{Light blue}: significant but
#'     \code{|log2FoldChange| <= lfc_threshold} (small change)
#'   \item \strong{Grey}: not significant
#' }
#'
#' Gene counts for four significant zones are annotated directly on the
#' plot near the top of each zone.
#'
#' @param res A \code{data.frame} (or tibble) of DESeq2 results. Must contain
#'   columns \code{log2FoldChange} and \code{padj}. Typically produced by
#'   \code{DESeq2::results()} or read from a saved CSV.
#' @param lfc_threshold Numeric. The log2 fold-change threshold used to
#'   separate "large" from "small" significant changes. Default is \code{1}.
#' @param padj_cutoff Numeric. Adjusted p-value cutoff for significance.
#'   Default is \code{0.05}.
#' @param col_up Character. Color for significant genes with
#'   \code{|log2FoldChange| > lfc_threshold}. Default is \code{"#E41A1C"}
#'   (red).
#' @param col_low Character. Color for significant genes with
#'   \code{|log2FoldChange| <= lfc_threshold}. Default is \code{"#6BAED6"}
#'   (light blue).
#' @param col_ns Character. Color for non-significant genes.
#'   Default is \code{"grey70"}.
#' @param point_size Numeric. Size of the points. Default is \code{1.2}.
#' @param point_alpha Numeric. Transparency of the points (0--1).
#'   Default is \code{0.7}.
#' @param title Character. Plot title. Default is \code{"Volcano Plot"}.
#' @param xlab Character or expression. X-axis label. Default is
#'   \code{"log2(Fold Change)"}.
#' @param ylab Expression or character. Y-axis label. Default uses
#'   \code{expression(-log[10]~italic(FDR))}.
#' @param xlim Numeric vector of length 2 giving the x-axis limits
#'   (e.g. \code{c(-5, 5)}), or \code{NULL} for automatic limits.
#'   Default is \code{NULL}.
#' @param ylim Numeric vector of length 2 giving the y-axis limits
#'   (e.g. \code{c(0, 50)}), or \code{NULL} for automatic limits.
#'   Default is \code{NULL}.
#' @param label_genes Character vector of gene symbols to label on the plot,
#'   or \code{NULL} (no labels). Default is \code{NULL}.
#' @param gene_col Character. Column name in \code{res} that holds gene
#'   symbols, used when \code{label_genes} is not NULL.
#'   Default is \code{"symbol"}.
#'
#' @return A \code{ggplot} object (invisibly). The plot is drawn as a side
#'   effect.
#'
#' @examples
#' \dontrun{
#' res <- read.csv("EUC021_LFC_AL.csv")
#' vizVolcano(res, lfc_threshold = 1)
#' vizVolcano(res, lfc_threshold = 0.5, padj_cutoff = 0.01)
#' vizVolcano(res, lfc_threshold = 1, xlim = c(-6, 6), ylim = c(0, 40))
#' }
#'
#' @importFrom stats complete.cases
#' @export
vizVolcano <- function(
  res,
  lfc_threshold = 1,
  padj_cutoff = 0.05,
  col_up = "#E41A1C",
  col_low = "#6BAED6",
  col_ns = "grey70",
  point_size = 1.2,
  point_alpha = 0.7,
  title = "Volcano Plot",
  xlab = "log2(Fold Change)",
  ylab = expression(-log[10]~italic(FDR)),
  xlim = NULL,
  ylim = NULL,
  label_genes = NULL,
  gene_col = "symbol"
) {

  # ── Check that ggplot2 is available ────────────────────────────────────────
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for vizVolcano(). ",
      "Install it with: install.packages('ggplot2')"
    )
  }

  # ── Validate input ────────────────────────────────────────────────────────
  if (!is.data.frame(res)) {
    stop("'res' must be a data.frame of DESeq2 results.")
  }

  required_cols <- c("log2FoldChange", "padj")
  missing_cols <- setdiff(required_cols, colnames(res))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s) in 'res': ",
      paste(missing_cols, collapse = ", "),
      ". Expected DESeq2 results with 'log2FoldChange' and 'padj'."
    )
  }

  # ── Prepare data ──────────────────────────────────────────────────────────
  # Remove rows with NA in key columns
  df <- res[complete.cases(res[, c("log2FoldChange", "padj")]), ]

  # Compute -log10(padj)
  df$neg_log10_padj <- -log10(df$padj)

  # ── Classify genes into 3 color categories ────────────────────────────────
  df$category <- ifelse(
    df$padj >= padj_cutoff,
    "NS",
    ifelse(
      abs(df$log2FoldChange) > lfc_threshold,
      "Sig_High",
      "Sig_Low"
    )
  )

  df$category <- factor(df$category, levels = c("Sig_High", "Sig_Low", "NS"))

  # Named colour vector for scale_color_manual
  color_map <- c(
    "Sig_High" = col_up,
    "Sig_Low"  = col_low,
    "NS"       = col_ns
  )

  # ── Count significant genes in 4 zones ────────────────────────────────────
  sig <- df[df$padj < padj_cutoff, ]
  n_down_high <- sum(sig$log2FoldChange < -lfc_threshold)
  n_down_low  <- sum(sig$log2FoldChange >= -lfc_threshold & sig$log2FoldChange < 0)
  n_up_low    <- sum(sig$log2FoldChange >= 0 & sig$log2FoldChange <= lfc_threshold)
  n_up_high   <- sum(sig$log2FoldChange > lfc_threshold)
  n_ns        <- sum(df$category == "NS")

  # ── Determine axis limits for annotation placement ────────────────────────
  # Use user-supplied limits if given, otherwise compute from data
  if (!is.null(xlim)) {
    x_lo <- xlim[1]
    x_hi <- xlim[2]
  } else {
    x_lo <- min(df$log2FoldChange, na.rm = TRUE)
    x_hi <- max(df$log2FoldChange, na.rm = TRUE)
  }
  if (!is.null(ylim)) {
    y_hi <- ylim[2]
  } else {
    y_hi <- max(df$neg_log10_padj, na.rm = TRUE)
  }

  # Zone midpoints (x) for annotation labels
  x_mid_down_high <- (x_lo + (-lfc_threshold)) / 2
  x_mid_down_low  <- ((-lfc_threshold) + 0) / 2
  x_mid_up_low    <- (0 + lfc_threshold) / 2
  x_mid_up_high   <- (lfc_threshold + x_hi) / 2

  # Y position for count labels: near top of plot
  y_label <- y_hi * 0.95

  # Build annotation data.frame
  count_labels <- data.frame(
    x = c(x_mid_down_high, x_mid_down_low, x_mid_up_low, x_mid_up_high),
    y = rep(y_label, 4),
    label = as.character(c(n_down_high, n_down_low, n_up_low, n_up_high)),
    col = c(col_up, col_low, col_low, col_up),
    stringsAsFactors = FALSE
  )

  # ── Build the plot ────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$log2FoldChange,
      y = .data$neg_log10_padj,
      color = .data$category
    )
  ) +
    ggplot2::geom_point(
      size = point_size,
      alpha = point_alpha
    ) +
    ggplot2::scale_color_manual(
      values = color_map,
      labels = c(
        "Sig_High" = paste0("|LFC| > ", lfc_threshold),
        "Sig_Low"  = paste0("|LFC| \u2264 ", lfc_threshold),
        "NS"       = "NS"
      )
    ) +
    # Vertical LFC threshold lines
    ggplot2::geom_vline(
      xintercept = c(-lfc_threshold, lfc_threshold),
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +
    # Vertical center line at x = 0
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +
    # Horizontal significance line
    ggplot2::geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +
    # Zone count annotations
    ggplot2::annotate(
      "text",
      x = count_labels$x,
      y = count_labels$y,
      label = count_labels$label,
      color = count_labels$col,
      fontface = "bold",
      size = 4.5
    ) +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab,
      color = paste0("padj < ", padj_cutoff)
    ) +
    # ── cowplot-style theme ────────────────────────────────────────────────
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      # Title
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      # Axes: thick lines, outward ticks
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = ggplot2::unit(4, "pt"),
      # No grid at all
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      # Legend
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold")
    )

  # ── Apply axis limits if provided ─────────────────────────────────────────
  if (!is.null(xlim)) {
    p <- p + ggplot2::xlim(xlim)
  }
  if (!is.null(ylim)) {
    p <- p + ggplot2::ylim(ylim)
  }

  # ── Optional gene labels ─────────────────────────────────────────────────
  if (!is.null(label_genes) && gene_col %in% colnames(df)) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      warning(
        "Package 'ggrepel' is needed for gene labels. ",
        "Install it with: install.packages('ggrepel'). ",
        "Falling back to geom_text()."
      )
      label_df <- df[df[[gene_col]] %in% label_genes, ]
      p <- p + ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(label = .data[[gene_col]]),
        size = 3,
        vjust = -0.8,
        show.legend = FALSE
      )
    } else {
      label_df <- df[df[[gene_col]] %in% label_genes, ]
      p <- p + ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(label = .data[[gene_col]]),
        size = 3,
        max.overlaps = 20,
        show.legend = FALSE
      )
    }
  }

  # ── Print summary to console ──────────────────────────────────────────────
  message(
    sprintf("vizVolcano summary (padj < %g, |LFC| threshold = %g):", padj_cutoff, lfc_threshold),
    sprintf("\n  Down & |LFC| > %g : %d genes", lfc_threshold, n_down_high),
    sprintf("\n  Down & |LFC| <= %g: %d genes", lfc_threshold, n_down_low),
    sprintf("\n  Up   & |LFC| <= %g: %d genes", lfc_threshold, n_up_low),
    sprintf("\n  Up   & |LFC| > %g : %d genes", lfc_threshold, n_up_high),
    sprintf("\n  Non-significant    : %d genes", n_ns),
    sprintf("\n  Total plotted      : %d genes", nrow(df))
  )

  print(p)
  invisible(p)
}
