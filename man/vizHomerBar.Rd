\name{vizHomerBar}
\alias{vizHomerBar}
\title{Visualize HOMER Known Motif Results as a ComplexHeatmap Bar Plot}
\description{
Creates a publication-quality visualization of HOMER known motif discovery
results using \pkg{ComplexHeatmap}. The plot contains five columns:
\enumerate{
  \item Rank of the motifs
  \item Transcription factor names
  \item Motif logos (PNG images)
  \item Horizontal bar plot of -log(P-value)
  \item Horizontal bar plot of percent of targets
}
}
\usage{
vizHomerBar(homer_data, top_n = 20, pvalue_col = "#2166AC",
            target_col = "#B2182B",
            title = "HOMER Known Motif Enrichment",
            logo_width = grid::unit(2, "cm"), bar_width = grid::unit(4, "cm"),
            row_height = grid::unit(8, "mm"), ...)
}
\arguments{
  \item{homer_data}{A \code{\link{HomerData-class}} object.}
  \item{top_n}{Integer or NULL. Number of top motifs to display. If NULL, all
    motifs are shown. Default is 20.}
  \item{pvalue_col}{Character. Color for the -log(P-value) bar plot.
    Default is \code{"#2166AC"} (blue).}
  \item{target_col}{Character. Color for the target percentage bar plot.
    Default is \code{"#B2182B"} (red).}
  \item{title}{Character. Title for the plot. Default is
    \code{"HOMER Known Motif Enrichment"}.}
  \item{logo_width}{\code{grid::unit} object. Width of the motif logo column.
    Default is \code{unit(2, "cm")}.}
  \item{bar_width}{\code{grid::unit} object. Width of each bar plot column.
    Default is \code{unit(4, "cm")}.}
  \item{row_height}{\code{grid::unit} object. Height of each row.
    Default is \code{unit(8, "mm")}.}
  \item{...}{Additional arguments passed to \code{ComplexHeatmap::draw()}.}
}
\value{
A \code{ComplexHeatmap::Heatmap} object (invisibly). The plot is drawn as a
side effect.
}
\seealso{
  \code{\link{HomerData-class}}, \code{\link{readHomerDir}}
}
\examples{
\dontrun{
homer_data <- readHomerDir("/path/to/homer/output")
vizHomerBar(homer_data, top_n = 15)

# Customize colors
vizHomerBar(homer_data,
            top_n = 10,
            pvalue_col = "#4393C3",
            target_col = "#D6604D",
            title = "My Motif Enrichment")
}
}
