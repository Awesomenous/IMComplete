% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{PlotViolin}
\alias{PlotViolin}
\title{Plot Violin Plots for Features}
\usage{
PlotViolin(
  object,
  assay_name,
  features,
  split_by = NULL,
  ncol = 1,
  jitter = TRUE,
  pt_size = 0.01,
  alpha = 0.5,
  geom_violin_args = list(),
  geom_jitter_args = list()
)
}
\arguments{
\item{object}{A `SpatialExperiment` object.}

\item{assay_name}{A character string specifying the assay to use for plotting.}

\item{features}{A character vector of feature names to plot.}

\item{split_by}{An optional column name in `colData` to split the violins by.
Default is `NULL`.}

\item{ncol}{Number of columns for arranging plots in the grid. Default is `1`.}

\item{jitter}{Logical, whether to add jitter points. Default is `TRUE`.}

\item{pt_size}{Point size for jitter. Default is `0.01`.}

\item{alpha}{Alpha transparency for jitter points. Default is `0.5`.}

\item{geom_violin_args}{A list of additional arguments passed to `geom_violin`.}

\item{geom_jitter_args}{A list of additional arguments passed to `geom_jitter`.}
}
\value{
A `patchwork`ed `ggplot2` object.
}
\description{
Creates violin plots for specified features from a `SpatialExperiment` object.
Allows customization of aesthetics and jitter options.
}
\examples{
PlotViolin(object = spe, assay_name = "scale.data", features = c("Gene1", "Gene2"))
}
