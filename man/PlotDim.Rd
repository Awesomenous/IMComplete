% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{PlotDim}
\alias{PlotDim}
\title{Plot Dimensional Reduction Embeddings}
\usage{
PlotDim(
  object,
  reduction,
  prefixes = NULL,
  assay_name = "scale.data",
  var = NULL,
  subsample = NULL,
  dim = c(1, 2),
  seed = 123,
  ...
)
}
\arguments{
\item{object}{A `SpatialExperiment` object.}

\item{reduction}{A character string, either `"UMAP"` or `"PCA"`.}

\item{prefixes}{An optional vector of prefixes to filter clusters. Default is `NULL`.}

\item{assay_name}{A character string specifying the assay to use. Default is `"scale.data"`.}

\item{var}{A column in `colData` to color points by. Default is `NULL`.}

\item{subsample}{An optional integer specifying the number of cells to subsample. Default is `NULL`.}

\item{dim}{A numeric vector specifying which dimensions to plot. Default is `c(1, 2)`.}

\item{seed}{A numeric seed for reproducibility. Default is `123`.}

\item{...}{Additional arguments passed to `dittoDimPlot`.}
}
\value{
A `ggplot2` object.
}
\description{
Generates scatter plots for dimensional reduction embeddings such as PCA or UMAP.
}
\examples{
PlotDim(object = spe, reduction = "UMAP", var = "Cluster")
}
