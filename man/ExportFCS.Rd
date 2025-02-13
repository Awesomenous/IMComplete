% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{ExportFCS}
\alias{ExportFCS}
\title{Export Data to FCS Files}
\usage{
ExportFCS(
  object,
  assay_name = "standardised",
  export_path = "FlowJo",
  analysis_path = "analysis"
)
}
\arguments{
\item{object}{A `SpatialExperiment` object.}

\item{assay_name}{A character string specifying the assay to export. Default
is `"standardised"`.}

\item{export_path}{A character string specifying the directory to save the
exported FCS files. Default is `"FlowJo"`.}

\item{analysis_path}{A character string specifying the path to analysis data.
Default is `"analysis"`.}
}
\value{
None. FCS files are written to the specified directory.
}
\description{
This function exports assay data and metadata from a `SpatialExperiment`
object to individual FCS files, grouped by the `ImShort` column.
}
\examples{
ExportFCS(object = spe, assay_name = "standardised", export_path = "FlowJo")
}
