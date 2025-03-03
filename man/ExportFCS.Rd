% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{ExportFCS}
\alias{ExportFCS}
\title{Export Data to FCS Files}
\usage{
ExportFCS(
  object,
  extra_metadata = NULL,
  assay_name = "scale.data",
  export_path = "analysis/5_R_analysis/FlowJo"
)
}
\arguments{
\item{object}{A `SpatialExperiment` object.}

\item{extra_metadata}{A vector of character strings specifying additional metadata columns to include in the object.
Default is `NULL` (no additional columns).}

\item{assay_name}{A character string specifying the assay to export. Default
is `"scale.data"`.}

\item{export_path}{A character string specifying the directory to save the
exported FCS files. Default is `"FlowJo"`.}
}
\value{
None. FCS files are written to the specified directory.
}
\description{
This function exports assay data and metadata from a `SpatialExperiment`
object to individual FCS files, grouped by the `ImShort` column.
}
\examples{
ExportFCS(object = spe, assay_name = "scale.data", export_path = "FlowJo")
}
