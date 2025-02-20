#' Batch Correction for SpatialExperiment Object
#'
#' This function performs batch correction on a `SpatialExperiment` object
#' using `scater` for PCA and `harmony` for batch correction. It splits the
#' data by a specified column (e.g., `DonorID`), runs PCA on selected markers,
#' and applies `harmony` to integrate data across batches.
#'
#' @param object A `SpatialExperiment` object
#'   containing the data to be batch-corrected.
#' @param ncomponents An integer specifying the number of principal components (PCs)
#'   to compute for batch correction.
#' @param col A character string specifying the column in `colData` used to
#'   define batch labels. Default is `"DonorID"`.
#' @param markers A character vector of marker names to use for PCA.
#'   Default is all markers (`rownames(object)`).
#' @param assay_name A character string specifying the assay to use for PCA.
#'   Default is `"data"`.
#' @param seed An integer seed for reproducibility. Default is 123.
#' @return The input object with a new reduced dimension named `"batch_corrected"`,
#'   containing the batch-corrected PCA embeddings.
#' @export
#' @importFrom scater runPCA
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom BiocSingular ExactParam
#' @importFrom harmony RunHarmony
#' @examples
#' # Perform batch correction using DonorID and specific markers
#' spe <- BatchCorrect(
#'   object = spe,
#'   ncomponents = 30,
#'   col = "DonorID",
#'   markers = c("Marker1", "Marker2")
#' )
BatchCorrect <- function(object,
                         ncomponents,
                         col = "DonorID",
                         markers = rownames(object),
                         assay_name = "data",
                         seed = 123) {
    # Stop if requesting too many PCs
    if (ncomponents > length(markers)) {
        stop("More components requested than number of markers available.")
    }

    # Run PCA
    object <- scater::runPCA(
        object,
        ncomponents = ncomponents,
        subset_row = markers,
        exprs_values = assay_name,
        BSPARAM = BiocSingular::ExactParam()
    )

    # Perform harmony batch correction
    set.seed(seed)
    out <- harmony::RunHarmony(
        object,
        group.by.vars = col
    )

    # Match column order to the original object
    embeddings <- SingleCellExperiment::reducedDim(out, "HARMONY")[colnames(object), ]

    # Add embeddings to the original object
    SingleCellExperiment::reducedDim(object, "batch_corrected") <- embeddings

    return(object)
}
