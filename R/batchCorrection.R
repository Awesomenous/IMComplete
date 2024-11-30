#' Batch Correction for SpatialExperiment Object
#'
#' This function performs batch correction on a `SpatialExperiment` object
#' using the `Seurat` package. It splits the data by a specified column
#' (e.g., `DonorID`), performs PCA on each batch, finds integration anchors,
#' and integrates the data across batches.
#'
#' @param object A `SpatialExperiment` object
#'   containing the data to be batch-corrected.
#' @param col A character string specifying the column in `colData` used to
#'   split the data into batches. Default is `"DonorID"`.
#' @param markers A character vector of marker names to use for batch correction.
#'   Default is all markers (`rownames(object)`).
#' @param k An integer specifying the number of nearest neighbors to use when
#'   finding integration anchors. Default is 20.
#' @return The input object with a new reduced dimension named `"batch_corrected"`,
#'   containing the batch-corrected PCA embeddings.
#' @export
#' @importFrom Seurat as.Seurat AddMetaData ScaleData RunPCA SplitObject
#' @importFrom Seurat FindIntegrationAnchors IntegrateData DefaultAssay Embeddings
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData
#' @examples
#' # Perform batch correction using DonorID and specific markers
#' spe <- BatchCorrect(
#'   object = spe,
#'   col = "DonorID",
#'   markers = c("Marker1", "Marker2"),
#'   k = 30
#' )
BatchCorrect <- function(object,
                         col = "DonorID",
                         markers = rownames(object),
                         k = 20) {
    # Rename columns to ensure unique identifiers
    colnames(object) <- paste0("cell", object$uCellID)

    # Convert to Seurat object
    seurat_obj <- suppressWarnings(
        Seurat::as.Seurat(
            object,
            counts = "counts",
            data = "data"
        )
    )
    seurat_obj <- Seurat::AddMetaData(
        seurat_obj, as.data.frame(SummarizedExperiment::colData(object))
    )

    # Split Seurat object by 'col'
    seurat_list <- Seurat::SplitObject(seurat_obj, split.by = col)

    # Scale data and run PCA for each Seurat object
    seurat_list <- lapply(seurat_list, function(x) {
        x <- Seurat::ScaleData(x, features = markers, verbose = FALSE)
        x <- Seurat::RunPCA(
            x, features = markers, verbose = FALSE, approx = FALSE
        )
        return(x)
    })

    # Find integration anchors
    anchors <- Seurat::FindIntegrationAnchors(
        object.list = seurat_list,
        anchor.features = markers,
        reduction = "rpca",
        k.anchor = k
    )

    # Integrate data
    combined <- Seurat::IntegrateData(anchorset = anchors)
    Seurat::DefaultAssay(combined) <- "integrated"

    # Scale and run PCA on integrated data
    combined <- Seurat::ScaleData(combined, verbose = FALSE)
    combined <- Seurat::RunPCA(
        combined,
        npcs = min(30, length(markers)),
        verbose = FALSE,
        approx = FALSE
    )

    # Reorder columns if necessary
    if (!all.equal(colnames(object), colnames(combined))) {
        # Match column order to the original object
        combined <- combined[, match(colnames(object), colnames(combined))]
    }

    # Add PCA embeddings to the original object
    SingleCellExperiment::reducedDim(object, "batch_corrected") <-
        Seurat::Embeddings(combined, reduction = "pca")

    return(object)
}
