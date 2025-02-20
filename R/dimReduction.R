#' Perform Principal Component Analysis (PCA)
#'
#' This function performs PCA on a `SpatialExperiment` object.
#' PCA can be performed on all cells or a subset defined by prefixes.
#'
#' @param object A `SpatialExperiment` object.
#' @param ncomponents An integer specifying the number of principal components
#'   to compute.
#' @param prefixes A character vector of prefixes used to subset the cells for
#'   PCA. Default is `NULL` (all cells are used).
#' @param markers A character vector of marker names to include in the PCA.
#'   Default is all markers (`rownames(object)`).
#' @param assay_name A character string specifying the assay to use for PCA.
#'   Default is `"scale.data"`.
#' @param seed An integer seed for reproducibility. Default is 123.
#' @param ... Additional arguments passed to `scater::runPCA`.
#' @return The input object with PCA results stored in reduced dimensions:
#'   - `"PCA_full"`: PCA on all cells.
#'   - `"PCA_subset"`: PCA on a subset of cells (if `prefixes` are provided).
#' @export
#' @importFrom SummarizedExperiment colData assay
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom scater runPCA
#' @importFrom S4Vectors metadata
#' @importFrom BiocSingular ExactParam
#' @examples
#' # Perform PCA on all cells
#' spe <- PerformPCA(object = spe, ncomponents = 10)
#'
#' # Perform PCA on a subset of cells
#' spe <- PerformPCA(object = spe, ncomponents = 10, prefixes = c("B cell", "T cell"))
PerformPCA <- function(object,
                       ncomponents,
                       prefixes = NULL,
                       markers = rownames(object),
                       assay_name = "scale.data",
                       seed = 123,
                       ...) {
    set.seed(seed)

    # Check if requested components exceed the number of markers
    if (ncomponents > length(markers)) {
        stop("More components requested than number of markers available.")
    }

    if (is.null(prefixes)) {
        # Perform PCA on all cells
        object <- scater::runPCA(
            object,
            ncomponents = ncomponents,
            subset_row = markers,
            exprs_values = assay_name,
            BSPARAM = BiocSingular::ExactParam(),
            ...
        )
        SingleCellExperiment::reducedDim(object, "PCA_full") <-
            SingleCellExperiment::reducedDim(object, "PCA")
    } else {
        # Perform PCA on subset of cells with the specified prefixes
        cluster_annotation <- SummarizedExperiment::colData(object)$Cluster
        is_selected <- sapply(
            prefixes,
            function(prefix) grepl(paste0("^", prefix), cluster_annotation)
        )
        selected_cells <- apply(is_selected, 1, any)

        # Subset the object for PCA
        subset_object <- object[, selected_cells]
        subset_object <- scater::runPCA(
            subset_object,
            ncomponents = ncomponents,
            subset_row = markers,
            exprs_values = assay_name,
            ...
        )

        # Initialize PCA_subset with NA
        pca_subset <- matrix(NA, nrow = ncol(object), ncol = ncomponents)
        pca_subset[selected_cells, ] <-
            SingleCellExperiment::reducedDim(subset_object, "PCA")

        # Add the PCA_subset to the full object
        SingleCellExperiment::reducedDim(object, "PCA_subset") <- pca_subset

        # Record the prefixes used for PCA_subset calculation
        S4Vectors::metadata(object)$PCA_subset_prefixes <- prefixes
    }

    return(object)
}

#' Perform Uniform Manifold Approximation and Projection (UMAP)
#'
#' This function performs UMAP dimensionality reduction on a
#' `SpatialExperiment` object. UMAP can be performed on all cells or a subset
#' defined by prefixes.
#'
#' @param object A `SpatialExperiment` object.
#' @param prefixes A character vector of prefixes used to subset the cells for
#'   UMAP. Default is `NULL` (all cells are used).
#' @param markers A character vector of marker names to include in the UMAP.
#'   Default is all markers (`rownames(object)`).
#' @param assay_name A character string specifying the assay to use for UMAP.
#'   Default is `"scale.data"`.
#' @param seed An integer seed for reproducibility. Default is 123.
#' @param batch_corrected A logical value. If `TRUE`, the `"batch_corrected"`
#'   reduced dimension is used instead of the assay. Default is `FALSE`.
#' @param ... Additional arguments passed to `scater::runUMAP`.
#' @return The input object with UMAP results stored in reduced dimensions:
#'   - `"UMAP_full"`: UMAP on all cells.
#'   - `"UMAP_subset"`: UMAP on a subset of cells (if `prefixes` are provided).
#' @export
#' @importFrom SummarizedExperiment colData assay
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom scater runUMAP
#' @importFrom S4Vectors metadata
#' @importFrom dplyr filter
#' @examples
#' # Perform UMAP on all cells
#' spe <- PerformUMAP(object = spe, markers = c("CD3", "CD4"))
#'
#' # Perform UMAP on batch-corrected data for a subset of cells
#' spe <- PerformUMAP(object = spe, prefixes = c("B", "T"), batch_corrected = TRUE)
PerformUMAP <- function(object,
                        prefixes = NULL,
                        markers = rownames(object),
                        assay_name = "scale.data",
                        seed = 123,
                        batch_corrected = FALSE,
                        ...) {
    set.seed(seed)

    # Check for "batch_corrected" if batch_corrected is TRUE
    if (batch_corrected) {
        if (!"batch_corrected" %in%
            SingleCellExperiment::reducedDimNames(object)) {
            stop("Please run BatchCorrect before using batch_corrected = TRUE.")
        }
    }

    if (is.null(prefixes)) {
        # Perform UMAP on all cells
        if (batch_corrected) {
            object <- scater::runUMAP(
                object,
                dimred = "batch_corrected",
                subset_row = markers,
                name = "UMAP_full",
                ...
            )
        } else {
            object <- scater::runUMAP(
                object,
                exprs_values = assay_name,
                subset_row = markers,
                name = "UMAP_full",
                ...
            )
        }
    } else {
        # Perform UMAP on subset of cells with the specified prefixes
        cluster_annotation <- SummarizedExperiment::colData(object)$Cluster
        is_selected <- sapply(
            prefixes,
            function(prefix) grepl(paste0("^", prefix), cluster_annotation)
        )
        selected_cells <- apply(is_selected, 1, any)

        # Subset the object
        subset_object <- object[, selected_cells]

        # Perform UMAP on the subset
        if (batch_corrected) {
            subset_object <- scater::runUMAP(
                subset_object,
                dimred = "batch_corrected",
                subset_row = markers,
                ...
            )
        } else {
            subset_object <- scater::runUMAP(
                subset_object,
                subset_row = markers,
                exprs_values = assay_name,
                ...
            )
        }

        # Initialize UMAP_subset with NA
        umap_subset <- matrix(NA, nrow = ncol(object), ncol = 2)
        umap_subset[selected_cells, ] <-
            SingleCellExperiment::reducedDim(subset_object, "UMAP")

        # Add the UMAP_subset to the full object
        SingleCellExperiment::reducedDim(object, "UMAP_subset") <- umap_subset

        # Record the prefixes used for UMAP_subset calculation
        S4Vectors::metadata(object)$UMAP_subset_prefixes <- prefixes
    }

    return(object)
}
