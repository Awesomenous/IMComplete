#' Cluster Cells Based on Marker Expression or Dimensional Reduction
#'
#' This function performs clustering on cells in a `SpatialExperiment` object
#' using either marker expression or a batch-corrected reduced dimension.
#'
#' @param object A `SpatialExperiment` object.
#' @param k An integer specifying the number of nearest neighbors to use for clustering.
#' @param prefix A character string prefix added to the front of cluster names
#'   (e.g., `prefix = "T cell"` will produce clusters named `"T cell_1"`, `"T cell_2`, etc.)
#' @param markers A character vector of markers to use for clustering. Default
#'   is all markers (`rownames(object)`).
#' @param filter_condition An optional condition for filtering cells prior to
#'   clustering, provided as a tidy evaluation formula (e.g., `celltype == "T cell"`).
#' @param assay_name A character string specifying the assay to use for clustering if
#'   `batch_corrected = FALSE`. Default is `"scale.data"`.
#' @param seed An integer seed for reproducibility. Default is 123.
#' @param batch_corrected A logical value. If `TRUE`, the `"batch_corrected"`
#'   reduced dimension is used instead of marker expression. Default is `FALSE`.
#' @return The input object with an updated `"Cluster"` column in its `colData`.
#'   If cells are already assigned to a cluster, that cluster identity will be
#'   overwritten with the new clustering results.
#' @export
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr filter select left_join mutate coalesce
#' @importFrom igraph membership
#' @importFrom Rphenograph Rphenograph
#' @importFrom rlang enquo
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' # Filter for T cells, perform clustering with 30 nearest neighbours, and name
#' # resulting clusters with the "T cell" prefix
#' spe <- ClusterCells(
#'     object = spe,
#'     k = 30,
#'     prefix = "T cell",
#'     filter_condition = celltype == "T cell"
#' )
#'
#' # Perform clustering of B cells using 20 nearest neighbours and
#' # batch-corrected reduced dimension instead of assay data
#' spe <- ClusterCells(
#'     object = spe,
#'     k = 20,
#'     prefix = "B cell",
#'     filter_condition = celltype == "B cell",
#'     batch_corrected = TRUE
#' )
ClusterCells <- function(object,
                         k,
                         prefix,
                         markers = rownames(object),
                         filter_condition = NULL,
                         assay_name = "scale.data",
                         seed = 123,
                         batch_corrected = FALSE) {

    # Make copy of original object
    object_ori <- object

    # Acquire metadata
    col_data <- as.data.frame(SummarizedExperiment::colData(object))

    # Filter object based on the filter_condition
    if (!missing(filter_condition)) {
        # Filter rows based on the condition within colData
        condition <- rlang::enquo(filter_condition)
        filtered_rows <- col_data %>%
            dplyr::filter(!!condition)

        if (nrow(filtered_rows) == 0) {
            stop("No rows match the filter condition.")
        }

        # Subset the object based on the filtered rows
        object <- object[, colnames(object) %in% rownames(filtered_rows)]
    }

    # Perform clustering
    if (batch_corrected) {
        if (!"batch_corrected" %in%
            SingleCellExperiment::reducedDimNames(object)) {
            stop("Please run BatchCorrect first using batch_corrected = TRUE.")
        }
        mat <- SingleCellExperiment::reducedDim(object, "batch_corrected")
    } else {
        mat <- t(SummarizedExperiment::assay(object, assay_name))[, markers]
    }

    set.seed(seed)
    out <- Rphenograph::Rphenograph(mat, k)

    # Generate cluster labels
    clusters <- factor(igraph::membership(out[[2]]))
    cluster_labels <- paste0(prefix, "_", clusters)

    # Prepare the new cluster data for merging
    new_clusters <- data.frame(
        uCellID = colnames(object),
        Cluster = as.factor(cluster_labels)
    )

    # Merge new clusters into the original object
    col_data_ori <- as.data.frame(SummarizedExperiment::colData(object_ori))
    col_data_ori$uCellID <- rownames(col_data_ori)

    if ("Cluster" %in% colnames(col_data_ori)) {
        # Update existing Cluster column
        col_data_ori <- col_data_ori %>%
            dplyr::left_join(
                new_clusters, by = "uCellID", suffix = c("_old", "_new")
            ) %>%
            dplyr::mutate(
                Cluster = dplyr::coalesce(.data$Cluster_new, .data$Cluster_old)
            ) %>%
            dplyr::select(-"Cluster_old", -"Cluster_new")
    } else {
        # Add new Cluster column
        col_data_ori <- col_data_ori %>%
            dplyr::left_join(new_clusters, by = "uCellID")
    }

    # Update colData in the original object
    SummarizedExperiment::colData(object_ori) <- col_data_ori %>%
        dplyr::select(-"uCellID") %>%
        S4Vectors::DataFrame()

    # Return the updated object
    return(object_ori)
}

#' Annotate Clusters with User-Defined Labels
#'
#' This function assigns user-defined labels to clusters in a
#' `SpatialExperiment` object based on a mapping.
#'
#' @param object A `SpatialExperiment` object.
#' @param label_mapping A named vector where names correspond to existing cluster
#'   labels and values correspond to the new labels.
#' @return The input object with a new `"LabelledCluster"` column in its `colData`.
#' @export
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr recode mutate
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' # Annotate clusters
#' label_map <- c(
#'     "T cell_1" = "Treg",
#'     "T cell_2" = "CD4+ T cell",
#'     "T cell_3" = "CD8+ T cell"
#' )
#' spe <- AnnotateClusters(object = spe, label_mapping = label_map)
AnnotateClusters <- function(object, label_mapping) {
    # Recode the Cluster column based on the provided mapping
    recoded <- dplyr::recode(
        SummarizedExperiment::colData(object)[["Cluster"]],
        !!!label_mapping
    )

    # Extract and modify metadata
    metadata <- SummarizedExperiment::colData(object) %>%
        as.data.frame()
    metadata[["LabelledCluster"]] <- recoded

    # Replace matching labels with NA
    metadata <- metadata %>%
        dplyr::mutate(
            LabelledCluster = dplyr::if_else(
                as.character(.data$LabelledCluster) == as.character(.data$Cluster),
                NA_character_,
                .data$LabelledCluster
            )
        )

    # Add the new LabelledCluster column back to the object
    SummarizedExperiment::colData(object)[["LabelledCluster"]] <-
        metadata[["LabelledCluster"]]

    return(object)
}
