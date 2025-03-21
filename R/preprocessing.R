#' Create a SpatialExperiment Object
#'
#' This function creates a `SpatialExperiment` object by loading, filtering,
#' and processing image, cell, and panel data. It generates spatial coordinates,
#' and assigns colour palettes for visualizations.
#'
#' @param analysis_path A character string indicating the path to the analysis folder.
#'   Default is "analysis".
#' @param raw_path A character string indicating the path to the raw data folder.
#'   Default is "raw".
#' @return A `SpatialExperiment` object with processed assays, metadata, and
#'   spatial coordinates.
#' @export
#' @importFrom utils read.csv write.csv
#' @importFrom dplyr mutate select filter rename left_join group_by row_number distinct
#' @importFrom stringr str_remove str_extract str_replace
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData assay
#' @importFrom RColorBrewer brewer.pal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' spe <- CreateSpatialExperiment(
#'   analysis_path = "path/to/analysis",
#'   raw_path = "path/to/raw"
#' )
CreateSpatialExperiment <- function(analysis_path = "analysis", raw_path = "raw") {
    # Load CSVs
    cells <- utils::read.csv(
        file.path(analysis_path, "4_pyprofiler_output", "cell.csv")
    )
    panel <- utils::read.csv("panel.csv")
    image <- utils::read.csv("image.csv")

    # Filter rows and columns
    panel <- panel %>%
        dplyr::filter(.data$Full == 1)

    # Join image data with cell data
    dt <- dplyr::left_join(
        cells,
        image,
        by = dplyr::join_by("Image")
    )

    # Split data into marker intensities, metadata and co-ordinates
    markers <- panel[["Target"]]
    meta <- setdiff(colnames(dt), markers)

    counts <- dt[, markers]
    metadata <- dt[, setdiff(meta, c("X","Y"))] %>%
        # Change variables to 'factor' type
        dplyr::mutate(
            Image = as.factor(.data$Image),
            ImageID = as.factor(.data$ImageID),
            DonorID = as.factor(.data$DonorID)
        )
    coords <- dt[, c("X","Y")]

    # Create table of data for each marker
    rownames(panel) <- panel[["Target"]]
    markerData <- dplyr::select(panel, -"Target")
    markerData <- markerData[markers, ]

    # Create spatial object
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = t(counts)),
        colData = metadata,
        rowData = markerData,
        spatialCoords = as.matrix(coords),
        sample_id = as.character(metadata[["ImageID"]])
    )

    ### Assign colour palettes
    colour_vectors <- list()
    # DonorID
    donor_ids <- unique(spe[["DonorID"]])
    # Create a color vector that repeats if there are more than 12 unique DonorIDs
    donor_colors <- RColorBrewer::brewer.pal(12, name = "Paired")[
        as.numeric(donor_ids) %% 12 + 1
    ]
    # Set the names of the colors to match the unique DonorIDs
    DonorID <- stats::setNames(donor_colors, donor_ids)
    colour_vectors[["DonorID"]] <- DonorID

    # ImageID
    image_ids <- unique(spe[["ImageID"]])
    # Create a color vector that repeats if there are more than 12 unique ImageIDs
    image_colors <- RColorBrewer::brewer.pal(12, name = "Paired")[
        as.numeric(image_ids) %% 12 + 1
    ]
    # Set the names of the colors to match the unique ImageIDs
    ImageID <- stats::setNames(image_colors, image_ids)
    colour_vectors[["ImageID"]] <- ImageID

    # Merge colour palettes back into 'spe' object
    S4Vectors::metadata(spe)[["colour_vectors"]] <- colour_vectors
    colnames(spe) <- paste0(spe[["ImageID"]], "_", spe[["CellID"]])

    return(spe)
}

#' Normalise Data in a SpatialExperiment Object
#'
#' This function applies an arcsinh-transformation to the raw `"counts"` assay,
#' before standardising the data across channels within each image by
#' converting them to z-scores.
#'
#' @param object A `SpatialExperiment` object.
#' @return The input object with arcsinh-transformed data added as a new `"data"`
#'   assay, and z-standardised data added as a new `"scale.data"` assay.
#' @export
#' @importFrom SummarizedExperiment assay colData
#' @examples
#' spe <- NormaliseData(spe)
NormaliseData <- function(object) {
    # Apply asinh transformation to counts and store in "data" assay
    SummarizedExperiment::assay(object, "data") <- asinh(
        SummarizedExperiment::assay(object, "counts")
    )

    # Get unique images and channels
    images <- SummarizedExperiment::colData(object)$Image
    unique_images <- unique(images)
    exprs_matrix <- SummarizedExperiment::assay(object, "data")
    channels <- rownames(exprs_matrix)

    # Prepare a matrix to store standardised expressions
    base_matrix <- matrix(
        NA,
        nrow = nrow(exprs_matrix),
        ncol = ncol(exprs_matrix)
    )
    rownames(base_matrix) <- channels
    colnames(base_matrix) <- colnames(exprs_matrix)

    # Iterate over each image and channel to scale
    for (image in unique_images) {
        image_indices <- which(images == image)
        for (channel in channels) {
            channel_data <- exprs_matrix[channel, image_indices]
            scaled_data <- base::scale(channel_data)
            base_matrix[channel, image_indices] <- scaled_data
        }
    }

    # Assign the standardised data back to the assay
    SummarizedExperiment::assay(object, "scale.data") <- base_matrix

    return(object)
}
