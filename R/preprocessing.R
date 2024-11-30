#' Create Sorted Vector
#'
#' This function generates a sorted vector of integers from `1` to `n`.
#' It sorts the integers as characters using the `radix` method for
#' consistent lexical ordering.
#'
#' @param n An integer specifying the maximum value in the sequence to be sorted.
#' @return A sorted integer vector from `1` to `n`.
#' @examples
#' CreateSortedVector(10)
CreateSortedVector <- function(n) {
    char_numbers = as.character(1:n)
    sorted_char_numbers = sort(char_numbers, method = "radix")
    sorted_numbers = as.integer(sorted_char_numbers)
    return(sorted_numbers)
}

#' Create a SpatialExperiment Object
#'
#' This function creates a `SpatialExperiment` object by loading, filtering,
#' and processing image, cell, and panel data. It supports the removal of
#' specified metals, generates spatial coordinates, and assigns colour palettes
#' for visualizations.
#'
#' @param to_remove A character vector specifying metal tags to remove from the panel.
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
#' @examples
#' spe <- CreateSpatialExperiment(
#'   to_remove = c("Cd113", "Pd105"),
#'   analysis_path = "path/to/analysis",
#'   raw_path = "path/to/raw"
#' )
CreateSpatialExperiment <- function(to_remove,
                                    analysis_path = "analysis",
                                    raw_path = "raw") {
    # Load CSVs
    cells <- utils::read.csv(
        file.path(analysis_path, "5_cellprofiler_output", "cell.csv")
    )
    panel <- utils::read.csv(
        file.path(raw_path, "panel.csv")
    )
    image <- utils::read.csv(
        file.path(analysis_path, "5_cellprofiler_output", "image.csv")
    )

    # Filter rows and columns
    image <- image %>%
        dplyr::mutate(
            Image = stringr::str_remove(FileName_FullStack , "_full\\.tiff"),
            ROI = as.integer(
                stringr::str_extract(FileName_FullStack , "(?<=a)\\d+(?=_ac)")
            ) # <- NOT GENERALISED
        ) %>%
        dplyr::select(-FileName_FullStack)
    panel <- panel %>%
        dplyr::filter(trimws(Target) != "") %>%
        dplyr::rename(Metal = Metal.Tag) %>%
        dplyr::filter(!(Metal %in% to_remove))

    # Join image data with cell data
    cellsCombined <- dplyr::left_join(
        cells,
        image,
        by = dplyr::join_by(ImageNumber)
    )

    # Add 'ImageShort' column to dataframe <- NOT GENERALISED
    cellsCombined <- cellsCombined %>%
        dplyr::mutate(
            ImShort = stringr::str_extract(Image, "(?<=_)[^_]+_s0_a\\d+")
        ) %>%
        dplyr::mutate(
            ImShort = stringr::str_replace(ImShort, "_s0_a", "_")
        )

    # Define old column names and what to change them to
    rename_vec <- c(
        "Image" = "Image",
        "ImShort" = "ImShort",
        "ImageNumber" = "ImageID",
        "ROI" = "ROI",
        "ObjectNumber" = "CellID",
        "AreaShape_Area" = "Area",
        "AreaShape_Center_X" = "X",
        "AreaShape_Center_Y" = "Y",
        "DonorID" = "DonorID",
        "Condition" = "Condition"
    )

    # Change names in 'cellsCombined' from old names to new names
    names(cellsCombined) <- ifelse(
        names(cellsCombined) %in% names(rename_vec),
        rename_vec[names(cellsCombined)],
        names(cellsCombined)
    )

    # Name marker columns
    markers <- panel[,"Target"]
    ordered_markers <- markers[CreateSortedVector(length(markers))]
    colnames(cellsCombined)[6:(length(markers)+5)] <- ordered_markers

    # Keep only desired columns
    meta <- unname(sapply(strsplit(rename_vec, " = "), '[', 1))
    keep <- c(meta, ordered_markers)
    dt <- cellsCombined %>% dplyr::select(dplyr::all_of(keep))

    # Split data into marker intensities, metadata and co-ordinates
    counts <- dt[, ordered_markers] * 65535
    metadata <- dt[, setdiff(meta, c("X","Y"))]
    coords <- dt[, c("X","Y")]

    # Create table of data for each marker
    rownames(panel) <- panel$Target
    markerData <- dplyr::select(panel, -Target)
    markerData <- markerData[ordered_markers, ]

    # Create spatial object
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = t(counts)),
        colData = metadata,
        rowData = markerData,
        spatialCoords = as.matrix(coords),
        sample_id = as.character(metadata$ImageID)
    )

    # Change variables to 'factor' type
    spe$Image <- as.factor(spe$Image)
    spe$ImageID <- as.factor(spe$ImageID)
    spe$DonorID <- as.factor(spe$DonorID)

    # Add a universal CellID column
    SummarizedExperiment::colData(spe)$uCellID <- 1:length(spe$CellID)
    uIDKey <- as.data.frame(SummarizedExperiment::colData(spe))
    uIDKey <- uIDKey %>%
        dplyr::group_by(Image) %>%
        dplyr::filter(dplyr::row_number() == 1) %>%
        dplyr::select(Image, ImShort, ImageID, uCellID)

    # Write a CSV matching ImageIDs with uCellIDs
    utils::write.csv(
        uIDKey,
        file.path(analysis_path, "Image_uCellID_key.csv"),
        row.names = FALSE
    )

    ### Assign colour palettes
    color_vectors <- list()
    # DonorID
    donor_ids <- unique(spe$DonorID)
    # Create a color vector that repeats if there are more than 12 unique DonorIDs
    donor_colors <- RColorBrewer::brewer.pal(12, name = "Paired")[
        as.numeric(donor_ids) %% 12 + 1
    ]
    # Set the names of the colors to match the unique DonorIDs
    DonorID <- setNames(donor_colors, donor_ids)
    color_vectors$DonorID <- DonorID

    # ImageID
    image_ids <- unique(spe$ImageID)
    # Create a color vector that repeats if there are more than 12 unique ImageIDs
    image_colors <- RColorBrewer::brewer.pal(12, name = "Paired")[
        as.numeric(image_ids) %% 12 + 1
    ]
    # Set the names of the colors to match the unique ImageIDs
    ImageID <- stats::setNames(image_colors, image_ids)
    color_vectors$ImageID <- ImageID

    # Merge colour palettes back into 'spe' object
    S4Vectors::metadata(spe)$color_vectors <- color_vectors
    colnames(spe) <- paste0(spe$ImageID, "_", spe$CellID)

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
#'   assay, and z-standardised data added as a new `"standardised"` assay.
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
    SummarizedExperiment::assay(object, "standardised") <- base_matrix

    return(object)
}
