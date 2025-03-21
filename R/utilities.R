#' Export Data to FCS Files
#'
#' This function exports assay data and metadata from a `SpatialExperiment`
#' object to individual FCS files, grouped by the `ImShort` column.
#'
#' @param object A `SpatialExperiment` object.
#' @param extra_metadata A vector of character strings specifying additional metadata columns to include in the object.
#'  Default is `NULL` (no additional columns).
#' @param assay_name A character string specifying the assay to export. Default
#'   is `"scale.data"`.
#' @param export_path A character string specifying the directory to save the
#'   exported FCS files. Default is `"FlowJo"`.
#' @return None. FCS files are written to the specified directory.
#' @export
#' @importFrom utils read.csv txtProgressBar setTxtProgressBar
#' @importFrom dplyr select rename left_join filter mutate group_by ungroup
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom flowCore flowFrame write.FCS
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' ExportFCS(object = spe, assay_name = "scale.data", export_path = "FlowJo")
ExportFCS <- function(object,
                      extra_metadata = NULL,
                      assay_name = "scale.data",
                      export_path = "analysis/5_R_analysis/FlowJo") {

    # Extract metadata
    metadata_list <- c("ImageID", "CellID", "Area", "ImShort")
    if (!is.null(extra_metadata)) {
        metadata_list <- unique(c(metadata_list, extra_metadata[extra_metadata %in% colnames(SummarizedExperiment::colData(object))])) # COME BACK HERE AND CHECK IF CORRECT
    }

    metadata <- as.data.frame(SummarizedExperiment::colData(object)) %>%
        dplyr::select(any_of(metadata_list))

    # CHECK THAT EXTRA_METADATA IS NUMERIC

    # Extract assay data
    assayData <- as.data.frame(
        t(SummarizedExperiment::assay(object, assay_name))
    )

    # Combine metadata and assay data
    df <- cbind(
        metadata,
        as.data.frame(SpatialExperiment::spatialCoords(object)),
        assayData
    )
    df$ImageID <- as.numeric(df[["ImageID"]])
    df <- df %>%
        dplyr::group_by(.data$ImageID) %>%
        dplyr::mutate(Y = max(.data$Y) - .data$Y) %>% # CHECK IF IMAGE IS INVERTED
        dplyr::ungroup()

    # Ensure export directory exists
    export_dir <- paste0(export_path, "/Export")
    dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

    # Export separate FCS files for each unique ImShort value
    unique_ImShort <- unique(df[["ImShort"]])
    n_files <- length(unique_ImShort)

    # Initialize progress bar
    pb <- utils::txtProgressBar(min = 0, max = n_files, style = 3)

    for (i in seq_along(unique_ImShort)) {
        im_short <- unique_ImShort[i]
        subset_df <- df %>%
            dplyr::filter(.data$ImShort == im_short) %>%
            dplyr::select(dplyr::any_of(
                setdiff(colnames(df), "ImShort")
            ))

        # Convert to matrix and ensure numeric
        data_mat <- apply(as.matrix(subset_df), 2, as.numeric)
        colnames(data_mat) <- make.names(colnames(subset_df))

        # Create flowFrame
        ff <- flowCore::flowFrame(data_mat)

        # Write FCS file with ImShort value in the name
        fcs_file_name <- paste0(export_path, "/", im_short, ".fcs")
        invisible(flowCore::write.FCS(ff, file = fcs_file_name))

        # Update progress bar
        utils::setTxtProgressBar(pb, i)
    }

    # Close progress bar
    close(pb)
}

#' Create Mantis Project Folder with Data
#'
#' This function generates a folder structure compatible with Mantis, including
#' TIFF images, segmentation masks, and metadata.
#'
#' @param object A `SpatialExperiment` object.
#' @param img_suffix A character string specifying the suffix added to the end of the full image names.
#'  Default is `"_full.tiff"`.
#' @param mask_suffix A character string specifying the suffix added to the end of the mask image names.
#'  Default is `"_CpSeg_mask.tif"`.
#' @param assay_name A character string specifying the assay to use for TIFF generation.
#'   Default is `"scale.data"`.
#' @param analysis_path A character string specifying the path to analysis data.
#'   Default is `"analysis"`.
#' @param mantis_path A character string specifying the path for the Mantis project.
#'   Default is `"MantisProject"`.
#' @return None. Files are written to the specified Mantis folder.
#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @importFrom SummarizedExperiment assay colData
#' @importFrom ijtiff write_tif
#' @importFrom dplyr filter distinct pull select bind_rows
#' @importFrom tiff readTIFF
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' CreateMantisFolder(object = spe, raw_path = "raw", analysis_path = "analysis")
CreateMantisFolder <- function(object,
                               img_suffix = "_full.tiff",
                               mask_suffix = "_CpSeg_mask.tif",
                               assay_name = "scale.data",
                               analysis_path = "analysis",
                               mantis_path = "analysis/5_R_analysis/MantisProject") {

    # Create 'MantisProject' folder
    dir.create(mantis_path, showWarnings = FALSE)

    # Load image names
    image_df <- utils::read.csv("image.csv")
    image_list <- image_df[["Image"]]
    imshort_list <- image_df[["ImShort"]]

    # Obtain panel
    panel <- utils::read.csv("panel.csv") %>%
        dplyr::filter(.data$Full == 1)

    # Extract the full assay data (markers) and transpose it
    markers <- rownames(object)
    marker_data <- as.data.frame(
        t(SummarizedExperiment::assay(object, assay_name))
    )

    # Add metadata as columns to marker_data
    metadata <- as.data.frame(SummarizedExperiment::colData(object))
    marker_data <- cbind(metadata, marker_data)

    # Progress bar for generating individual channel images
    pb <- utils::txtProgressBar(min = 0, max = length(image_list), style = 3)
    for (x in seq_along(image_list)) {
        # Create a new folder for the current image
        dir.create(file.path(mantis_path, imshort_list[x]), showWarnings = FALSE)

        # Read TIFF and CSV associated with the current image
        cur_tiff <- suppressWarnings(
            tiff::readTIFF(
                file.path(analysis_path, "2_cleaned", paste0(image_list[x], img_suffix)),
                all = TRUE
            )
        )

        # Save each individual channel in the TIFF as a separate image
        for (channel in seq_along(cur_tiff)) {
            marker <- panel$Target[channel]
            filename <- file.path(
                mantis_path,
                imshort_list[x],
                paste0(marker, ".tiff")
            )
            mat <- cur_tiff[[channel]] * 65535
            ijtiff::write_tif(mat, filename, overwrite = TRUE, msg = FALSE)
        }

        # Copy segmentation mask
        file.copy(
            file.path(analysis_path, "3_segmentation", "3c_cellpose_mask", paste0(image_list[x], mask_suffix)),
            file.path(mantis_path, imshort_list[x], "SegmentationFile.tif")
        )

        # Filter combined data to keep only rows corresponding to image
        filtered_data <- marker_data[marker_data$ImShort == imshort_list[x], ]

        # Select only relevant columns
        selected_columns <- c("CellID", markers)
        final_data <- filtered_data[, selected_columns]

        # Rename the remaining columns if necessary
        colnames(final_data)[
            colnames(final_data) %in% c("CellID")
        ] <- c("Segment ID")

        # Save the final filtered and selected data to a CSV file
        utils::write.csv(
            final_data,
            file.path(mantis_path, imshort_list[x], "SegmentFeatures.csv"),
            row.names = FALSE
        )

        utils::setTxtProgressBar(pb, x)
    }
    close(pb)
}

#' Import Data from FlowJo CSV Files
#'
#' This function imports FlowJo annotations into a `SpatialExperiment` object
#' and updates the metadata.
#'
#' @param object A `SpatialExperiment` object.
#' @param FlowJo_path A character string specifying the path to the FlowJo export directory.
#' @return The input object with updated `FlowJo_celltype` annotations in its `colData`.
#' @export
#' @importFrom utils read.csv
#' @importFrom dplyr select bind_rows left_join
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' spe <- ImportFlowJoData(object = spe, FlowJo_path = "FlowJo")
ImportFlowJoData <- function(object, FlowJo_path = "analysis/5_R_analysis/FlowJo") {
    # List exported CSVs from FlowJo
    csvs <- list.files(
        path = file.path(FlowJo_path, "Export"),
        pattern = "*.csv",
        full.names = TRUE
    )

    # Helper function to read and annotate each CSV with its cell type
    add_celltype <- function(x) {
        df <- utils::read.csv(x, header = TRUE)
        celltype <- sub("\\.csv$", "", basename(x))
        df[["FlowJo_celltype"]] <- rep(celltype, nrow(df))
        return(df)
    }

    # Combine CSVs from FlowJo
    csvs <- lapply(csvs, add_celltype)
    flowjo_df <- dplyr::bind_rows(csvs) %>%
        dplyr::select("ImageID", "CellID", "FlowJo_celltype") %>%
        dplyr::mutate(uCellID = paste0(.data$ImageID, "_", .data$CellID))
    flowjo_df[["FlowJo_celltype"]] <- as.factor(flowjo_df[["FlowJo_celltype"]])

    # Extract existing metadata
    metadata <- as.data.frame(SummarizedExperiment::colData(object)) %>%
        dplyr::mutate(uCellID = paste0(.data$ImageID, "_", .data$CellID))

    # Remove old FlowJo annotations if they exist
    if ("FlowJo_celltype" %in% names(SummarizedExperiment::colData(object))) {
        metadata <- metadata %>%
            dplyr::select(-"FlowJo_celltype")
    }

    # Merge FlowJo data with existing metadata
    flowjo_df <- dplyr::left_join(metadata, flowjo_df, by = "uCellID")

    # Remove duplicates by keeping the first label for each cell
    num_duplicates <- sum(duplicated(flowjo_df[["uCellID"]]))
    message(paste(
        num_duplicates, "cells found in multiple gates. First instance chosen."
    ))
    flowjo_df <- flowjo_df[!duplicated(flowjo_df[["uCellID"]]), ]

    # Add FlowJo annotations to the object
    SummarizedExperiment::colData(object)[["FlowJo_celltype"]] <-
        flowjo_df[["FlowJo_celltype"]]

    return(object)
}

#' Export Metadata to Mantis-Compatible CSV
#'
#' This function extracts and exports metadata from a `SpatialExperiment` object
#' to a CSV file compatible with Mantis.
#'
#' @param object A `SpatialExperiment` object.
#' @param column A character string specifying the column in `colData` to export.
#' @param file_path A character string specifying the file path for the output CSV.
#' @return None. A CSV file is written to the specified path.
#' @export
#' @importFrom utils write.table
#' @importFrom dplyr select any_of
#' @importFrom stats na.omit
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr %>%
#' @examples
#' ExportToMantis(object = spe, column = "Cluster", file_path = "mantis_export.csv")
ExportToMantis <- function(object, column, file_path) {
    # Extract relevant data from colData
    to_export <- SummarizedExperiment::colData(object) %>%
        as.data.frame() %>%
        dplyr::select(dplyr::any_of(c("ImShort", "CellID", column))) %>%
        stats::na.omit()

    # Split data by ImShort
    export_by_image <- split(to_export, to_export$ImShort)

    # Function to extract cell types for each image
    extract_celltypes <- function(img) {
        df <- export_by_image[[img]]
        return(as.data.frame(df))
    }

    # Apply the extraction function and combine results
    export <- lapply(names(export_by_image), extract_celltypes) %>%
        dplyr::bind_rows()

    # Write the combined data to a file
    utils::write.table(
        export,
        file = file_path,
        sep = ",",
        col.names = FALSE,
        row.names = FALSE
    )
}

#' Retrieve a Set of Markers for Analysis
#'
#' This function extracts a set of markers from panel.csv based on a specified column.
#' The extracted markers can then be used as input for other functions in downstream analyses.
#'
#' @param colname A character string specifying the column in panel.csv to filter by.
#'                Default is "Full".
#' @return A character vector containing the names of the selected markers.
#' @export
#' @importFrom utils read.csv
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @examples
#' # Extract markers from the "Full" column
#' markers <- UseMarkers(colname = "Full")
#'
#' # Extract markers based on a different column
#' markers_subset <- UseMarkers(colname = "Subset")
UseMarkers <- function(colname = "Full") {
    # Load panel data
    panel <- utils::read.csv("panel.csv")

    # Filter panel data based on given 'colname'
    panel <- panel %>%
        dplyr::filter(.data[[colname]] == 1)

    # Extract markers
    markers <- panel[["Target"]]

    # Return markers
    return(markers)
}


