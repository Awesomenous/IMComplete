#' Export Data to FCS Files
#'
#' This function exports assay data and metadata from a `SpatialExperiment`
#' object to individual FCS files, grouped by the `ImShort` column.
#'
#' @param object A `SpatialExperiment` object.
#' @param assay_name A character string specifying the assay to export. Default
#'   is `"standardised"`.
#' @param export_path A character string specifying the directory to save the
#'   exported FCS files. Default is `"FlowJo"`.
#' @param analysis_path A character string specifying the path to analysis data.
#'   Default is `"analysis"`.
#' @return None. FCS files are written to the specified directory.
#' @export
#' @importFrom utils read.csv txtProgressBar setTxtProgressBar
#' @importFrom dplyr select rename left_join filter mutate group_by ungroup
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom flowCore flowFrame write.FCS
#' @examples
#' ExportFCS(object = spe, assay_name = "standardised", export_path = "FlowJo")
ExportFCS <- function(object,
                      assay_name = "standardised",
                      export_path = "FlowJo",
                      analysis_path = "analysis") {

    # Read image data
    image <- utils::read.csv(
        file.path(analysis_path, "5_cellprofiler_output", "Image.csv")
    ) %>%
        dplyr::select(ImageNumber) %>%
        dplyr::rename(ImageID = ImageNumber)

    # Extract metadata
    metadata <- as.data.frame(SummarizedExperiment::colData(object)) %>%
        dplyr::select(ImageID, uCellID, CellID, Area, ImShort)

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
    df$ImageID <- as.numeric(df$ImageID)
    df <- df %>%
        dplyr::left_join(image, by = "ImageID") %>%
        dplyr::group_by(ImageID) %>%
        dplyr::mutate(Y = max(Y) - Y) %>%
        dplyr::ungroup()

    # Ensure export directory exists
    export_dir <- paste0(export_path, "/Export")
    dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

    # Export separate FCS files for each unique ImShort value
    unique_ImShort <- unique(df$ImShort)
    n_files <- length(unique_ImShort)

    # Initialize progress bar
    pb <- utils::txtProgressBar(min = 0, max = n_files, style = 3)

    for (i in seq_along(unique_ImShort)) {
        im_short <- unique_ImShort[i]
        subset_df <- df %>%
            dplyr::filter(ImShort == im_short) %>%
            dplyr::select(dplyr::any_of(colnames(assayData)))

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
#' @param assay_name A character string specifying the assay to use for TIFF generation.
#'   Default is `"standardised"`.
#' @param raw_path A character string specifying the path to raw data. Default is `"raw"`.
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
#' @examples
#' CreateMantisFolder(object = spe, raw_path = "raw", analysis_path = "analysis")
CreateMantisFolder <- function(object,
                               assay_name = "standardised",
                               raw_path = "raw",
                               analysis_path = "analysis",
                               mantis_path = "MantisProject") {

    # Create 'MantisProject' folder
    dir.create(mantis_path, showWarnings = FALSE)

    # List TIFFs and CSVs in the '1c_full_images' folder
    tiffs <- list.files(
        file.path(analysis_path, "1c_full_images"), pattern = "\\.tiff$"
    )
    csvs <- list.files(
        file.path(analysis_path, "1c_full_images"), pattern = "\\.csv$"
    )

    metadata <- as.data.frame(SummarizedExperiment::colData(object))

    # Progress bar for generating individual channel images
    pb1 <- utils::txtProgressBar(min = 0, max = length(tiffs), style = 3)
    for (x in seq_along(tiffs)) {
        # Extract the name of the current image from the filename
        cur_img_full <- gsub(".{10}$", "", tiffs[x])

        # Obtain corresponding 'ImShort' label
        cur_img <- metadata %>%
            dplyr::filter(Image == cur_img_full) %>%
            dplyr::distinct(ImShort) %>%
            dplyr::pull(ImShort)

        # Create a new folder for the current image
        dir.create(file.path(mantis_path, cur_img), showWarnings = FALSE)

        # Read TIFF and CSV associated with the current image
        cur_tiff <- suppressWarnings(
            tiff::readTIFF(
                file.path(analysis_path, "1c_full_images", tiffs[x]),
                all = TRUE
            )
        )
        cur_csv <- utils::read.csv(
            file.path(analysis_path, "1c_full_images", csvs[x]),
            header = FALSE
        )$V1

        # Save each individual channel in the TIFF as a separate image
        panel = utils::read.csv(file.path(raw_path, "panel.csv"))
        for (channel in seq_along(cur_tiff)) {
            metal_tag <- cur_csv[channel]
            marker <- panel$Target[panel$Metal.Tag == metal_tag]
            filename <- file.path(
                mantis_path,
                cur_img,
                paste0(marker, ".tiff")
            )
            mat <- cur_tiff[[channel]] * 65535
            ijtiff::write_tif(mat, filename, overwrite = TRUE, msg = FALSE)
        }
        utils::setTxtProgressBar(pb1, x)
    }
    close(pb1)

    # Progress bar for copying segmentation masks
    seg_list <- list.files(
        file.path(analysis_path, "3a_segmentation_masks"), pattern = "\\.tif$"
    )
    pb2 <- utils::txtProgressBar(min = 0, max = length(seg_list), style = 3)
    for (i in seq_along(seg_list)) {
        cur_dir <- seg_list[i]
        cur_img_full <- gsub(".{23}$", "", cur_dir)

        cur_img <- metadata %>%
            dplyr::filter(Image == cur_img_full) %>%
            dplyr::distinct(ImShort) %>%
            dplyr::pull(ImShort)

        new_dir <- file.path(mantis_path, cur_img, "SegmentationFile.tif")
        file.copy(
            file.path(analysis_path, "3a_segmentation_masks", cur_dir),
            new_dir
        )
        utils::setTxtProgressBar(pb2, i)
    }
    close(pb2)

    # Extract the full assay data (markers) and transpose it
    markers <- rownames(object)
    marker_data <- as.data.frame(
        t(SummarizedExperiment::assay(object, assay_name))
    )

    # Add metadata as columns to marker_data
    metadata <- as.data.frame(SummarizedExperiment::colData(object))
    marker_data <- cbind(metadata, marker_data)

    # Progress bar for writing segment features
    unique_imgs <- unique(SummarizedExperiment::colData(object)$ImShort)
    pb3 <- utils::txtProgressBar(min = 0, max = length(unique_imgs), style = 3)
    for (j in seq_along(unique_imgs)) {
        img <- unique_imgs[j]

        # Filter combined data to keep only rows corresponding to images
        filtered_data <- marker_data[marker_data$ImShort %in% c(img), ]

        # Select only relevant columns
        selected_columns <- c("ImShort", "CellID", markers)
        final_data <- filtered_data[, selected_columns]

        # Check if all entries in the 'ImShort' column are the same
        if (length(unique(final_data$ImShort)) == 1) {
            final_data$ImShort <- NULL
        }

        # Rename the remaining columns if necessary
        colnames(final_data)[
            colnames(final_data) %in% c("CellID")
        ] <- c("Segment ID")

        # Save the final filtered and selected data to a CSV file
        utils::write.csv(
            final_data,
            file.path(mantis_path, img, "SegmentFeatures.csv"),
            row.names = FALSE
        )
        utils::setTxtProgressBar(pb3, j)
    }
    close(pb3)
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
#' @examples
#' spe <- ImportFlowJoData(object = spe, FlowJo_path = "FlowJo")
ImportFlowJoData <- function(object, FlowJo_path) {
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
        df$FlowJo_celltype <- rep(celltype, nrow(df))
        return(df)
    }

    # Combine CSVs from FlowJo
    csvs <- lapply(csvs, add_celltype)
    flowjo_df <- dplyr::bind_rows(csvs) %>%
        dplyr::select(uCellID, FlowJo_celltype)
    flowjo_df$FlowJo_celltype <- as.factor(flowjo_df$FlowJo_celltype)

    # Extract existing metadata
    metadata <- as.data.frame(SummarizedExperiment::colData(object))

    # Remove old FlowJo annotations if they exist
    if ("FlowJo_celltype" %in% names(SummarizedExperiment::colData(object))) {
        metadata <- metadata %>%
            dplyr::select(-FlowJo_celltype)
    }

    # Merge FlowJo data with existing metadata
    flowjo_df <- dplyr::left_join(metadata, flowjo_df, by = "uCellID")

    # Remove duplicates by keeping the first label for each cell
    flowjo_df <- flowjo_df[!duplicated(flowjo_df$uCellID), ]

    # Add FlowJo annotations to the object
    SummarizedExperiment::colData(object)$FlowJo_celltype <-
        flowjo_df$FlowJo_celltype

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
