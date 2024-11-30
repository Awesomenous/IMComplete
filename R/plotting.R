#' Plot Violin Plots for Features
#'
#' Creates violin plots for specified features from a `SpatialExperiment` object.
#' Allows customization of aesthetics and jitter options.
#'
#' @param object A `SpatialExperiment` object.
#' @param assay_name A character string specifying the assay to use for plotting.
#' @param features A character vector of feature names to plot.
#' @param split_by An optional column name in `colData` to split the violins by.
#'   Default is `NULL`.
#' @param ncol Number of columns for arranging plots in the grid. Default is `1`.
#' @param jitter Logical, whether to add jitter points. Default is `TRUE`.
#' @param pt_size Point size for jitter. Default is `0.01`.
#' @param alpha Alpha transparency for jitter points. Default is `0.5`.
#' @param geom_violin_args A list of additional arguments passed to `geom_violin`.
#' @param geom_jitter_args A list of additional arguments passed to `geom_jitter`.
#' @return A `patchwork`ed `ggplot2` object.
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_violin geom_jitter theme_classic theme element_text element_blank labs
#' @importFrom dplyr filter
#' @importFrom SummarizedExperiment colData assay
#' @importFrom patchwork wrap_plots
#' @examples
#' PlotViolin(object = spe, assay_name = "standardised", features = c("Gene1", "Gene2"))
PlotViolin <- function(object,
                       assay_name,
                       features,
                       split_by = NULL,
                       ncol = 1,
                       jitter = TRUE,
                       pt_size = 0.01,
                       alpha = 0.5,
                       geom_violin_args = list(),
                       geom_jitter_args = list()) {

    # Extract metadata
    metadata <- SummarizedExperiment::colData(object) %>%
        as.data.frame()

    # Extract assay data and combine with metadata
    data <- SummarizedExperiment::assay(object, assay_name) %>%
        t() %>%
        as.data.frame() %>%
        cbind(metadata)

    # Create a list of violin plots
    plot_list <- lapply(features, function(feature) {
        if (!is.null(split_by)) {
            p <- data %>%
                dplyr::filter(!is.na(.data[[split_by]])) %>%
                ggplot2::ggplot(
                    ggplot2::aes_string(x = split_by, y = feature)
                ) +
                do.call(
                    ggplot2::geom_violin,
                    c(
                        list(
                            ggplot2::aes_string(fill = split_by),
                            scale = "width",
                            trim = TRUE,
                            show.legend = FALSE
                        ),
                        geom_violin_args
                    )
                ) +
                ggplot2::theme_classic() +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(
                        angle = 90, vjust = 0.5, hjust = 1
                    )
                )
        } else {
            p <- ggplot2::ggplot(
                data, ggplot2::aes_string(x = "1", y = feature)
            ) +
                do.call(
                    ggplot2::geom_violin,
                    c(
                        list(
                            scale = "width",
                            trim = TRUE,
                            fill = "#FF7F7F"
                        ),
                        geom_violin_args
                    )
                ) +
                ggplot2::theme_classic() +
                ggplot2::theme(
                    axis.title.x = ggplot2::element_blank(),
                    axis.text.x = ggplot2::element_blank(),
                    axis.ticks.x = ggplot2::element_blank()
                )
        }

        # Add jitter conditionally
        if (jitter) {
            p <- p + do.call(
                ggplot2::geom_jitter,
                c(
                    list(
                        height = 0,
                        size = pt_size,
                        alpha = alpha,
                        show.legend = FALSE
                    ),
                    geom_jitter_args
                )
            )
        }

        # Add common elements
        p <- p +
            ggplot2::labs(
                y = paste(feature, "Expression"),
                title = feature
            ) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(
                    hjust = 0.5, face = "bold", size = 16
                )
            )

        return(p)
    })

    # Combine plots into a grid
    final_plot <- patchwork::wrap_plots(plot_list, ncol = ncol)

    # Return the combined plot
    return(final_plot)
}

#' Plot Ridge Plots for Features
#'
#' Creates ridge plots for specified features from a `SpatialExperiment` object.
#'
#' @param object A `SpatialExperiment` object.
#' @param assay_name A character string specifying the assay to use for plotting.
#' @param features A character vector of feature names to plot.
#' @param split_by An optional column name in `colData` to split the ridges by.
#'   Default is `NULL`.
#' @param ncol Number of columns for arranging plots in the grid. Default is `1`.
#' @param geom_density_ridges_args A list of additional arguments passed to `geom_density_ridges`.
#' @return A `patchwork`ed `ggplot2` object.
#' @export
#' @importFrom ggplot2 ggplot aes_string labs theme_classic theme element_blank
#' @importFrom ggridges geom_density_ridges
#' @importFrom SummarizedExperiment colData assay
#' @importFrom dplyr filter
#' @importFrom patchwork wrap_plots
#' @examples
#' PlotRidge(object = spe, assay_name = "standardised", features = c("Gene1", "Gene2"))
PlotRidge <- function(object,
                      assay_name,
                      features,
                      split_by = NULL,
                      ncol = 1,
                      geom_density_ridges_args = list()) {

    # Extract metadata
    metadata <- SummarizedExperiment::colData(object) %>%
        as.data.frame()

    # Extract assay data and combine with metadata
    data <- SummarizedExperiment::assay(object, assay_name) %>%
        t() %>%
        as.data.frame() %>%
        cbind(metadata)

    # Create a list of ridge plots
    plot_list <- lapply(features, function(feature) {
        if (!is.null(split_by)) {
            p <- data %>%
                dplyr::filter(!is.na(.data[[split_by]])) %>%
                ggplot2::ggplot(
                    ggplot2::aes_string(
                        x = feature, y = split_by, fill = split_by
                    )
                ) +
                do.call(
                    ggridges::geom_density_ridges,
                    c(
                        list(show.legend = FALSE),
                        geom_density_ridges_args
                    )
                ) +
                ggplot2::theme_classic()
        } else {
            p <- ggplot2::ggplot(
                data, ggplot2::aes_string(x = feature, y = "1")
            ) +
                do.call(
                    ggridges::geom_density_ridges,
                    c(
                        list(fill = "#FF7F7F"),
                        geom_density_ridges_args
                    )
                ) +
                ggplot2::theme_classic() +
                ggplot2::theme(
                    axis.title.y = ggplot2::element_blank(),
                    axis.text.y = ggplot2::element_blank(),
                    axis.ticks.y = ggplot2::element_blank()
                )
        }

        # Add common elements
        p <- p +
            ggplot2::labs(x = paste(feature, "Expression"), title = feature) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(
                    hjust = 0.5, face = "bold", size = 16
                )
            )

        return(p)
    })

    # Combine plots into a grid
    final_plot <- patchwork::wrap_plots(plot_list, ncol = ncol)

    # Return the combined plot
    return(final_plot)
}

#' Plot Heatmap for Features
#'
#' Creates a heatmap of feature expression from a `SpatialExperiment` object.
#' Allows aggregation and annotation.
#'
#' @param object A `SpatialExperiment` object.
#' @param prefix A character string to filter clusters by prefix.
#' @param annot_by A column in `colData` used for annotation.
#' @param genes A character vector of gene names to include in the heatmap.
#'   Default is all row names.
#' @param assay_name A character string specifying the assay to use for plotting.
#'   Default is `"standardised"`.
#' @param disp_min Minimum value for display. Values below will be clipped. Default is `NULL`.
#' @param disp_max Maximum value for display. Values above will be clipped. Default is `NULL`.
#' @param seed A numeric seed for reproducibility. Default is `123`.
#' @param aggregate Logical, whether to aggregate data. Default is `FALSE`.
#' @param aggregate_by A column in `colData` to use for aggregation. Required if `aggregate = TRUE`.
#' @param ... Additional arguments passed to `dittoHeatmap`.
#' @return A `pheatmap` object.
#' @export
#' @importFrom SummarizedExperiment colData assay
#' @importFrom S4Vectors metadata
#' @importFrom scuttle aggregateAcrossCells
#' @importFrom dittoSeq dittoHeatmap dittoColors
#' @examples
#' PlotHeatmap(object = spe, prefix = "Cluster", annot_by = "DonorID", genes = c("Gene1", "Gene2"))
PlotHeatmap <- function(object,
                        prefix,
                        annot_by,
                        genes = rownames(object),
                        assay_name = "standardised",
                        disp_min = NULL,
                        disp_max = NULL,
                        seed = 123,
                        aggregate = FALSE,
                        aggregate_by = NULL,
                        ...) {
    set.seed(seed)

    # Ensure unique column names
    colnames(object) <- paste0(object$ImageID, "_", object$CellID)

    # Define a logical vector for valid rows
    valid_rows <-!is.na(SummarizedExperiment::colData(object)$Cluster) &
        grepl(
            paste0("^", prefix),
            SummarizedExperiment::colData(object)$Cluster
        )

    # Subset the object
    subset_object <- object[, valid_rows]

    # Sort clusters numerically
    cluster_levels <- levels(factor(subset_object$Cluster))
    sorted_levels <- cluster_levels[
        order(as.numeric(gsub(paste0(prefix, "_"), "", cluster_levels)))
    ]
    subset_object$Cluster <- factor(
        subset_object$Cluster, levels = sorted_levels
    )

    # Sample cells if necessary
    cell_count <- min(ncol(subset_object), 2000)
    cur_cells <- sample(seq_len(ncol(subset_object)), cell_count)
    subset_object <- subset_object[, cur_cells]

    # Clip the assay data
    mat <- SummarizedExperiment::assay(subset_object, assay_name)

    if (!is.null(disp_min)) {
        mat <- pmax(mat, disp_min)
    }
    if (!is.null(disp_max)) {
        mat <- pmin(mat, disp_max)
    }

    SummarizedExperiment::assay(subset_object, assay_name) <- mat

    # Set annotation colors
    colours <- list(
        "Cluster" = dittoSeq::dittoColors(1)[
            1:length(unique(subset_object$Cluster))
        ],
        "ImageID" = S4Vectors::metadata(object)$color_vectors$ImageID,
        "DonorID" = S4Vectors::metadata(object)$color_vectors$DonorID
    )

    # Aggregate if requested
    if (aggregate) {
        if (is.null(aggregate_by)) {
            stop("Please supply a variable to aggregate by.")
        }

        mean_aggregation <- scuttle::aggregateAcrossCells(
            x = as(subset_object, "SingleCellExperiment"),
            ids = SummarizedExperiment::colData(subset_object)$Cluster,
            statistics = "mean",
            use.assay.type = assay_name,
            subset.row = genes
        )

        p <- dittoSeq::dittoHeatmap(
            mean_aggregation,
            assay = assay_name,
            genes = genes,
            cluster_rows = FALSE,
            scale = "none",
            heatmap.colors = grDevices::colorRampPalette(
                c("#FF00FF", "black", "#FFFF00")
            )(100),
            annot.by = annot_by,
            annot.colors = colours[annot_by] %>% unname() %>% unlist(),
            ...
        )
        return(p)
    }

    # Generate heatmap without aggregation
    p <- dittoSeq::dittoHeatmap(
        subset_object,
        genes = genes,
        cluster_rows = FALSE,
        assay = assay_name,
        scale = "none",
        heatmap.colors = grDevices::colorRampPalette(
            c("#FF00FF", "black", "#FFFF00")
        )(100),
        annot.by = annot_by,
        annot.colors = colours[annot_by] %>% unname() %>% unlist(),
        ...
    )
    return(p)
}

#' Plot Dimensional Reduction Embeddings
#'
#' Generates scatter plots for dimensional reduction embeddings such as PCA or UMAP.
#'
#' @param object A `SpatialExperiment` object.
#' @param reduction A character string, either `"UMAP"` or `"PCA"`.
#' @param prefixes An optional vector of prefixes to filter clusters. Default is `NULL`.
#' @param assay_name A character string specifying the assay to use. Default is `"standardised"`.
#' @param var A column in `colData` to color points by. Default is `NULL`.
#' @param subsample An optional integer specifying the number of cells to subsample. Default is `NULL`.
#' @param dim A numeric vector specifying which dimensions to plot. Default is `c(1, 2)`.
#' @param seed A numeric seed for reproducibility. Default is `123`.
#' @param ... Additional arguments passed to `dittoDimPlot`.
#' @return A `ggplot2` object.
#' @export
#' @importFrom SummarizedExperiment colData assay
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom dittoSeq dittoDimPlot
#' @importFrom ggplot2 labs theme element_blank
#' @examples
#' PlotDim(object = spe, reduction = "UMAP", var = "Cluster")
PlotDim <- function(object,
                    reduction,
                    prefixes = NULL,
                    assay_name = "standardised",
                    var = NULL,
                    subsample = NULL,
                    dim = c(1, 2),
                    seed = 123,
                    ...) {
    set.seed(seed)

    # Ensure reduction type is valid
    reduction <- toupper(reduction)
    if (!(reduction %in% c("UMAP", "PCA"))) {
        stop("Please set 'reduction' to either 'UMAP' or 'PCA'.")
    }

    # Determine the appropriate reduction slot
    if (is.null(prefixes)) {
        slot <- paste0(reduction, "_full")
        if (!(slot %in% SingleCellExperiment::reducedDimNames(object))) {
            stop(
                paste0("Please run Perform", reduction," on the full data.")
            )
        }
    } else {
        slot <- paste0(reduction, "_subset")
        if (!(slot %in% SingleCellExperiment::reducedDimNames(object))) {
            stop(paste0(
                "Please run Perform",
                reduction,
                " with the specified prefixes."
            ))
        }

        # Validate prefixes for subset
        prefixes_match <- if (reduction == "UMAP") {
            setequal(
                S4Vectors::metadata(object)$UMAP_subset_prefixes, prefixes
            )
        } else {
            setequal(
                S4Vectors::metadata(object)$PCA_subset_prefixes, prefixes
            )
        }

        if (!prefixes_match) {
            stop(paste0("Please run Perform",
                        reduction,
                        " with the specified prefixes."))
        }
    }

    # Filter out rows with NA values in the reduction slot
    valid_cells <- !is.na(SingleCellExperiment::reducedDim(object, slot)[, 1])
    object <- object[, valid_cells]

    # Handle subsampling
    if (!is.null(subsample)) {
        if (subsample > ncol(object)) {
            stop("Subsample size exceeds number of available observations.")
        }
        subsample_indices <- sample(seq_len(ncol(object)), subsample)
        object <- object[, subsample_indices]
    }

    # Generate labels for axes
    x_label <- paste0(reduction, "_", dim[1])
    y_label <- paste0(reduction, "_", dim[2])

    # Generate the plot
    colnames(object) <- paste0("cell", seq_len(ncol(object)))
    plot <- dittoSeq::dittoDimPlot(
        object,
        var = var,
        reduction.use = slot,
        legend.title = var,
        assay = assay_name,
        dim.1 = dim[1],
        dim.2 = dim[2],
        ...
    ) +
        ggplot2::labs(x = x_label, y = y_label) +
        ggplot2::theme(
            plot.title = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
        )

    return(plot)
}

#' Plot Scatter of Two Features
#'
#' Creates a scatter plot of two features from a `SpatialExperiment` object.
#'
#' @param object A `SpatialExperiment` object.
#' @param feature1 A character string specifying the first feature (x-axis).
#' @param feature2 A character string specifying the second feature (y-axis).
#' @param assay_name1 A character string specifying the assay for `feature1`. Default is `"standardised"`.
#' @param assay_name2 A character string specifying the assay for `feature2`. Default is `"standardised"`.
#' @param colour_by A column in `colData` to color points by. Default is `NULL`.
#' @param remove_na Logical, whether to remove cells with `NA` in `colour_by`. Default is `TRUE`.
#' @return A `ggplot2` object.
#' @export
#' @importFrom SummarizedExperiment colData assay
#' @importFrom dittoSeq dittoScatterPlot
#' @importFrom ggplot2 theme element_blank
#' @examples
#' PlotScatter(object = spe, feature1 = "Gene1", feature2 = "Gene2", colour_by = "Cluster")
PlotScatter <- function(object,
                        feature1,
                        feature2,
                        assay_name1 = "standardised",
                        assay_name2 = "standardised",
                        colour_by = NULL,
                        remove_na = TRUE) {

    # Remove rows with NA values in the specified column, if requested
    if (remove_na && !is.null(colour_by) && (
        colour_by %in% names(SummarizedExperiment::colData(object))
    )
    ) {
        valid_cells <- !is.na(object[[colour_by]])
        object <- object[, valid_cells]
    }

    # Ensure unique column names
    colnames(object) <- paste0("cell", seq_len(ncol(object)))

    # Generate scatter plot
    plot <- dittoSeq::dittoScatterPlot(
        object,
        x.var = feature1,
        y.var = feature2,
        assay.x = assay_name1,
        assay.y = assay_name2,
        color.var = colour_by
    ) +
        ggplot2::theme(plot.title = ggplot2::element_blank())

    return(plot)
}
