#' Perform analysis on how the HIT Index value is changing on an exon-basis
#'
#' @param data_df input with control and test noted
#' @param output_location output location for plots
#' @param threshold for NA in HIT Index
#' @return figures about HIT Index and csv for differentially used exons (diff exon types)
#' @importFrom dplyr %>% full_join ungroup mutate group_by summarise filter arrange desc left_join cur_group_id
#' @importFrom pheatmap pheatmap
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom data.table fread setnames
#' @importFrom purrr reduce
#' @importFrom grid grid.newpage grid.draw
#' @importFrom tibble column_to_rownames
#'
#' @examples
#'
#' pdir <- system.file("extdata", package="SpliceImpactR")
#' dataDirectory <- paste0(pdir, "/rawData/")
#' test_group <- paste0(dataDirectory, c("test1", "test2", "test3"))
#' control_group <- paste0(dataDirectory, c("control1", "control2", "control3"))
#' data_df <- data.frame(
#'     sample_names = c(control_group, test_group),
#'     phenotype_names = c(
#'       rep("control", length(control_group)),
#'       rep("test", length(test_group))
#'      ),
#'    stringsAsFactors = FALSE
#'   )
#'   data_df$utc <- "control"
#'   data_df$utc[data_df$phenotype_names == unique(data_df$phenotype_names)[2]] <- "test"
#' compareHIT <- getHitCompare(data_df, output_location = NULL, threshold = .1)
#'
#' @export
#'
getHitCompare <- function(data_df, output_location, threshold = .1) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Please install circlize to use getHitCompare.")
  }
  sample_types <- list()
  for (i in seq_len(nrow(data_df))) {
    sample_types <- c(sample_types, list(c(data_df$sample_names[i], data_df$utc[i], data_df$phenotype_names[i])))
  }
  sample_list <- c(sample_types[which(unlist(lapply(sample_types, "[[", 2)) == "control")],
                   sample_types[which(unlist(lapply(sample_types, "[[", 2)) == "test")])
  u <- unlist(lapply(sample_list, "[[", 1))  # Extract the first element from each entry in sample_list

  merged_data <- sample_list %>%
    lapply(function(i) {
      # Find the position of i[1] in u
      j <- which(u == i[1])

      # Read the file
      data <- data.table::fread(
        paste0(i[1], ".exon"),
        select = c("gene", "exon", "nUP", "nDOWN", "HITindex")
      )

      data <- data[data$nUP > 10 | data$nDOWN > 10,c('gene', 'exon', 'HITindex')]

      # Rename HITindex column to include the 2nd element and its position in u
      data.table::setnames(data, "HITindex", paste0("HITindex_", i[2], "_", j))
      return(data)
    }) %>%
    purrr::reduce(function(df1, df2) {
      dplyr::full_join(df1, df2, by = c("gene", "exon"))
    })


  diHIT <- getDifferentialHIT(merged_data)
  diPlot <- diHIT_plots(diHIT)

  merged_data_df <- as.data.frame(merged_data)
  NAs_df <- is.na(merged_data_df)

  merged_data_df$NAcount_control <- rowMeans(NAs_df[,grep("control", colnames(NAs_df))])
  merged_data_df$NAcount_test <- rowMeans(NAs_df[,grep("test", colnames(NAs_df))])
  filtered_data <- merged_data_df[merged_data_df$NAcount_control <= threshold & merged_data_df$NAcount_test <= threshold,]


  heatmap_data <- filtered_data %>%
    tidyr::pivot_longer(
      cols = contains(c("control", "test")),
      names_to = "condition",
      values_to = "HITindex"
    ) %>%
    dplyr::group_by(gene, exon) %>%
    dplyr::mutate(row_id = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(condition = factor(condition, levels = c(
      grep("control", unique(condition), value = TRUE),
      grep("test", unique(condition), value = TRUE)
    )))


  heatmap_data <- heatmap_data %>%
    dplyr::mutate(
      status = dplyr::case_when(
        grepl("control", condition) ~ "control",
        grepl("test", condition) ~ "test",
        TRUE ~ NA_character_
      )
    )

  # Replace NAs with group mean
  heatmap_data_filled <- heatmap_data %>%
    dplyr::group_by(gene, exon, status) %>%
    dplyr::mutate(
      HITindex = ifelse(
        is.na(HITindex),
        mean(HITindex, na.rm = TRUE),  # Use group mean for NA values
        HITindex
      )
    ) %>%
    dplyr::ungroup()

  heatmap_data_filled <- as.data.frame(heatmap_data_filled)

  # Create a unique identifier for rows (gene_exon)
  heatmap_data_filled <- heatmap_data_filled %>%
    dplyr::mutate(gene_exon = paste(gene, exon, sep = "_"))

  # Compute the mean HITindex for control columns for sorting
  control_means <- heatmap_data_filled %>%
    dplyr::filter(grepl("control", condition)) %>%  # Select only control rows
    dplyr::group_by(gene_exon) %>%
    dplyr::summarise(mean_control = mean(HITindex, na.rm = TRUE), .groups = "drop")

  # Join the control means back to the heatmap data and sort by mean control
  heatmap_data_sorted <- heatmap_data_filled %>%
    dplyr::left_join(control_means, by = "gene_exon") %>%
    dplyr::arrange(dplyr::desc(mean_control))  # Sort rows by mean control values

  # Pivot the data to create a wide matrix-like structure
  heatmap_matrix <- heatmap_data_sorted %>%
    dplyr::select(gene_exon, condition, HITindex) %>%
    tidyr::pivot_wider(names_from = condition, values_from = HITindex) %>%
    tibble::column_to_rownames(var = "gene_exon") %>%
    as.matrix()

  # Remove columns with NA in their names (if applicable)
  heatmap_matrix <- heatmap_matrix[, !(grepl("NA", colnames(heatmap_matrix)))]

  # Define custom colors for the heatmap
  custom_colors_total <- circlize::colorRamp2(
    breaks = c(1, 0.8, 0.3, -0.3, -0.8, -1),  # Interval breaks
    colors = c("red", "turquoise", "grey", "grey", "orange", "purple")  # Corresponding colors
  )

  total_col_annotation <- data.frame(
    group = ifelse(grepl("control", colnames(heatmap_matrix)), "control", "test")
  )
  rownames(total_col_annotation) <- colnames(heatmap_matrix)
  control_test_colors <- c("control" = "deeppink4", "test" = "cadetblue4")

  if (length(data_df$sample_names[data_df$utc == 'test']) <= 10 &
      length(data_df$sample_names[data_df$utc == 'control']) <= 10) {
    # Generate the mean heatmap with custom colors
    totalHeatmap <- pheatmap::pheatmap(
      heatmap_matrix,
      cluster_rows = FALSE,  # Disable row clustering
      cluster_cols = FALSE,  # Disable column clustering
      color = custom_colors_total(seq(-1, 1, length.out = 100)),  # Gradient between control and test colors
      show_rownames = FALSE,  # Show row names (gene_exon)
      show_colnames = FALSE,  # Show column names (control and test)
      annotation_col = total_col_annotation,  # Add labels for control and test
      annotation_colors = list(group = control_test_colors), silent = TRUE # Apply custom colors

    )
  }
  heatmap_mean_data <- heatmap_data_filled %>%
    dplyr::mutate(gene_exon = paste(gene, exon, sep = "_"))

  # Compute the mean HITindex for control and test columns
  mean_heatmap_data <- heatmap_mean_data %>%
    dplyr::mutate(status = ifelse(grepl("control", condition), "control", "test")) %>%
    dplyr::group_by(gene_exon, status) %>%
    dplyr::summarise(mean_HITindex = mean(HITindex, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = status, values_from = mean_HITindex)

  # Sort rows by the mean of the control group
  mean_heatmap_data_sorted <- mean_heatmap_data %>%
    dplyr::arrange(dplyr::desc(control))  # Sort by control mean in descending order

  # Convert to matrix format for heatmap
  mean_heatmap_matrix <- mean_heatmap_data_sorted %>%
    tibble::column_to_rownames(var = "gene_exon") %>%
    as.matrix()

  # Define custom colors for the mean heatmap
  custom_colors_mean <- circlize::colorRamp2(
    breaks = c(1, 0.8, 0.3, -0.3, -0.8, -1),  # Interval breaks
    colors = c("red", "turquoise", "grey", "grey", "orange", "purple")  # Corresponding colors
  )


  control_test_colors <- c("control" = "deeppink4", "test" = "cadetblue4")
  mean_col_annotation <- data.frame(
    group = ifelse(grepl("control", colnames(mean_heatmap_matrix)), "control", "test")
  )
  rownames(mean_col_annotation) <- colnames(mean_heatmap_matrix)

  # Generate the mean heatmap with custom colors
  meanHeatmap <- pheatmap::pheatmap(
    mean_heatmap_matrix,
    cluster_rows = FALSE,  # Disable row clustering
    cluster_cols = FALSE,  # Disable column clustering
    color = custom_colors_mean(seq(-1, 1, length.out = 100)),  # Gradient between control and test colors
    show_rownames = FALSE,  # Show row names (gene_exon)
    show_colnames = FALSE,  # Show column names (control and test)
    annotation_col = mean_col_annotation,  # Add labels for control and test
    annotation_colors = list(group = control_test_colors), silent = TRUE  # Apply custom colors
  )
  if (!is.null(output_location)) {
    if (length(data_df$sample_names[data_df$utc == 'test']) <= 10 &
        length(data_df$sample_names[data_df$utc == 'control']) <= 10) {
      pdf(paste0(output_location, "HITheatmap.pdf"))
      grid::grid.newpage()  # Start a new page
      grid::grid.draw(totalHeatmap$gtable)  # Draw the pheatmap plot

      # Print the mean heatmap on the second page
      grid::grid.newpage()  # Start a new page
      grid::grid.draw(meanHeatmap$gtable)
      dev.off()
    } else {
      pdf(paste0(output_location, "HITheatmap.pdf"))
      grid::grid.newpage()  # Start a new page
      grid::grid.draw(meanHeatmap$gtable)
      dev.off()
    }

    pdf(paste0(output_location, "dotPlotHIT.pdf"))
    print(diPlot$dotPlot)
    dev.off()
    pdf(paste0(output_location, "volcanoPlotHIT.pdf"))
    print(diPlot$volcanoPlot)
    dev.off()

    write_csv(diHIT, paste0(output_location, "deltaHIT.csv"))
  }


  return(list(diPlot = diPlot,
              meanHeatmap = meanHeatmap,
              diHIT = diHIT))
}

#' get differential HIT
#'
#' @param merged_data data from internal working of getHITCompare
#' @return final_results a dataframe for HIT Index summary
#' @importFrom stats wilcox.test p.adjust
#' @importFrom dplyr left_join mutate select summarise group_by
#' @importFrom tidyr pivot_wider pivot_longer
#' @keywords internal
#'
#'
getDifferentialHIT <- function(merged_data) {
  long_data <- merged_data %>%
    tidyr::pivot_longer(
      cols = starts_with("HITindex"),
      names_to = "condition_sample",
      values_to = "HITindex"
    ) %>%
    dplyr::mutate(
      condition = ifelse(grepl("_control_", condition_sample), "control", "test")
    )

  avg_hitindex <- long_data %>%
    dplyr::group_by(gene, exon, condition) %>%
    dplyr::summarise(
      avg_HITindex = mean(HITindex, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = condition,
      values_from = avg_HITindex,
      names_prefix = "",
      values_fill = NA
    )

  avg_hitindex <- avg_hitindex %>%
    dplyr::mutate(delta_HIT = abs(test - control))

  wilcox_results <- long_data %>%
    dplyr::group_by(gene, exon) %>%
    dplyr::summarise(
      p_value = if (sum(condition == "control" & !is.na(HITindex)) >= 2 &&
                    sum(condition == "test" & !is.na(HITindex)) >= 2) {
        stats::wilcox.test(
          HITindex[condition == "control"],
          HITindex[condition == "test"],
          exact = FALSE
        )$p.value
      } else {
        NA
      },
      .groups = "drop"
    )

  wilcox_results <- wilcox_results %>%
    dplyr::mutate(fdr = stats::p.adjust(p_value, method = "fdr"))

  final_results <- avg_hitindex %>%
    dplyr::left_join(wilcox_results, by = c("gene", "exon"))

  final_results <- final_results %>%
    dplyr::select(gene, exon, control_HIT = control, test_HIT = test, delta_HIT, p_value, fdr)


  final_results <- final_results[!is.na(final_results$fdr),]
  return(final_results)
}

#' get HITindex plots
#'
#' @param final_results from getDifferentialHIT
#' @return figures showing hit index plots (dot plot, volcano)
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_minimal theme element_text geom_hline geom_vline scale_color_gradient2 labs geom_abline
#' @keywords internal
#'
diHIT_plots <- function(final_results) {
  volcanoPlot <- ggplot2::ggplot(final_results, ggplot2::aes(x = delta_HIT, y = -log10(fdr))) +
    ggplot2::geom_point(ggplot2::aes(color = fdr < 0.05 & delta_HIT > .5), alpha = 0.7, size = 2) +  # Highlight significant points
    ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +  # Red for significant, gray for others
    ggplot2::labs(
      x = expression(Delta ~ HIT ~ "Index"),
      y = expression(-log[10] ~ "(FDR)"),
      title = "Volcano Plot of HIT Changes"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = ggplot2::element_text(size = 12),
      legend.position = "none"
    ) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add threshold line for significance
    ggplot2::geom_vline(xintercept = .5, linetype = "dotted", color = "black")  # Add zero effect line

  dotPlot <- ggplot2::ggplot(final_results, ggplot2::aes(x = control_HIT, y = test_HIT)) +
    ggplot2::geom_point(ggplot2::aes(color = delta_HIT), alpha = 0.7, size = 2) +  # Color by delta_HIT for gradient effect
    ggplot2::scale_color_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      name = expression(Delta ~ HIT ~ "Index")
    ) +
    ggplot2::labs(
      x = "Average HIT in Control",
      y = "Average HIT in Test",
      title = "Dot Plot: Control vs Test HIT Averages"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 10),
      legend.position = "right"
    ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  # Diagonal reference line

  return(list(dotPlot = dotPlot,
              volcanoPlot = volcanoPlot))
}
