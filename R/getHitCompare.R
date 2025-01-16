#' get paired as events / exons
#'
#' @param data_df input with control and test noted
#' @param output_location output location for plots
#' @return figure showing proximal / distal shift and table of results
#' @importFrom dplyr %>% full_join ungroup mutate group_by summarise filter arrange desc left_join
#' @importFrom pheatmap pheatmap
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom data.table fread
#' @importFrom purrr reduce
#' @importFrom grid grid.newpage grid.draw
#' @importFrom circlize colorRamp2
#' @export
#'
getHitCompare <- function(data_df, output_location) {
  sample_types <- list()
  for (i in 1:nrow(data_df)) {
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
        select = c("gene", "exon", "HITindex")
      )

      # Rename HITindex column to include the 2nd element and its position in u
      setnames(data, "HITindex", paste0("HITindex_", i[2], "_", j))
      return(data)
    }) %>%
    purrr::reduce(function(df1, df2) {
      dplyr::full_join(df1, df2, by = c("gene", "exon"))
    })


  merged_data_df <- as.data.frame(merged_data)
  NAs_df <- is.na(merged_data_df)

  merged_data_df$NAcount_control <- rowMeans(NAs_df[,grep("control", colnames(NAs_df))])
  merged_data_df$NAcount_test <- rowMeans(NAs_df[,grep("test", colnames(NAs_df))])
  threshold <- 0.25
  filtered_data <- merged_data_df[merged_data_df$NAcount_control < threshold & merged_data_df$NAcount_test < threshold,]


  heatmap_data <- filtered_data %>%
    tidyr::pivot_longer(
      cols = contains(c("control", "test")),
      names_to = "condition",
      values_to = "HITindex"
    ) %>%
    dplyr::group_by(gene, exon) %>%
    dplyr::mutate(row_id = cur_group_id()) %>%
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
    column_to_rownames(var = "gene_exon") %>%
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

  # Generate the mean heatmap with custom colors
  totalHeatmap <- pheatmap::pheatmap(
    heatmap_matrix,
    cluster_rows = FALSE,  # Disable row clustering
    cluster_cols = FALSE,  # Disable column clustering
    color = custom_colors_total(seq(-1, 1, length.out = 100)),  # Gradient between control and test colors
    show_rownames = F,  # Show row names (gene_exon)
    show_colnames = F,  # Show column names (control and test)
    annotation_col = total_col_annotation,  # Add labels for control and test
    annotation_colors = list(group = control_test_colors)  # Apply custom colors
  )

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
    column_to_rownames(var = "gene_exon") %>%
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
    show_rownames = F,  # Show row names (gene_exon)
    show_colnames = F,  # Show column names (control and test)
    annotation_col = mean_col_annotation,  # Add labels for control and test
    annotation_colors = list(group = control_test_colors)  # Apply custom colors
  )


  pdf(paste0(output_location, "Foreground/", "HITheatmap.pdf"))
  grid::grid.newpage()  # Start a new page
  grid::grid.draw(totalHeatmap$gtable)  # Draw the pheatmap plot

  # Print the mean heatmap on the second page
  grid::grid.newpage()  # Start a new page
  grid::grid.draw(meanHeatmap$gtable)
  dev.off()

  return(list(totalHeatmap = totalHeatmap,
              meanHeatmap = meanHeatmap))
}
