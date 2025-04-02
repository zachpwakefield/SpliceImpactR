#' identify any global shift in proximal or distal exon use
#'
#' @param type event type
#' @param exon_pairs recent dataframe made from getPaired
#' @param ep_supp recent dataframe made from getPaired
#' @param output_location location for saving
#' @return figure showing proximal / distal shift and table of results
#' @importFrom dplyr %>% mutate select distinct left_join case_when
#' @importFrom ggplot2 ggplot aes geom_segment geom_text scale_color_manual coord_cartesian element_blank theme labs theme_minimal element_text
#' @importFrom grid arrow
#' @importFrom tidyr separate
#' @importFrom grid unit
#'
#' @examples
#' pdir <- system.file("extdata", package="SpliceImpactR")
#' dataDirectory <- paste0(pdir, "/")
#' test_group <- paste0(dataDirectory, "rawData/", c("test1","test2", "test3"))
#' control_group <- paste0(dataDirectory, "rawData/", c("control1", "control2", "control3"))
#' transDF <- readr::read_csv(paste0(dataDirectory, "transcripts_limited_transDF.csv"))
#' c_trans <- readr::read_lines(paste0(dataDirectory, "transcripts_limited_c_trans.csv"))
#'
#' transcripts_sample <- list(transDF = transDF,
#'                            c_trans = c_trans)
#'
#' gtf_sample <- list(gtf = readr::read_csv(paste0(dataDirectory, "gtf_limited.csv")),
#'             transcript_gtf = readr::read_csv(paste0(dataDirectory, "transcript_gtf_limited.csv")))
#' translations_sample <- readr::read_lines(paste0(dataDirectory, "translations_limited.csv"))
#' biomart_data_sample <- readr::read_csv(paste0(dataDirectory, "biomart_data_sample.csv"))
#'
#'
#' result <- differential_inclusion_HITindex(test_names = test_group,
#'                                           control_names = control_group,
#'                                           et = "AFE",
#'                                           outlier_threshold = "Inf",
#'                                           minReads = 10,
#'                                           min_prop_samples = 0,
#'                                           chosen_method = "qbGLM"
#'                                           )
#'
#' fg <- getForeground(input = result,
#'                             test_names = test_group,
#'                             control_names = control_group,
#'                             thresh = .1,
#'                             fdr = .05,
#'                             mOverlap = .1,
#'                             exon_type = "AFE",
#'                             output_location = NULL,
#'                             cores = 1,
#'                             gtf = gtf_sample,
#'                             max_zero_prop = 1,
#'                             min_prop_samples = 0,
#'                             translations = translations_sample)
#' library(msa)
#' pfg <- getPaired(foreground = fg$proBed,
#'           et = "AFE",
#'           nucleotides = transcripts_sample,
#'           newGTF = gtf_sample,
#'           cores = 1,
#'           output_location = NULL,
#'           saveAlignments = FALSE,
#'           exon_data = biomart_data_sample)
#'
#' proxShift <- getProximalShift("AFE", pfg$exon_pairs, pfg$paired_proBed, output_location = NULL)
#' @export
#'
getProximalShift <- function(type, exon_pairs, ep_supp, output_location = NULL) {
  ep_supp <- ep_supp %>% dplyr::select(gene, strand) %>% dplyr::distinct()
  exon_pairs <- dplyr::left_join(exon_pairs, ep_supp) %>%
    tidyr::separate(pos_exon_id, into = c('pos_chr', 'pos_start', 'pos_stop'), sep = "[:\\-]") %>%
    tidyr::separate(neg_exon_id, into = c('neg_chr', 'neg_start', 'neg_stop'), sep = "[:\\-]") %>%
    dplyr::mutate(
      neg_start = as.numeric(neg_start),
      pos_start = as.numeric(pos_start),
      type = type,
      label = dplyr::case_when(
        # Conditions for "distal"
        type == "AFE" & strand == "+" & pos_start < neg_start ~ "distal",
        type == "AFE" & strand == "-" & pos_start > neg_start ~ "distal",
        type == "ALE" & strand == "+" & pos_start > neg_start ~ "distal",
        type == "ALE" & strand == "-" & pos_start < neg_start ~ "distal",
        # Default to "proximal"
        TRUE ~ "proximal"
      )
    )

  data <- data.frame(table(exon_pairs$label))
  colnames(data) <- c('group', 'size')
  data$type <- type

  if (max(data$size) > 30) {
    data$size_adjusted <- c(data$size)*(30/max(data$size))
  } else {
    data$size_adjusted <- data$size
  }

  shift <- 10
  # Create start and end points for arrows
  data2 <- data %>%
    dplyr::mutate(
      # Reverse arrow direction for ALE
      x_start = ifelse(type == "ALE",
                       ifelse(group == "proximal", -size_adjusted / 2 + shift, size_adjusted / 2 + shift),
                       ifelse(group == "proximal", size_adjusted / 2 + shift, -size_adjusted / 2 + shift)),
      x_end = ifelse(type == "ALE",
                     ifelse(group == "proximal", size_adjusted / 2 + shift, -size_adjusted / 2 + shift),
                     ifelse(group == "proximal", -size_adjusted / 2 + shift, size_adjusted / 2 + shift)),
      y = ifelse(group == "proximal", 1, 1.1)  # Adjust y positions
    )
  x_limits <- c(min(data2$x_end) - 15, max(data2$x_start) + 15)  # Add padding to fit the full arrow

  proxPlot <- ggplot2::ggplot(data2) +
    # Add arrows
    ggplot2::geom_segment(
      ggplot2::aes(x = x_start, xend = x_end, y = y, yend = y, color = group),
      arrow = grid::arrow(length = unit(0.2, "cm")),
      size = 1
    ) +
    # Add labels to the side of arrows
    ggplot2::geom_text(
      ggplot2::aes(
        x = ifelse(group == "proximal", min(x_end, x_start) - 10, min(x_end, x_start) - 10),
        y = y,
        label = paste0(group, " (", size, ")"),  # Append counts to labels
        color = group
      ),
      size = 5, hjust = 0.5
    ) +
    # Customize colors
    ggplot2::scale_color_manual(values = c("distal" = "deeppink4", "proximal" = "cadetblue4")) +
    # Fix y-axis to avoid relative positioning issues
    ggplot2::coord_cartesian(xlim = x_limits, ylim = c(0.7, 1.2)) +  # Set y-axis limits to fix alignment
    # Customize plot appearance
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0("Proximal / Distal shift for ", type)) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin=grid::unit(c(10,10,10,10), "pt"),
      plot.title = ggplot2::element_text(hjust = 0.5, vjust = -4),
      legend.position = "none"  # Remove legend
    )

  if (!is.null(output_location)) {
    pdf(paste0(output_location, "pairedOutput/", "proximalShiftPlot.pdf"), height = 2)
    print(proxPlot)
    dev.off()
  }


  return(list(proxPlot = proxPlot,
              proximalShift = table(exon_pairs$label)))
}

