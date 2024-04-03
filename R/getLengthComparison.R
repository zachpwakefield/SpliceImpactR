#' Get comparison between length across isoform swaps due to AS
#'
#' @param paired_df output df from getPaired
#' @param output_location location where everything is being saved
#' @return 1 combined plot.
#' @examples
#' getOverviewComparison(paired_combined_rows, output_path)
getLengthComparison <- function(paired_df, output_location) {

  paired_df$protLength <- nchar(paired_df$prot)
  paired_df$protLength[paired_df$prot == "none"] <- 0
  pc <- unlist(lapply(seq(1, nrow(paired_df), 2), function(x) {
    if (paired_df$protLength[x+1] != 0 & paired_df$protLength[x] != 0 ) {
      c(x, x+1)
    }
  }
  ))

  proteinLength <- data.frame(protLength = paired_df$protLength[pc],
                              type = rep(c("test", "control"), length(pc)/2))

  posneg_change_in_length <- unlist(lapply(seq(1, nrow(paired_df), by = 2), function(x) {
    if (paired_df$protLength[x] != 0 & paired_df$protLength[x+1] != 0) {
      paired_df$protLength[x+1]-paired_df$protLength[x]
    }}))

  proteinLength$deltaLength <- posneg_change_in_length
  whichPC <- unlist(lapply(seq(1, nrow(paired_df), by = 2), function(x) {
    if (paired_df$protLength[x] == 0 & paired_df$protLength[x+1] == 0) {
      "both"
    } else if (paired_df$protLength[x] == 0) {
      "test"
    } else if (paired_df$protLength[x+1] == 0) {
      "control"
    } else {
      "pc"
    }
  }))

  dfPC <- data.frame(count = as.numeric(table(whichPC)),
                     type = names(table(whichPC)))

  lengthPlot <- ggplot2::ggplot(proteinLength, ggplot2::aes(x = .data$type, y = .data$protLength, fill = .data$type)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values=c("brown", "chartreuse4"), breaks = c("control", "test")) + ggplot2::xlab("Group") +
    ggplot2::ylab("Protein Length") +
    ggplot2::theme(legend.position = "none")

  changeDistribution_temp <- ggplot2::ggplot(proteinLength, ggplot2::aes(x = .data$deltaLength)) +
    ggplot2::geom_density(fill = "deeppink4", alpha = .4) +
    ggplot2::theme_bw()

  maxVal <- round(max(ggplot2::ggplot_build(changeDistribution_temp)$data[[1]]$y), 3)

  changeDistribution <- changeDistribution_temp + ggplot2::scale_y_continuous(breaks=seq(0,maxVal, maxVal)) + ggplot2::coord_flip() + ggplot2::ylab("Density") +
    ggplot2::xlab("Change in Protein Length")


  proteinCodingPlot <- ggplot2::ggplot(dfPC, ggplot2::aes(x = .data$count, y = .data$type, fill = .data$type)) + geom_bar(stat="identity") + theme_bw() +
    ggplot2::scale_fill_manual(values=c("brown", "deeppink4", "chartreuse4", "azure4"), breaks = c("pc", "control", "test", 'both'))

  comb_plot <- ggpubr::ggarrange(lengthPlot, changeDistribution, proteinCodingPlot, nrow = 1, widths = c(3, 2.5, 3))

  pdf(paste0(output_location, "pairedOutput/", "paired_length_comparison_plots.pdf"), height = 8, width = 12)
  print(comb_plot)
  dev.off()

  return(comb_plot)
}
