#' Get comparison between length across isoform swaps due to AS
#'
#' @param data_df df with paths, phenotype names, class
#' @param paired_df from getPaired df
#' @param output_location location where everything is being saved
#'
#' @return 1 combined plot.
#'
#' @description
#' Uses the paired results to identify any potential changes in lengths across phenotype
#'
#' @importFrom ggpubr ggpaired ggarrange stat_compare_means
#' @importFrom ggplot2 xlab ylab aes ggplot element_blank geom_density element_line theme geom_density scale_y_continuous scale_fill_manual coord_flip ggplot_build geom_bar theme_bw
#' @export
getLengthComparison <- function(data_df, paired_df, output_location) {

  paired_df$protLength <- nchar(paired_df$prot)
  paired_df$protLength[paired_df$prot == "none"] <- 0
  pc <- unlist(lapply(seq(1, nrow(paired_df), 2), function(x) {
    if (paired_df$protLength[x+1] != 0 & paired_df$protLength[x] != 0 ) {
      c(x, x+1)
    }
  }
  ))



  proteinLength <- data.frame(protLength = paired_df$protLength[pc],
                              type = rep(c(unique(data_df$phenotype_names[data_df$utc == "test"]),
                                           unique(data_df$phenotype_names[data_df$utc == "control"])), length(pc)/2))

  posneg_change_in_length <- unlist(lapply(seq(1, nrow(paired_df), by = 2), function(x) {
    if (paired_df$protLength[x] != 0 & paired_df$protLength[x+1] != 0) {
      paired_df$protLength[x]-paired_df$protLength[x+1]
    }}))

  proteinLength$deltaLength <- posneg_change_in_length
  whichPC <- unlist(lapply(seq(1, nrow(paired_df), by = 2), function(x) {
    if (paired_df$protLength[x] == 0 & paired_df$protLength[x+1] == 0) {
      "noPC"
    } else if (paired_df$protLength[x] == 0) {
      unique(data_df$phenotype_names[data_df$utc == "test"])
    } else if (paired_df$protLength[x+1] == 0) {
      unique(data_df$phenotype_names[data_df$utc == "control"])
    } else {
      "bothPC"
    }
  }))
  proteinLengthPaired <- data.frame(test = as.integer(paired_df$protLength[pc][seq(1, length(pc), by = 2)]),
                                    cont = as.integer(paired_df$protLength[pc][seq(2, length(pc), by = 2)]))
  names(proteinLengthPaired) <- c(unique(data_df$phenotype_names[data_df$utc == "test"]), unique(data_df$phenotype_names[data_df$utc == "control"]))
  pairedLengthPlot <- ggpubr::ggpaired(proteinLengthPaired,
                                       cond1 = unique(data_df$phenotype_names[data_df$utc == "control"]),
                                       cond2 = unique(data_df$phenotype_names[data_df$utc == "test"]),
                                       line.color = "grey", line.size = ifelse(nrow(proteinLengthPaired) > 20, 0, 0.4), point.size = ifelse(nrow(proteinLengthPaired) > 20, 0, 1.2),
                                       fill = "condition")+
    ggpubr::stat_compare_means(paired = TRUE) +
    ggplot2::scale_fill_manual(values=c("brown", "chartreuse4")) +
    ggplot2::theme(legend.position = "none")

  dfPC <- data.frame(count = as.numeric(table(whichPC)),
                     type = names(table(whichPC)))

  changeDistribution_temp <- ggplot2::ggplot(proteinLength, ggplot2::aes(x = .data$deltaLength)) +
    ggplot2::geom_density(fill = "deeppink4", alpha = .4) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black"))

  maxVal <- round(max(ggplot2::ggplot_build(changeDistribution_temp)$data[[1]]$y), 3)

  changeDistribution <- changeDistribution_temp + ggplot2::scale_y_continuous(breaks=seq(0,maxVal, maxVal)) + ggplot2::coord_flip() + ggplot2::ylab("Density") +
    ggplot2::xlab("Change in Protein Length") +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black"))


  proteinCodingPlot <- ggplot2::ggplot(dfPC, ggplot2::aes(x = .data$count, y = .data$type, fill = .data$type)) + ggplot2::geom_bar(stat="identity") + ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values=c("brown", "deeppink4", "chartreuse4", "azure4"), breaks = c("bothPC",
                                                                                                   unique(data_df$phenotype_names[data_df$utc == "control"]),
                                                                                                   unique(data_df$phenotype_names[data_df$utc == "test"]),
                                                                                                   'noPC'))

  comb_plot <- ggpubr::ggarrange(pairedLengthPlot, changeDistribution, proteinCodingPlot, nrow = 1, widths = c(3, 2.5, 3))

  pdf(paste0(output_location, "pairedOutput/", "paired_length_comparison_plots.pdf"), height = 8, width = 12)
  print(comb_plot)
  dev.off()

  return(comb_plot)
}
