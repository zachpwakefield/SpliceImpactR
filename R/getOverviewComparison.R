#' Get comparison between phenotypes: count of AS event, count of AS event per gene, ECDF compare
#'
#' @param data_list A list containing dataframes of each output file from HIT Index or rMATS as index 1 and "control" or "test" as index 2 for each position in the list
#' @param sample_list sample_types_sorted from di function
#' @param exon_type type of AS being investigated
#' @return 3 individual plots and 1 combined plot.
#' @examples
#' getOverviewComparison(load_output, sample_types_sorted, "AFE")
getOverviewComparison <- function(data_list, sample_list, exon_type) {

  ## Number of AS across phenotype
  nASE <- lapply(data_list, function(x) nrow(x[x[,grep("PSI", colnames(x))] > 0,]))
  dfCount <- data.frame(AScount = unlist(nASE),
                        type = unlist(lapply(sample_list, "[[", 2)))
  if (length(sample_list) > 8) {
    p1 <- ggplot2::ggplot(dfCount, ggplot2::aes(x = type, y = AScount, fill = type)) +
      ggplot2::geom_dotplot(binaxis='y', stackdir='center') + ggplot2::theme_bw()+
      ggplot2::scale_fill_manual(values=c("brown", "chartreuse4"))+ ggplot2::xlab("Group") +
      ggplot2::ylab("Number of Events") +ggplot2::geom_violin(fill = NA)
  } else {
    p1 <- ggplot2::ggplot(dfCount, ggplot2::aes(x = type, y = AScount, fill = type)) +
      ggplot2::geom_dotplot(binaxis='y', stackdir='center') + ggplot2::theme_bw()+
      ggplot2::scale_fill_manual(values=c("brown", "chartreuse4"))+ ggplot2::xlab("Group") +
      ggplot2::ylab("Number of Events")
  }


  ## Number of AS per gene across phenotype
  nASpg <- lapply(data_list, function(x) mean(as.numeric(table(x$gene ))))
  dfASpg <- data.frame(ASpg = unlist(nASpg),
                       type = unlist(lapply(sample_list, "[[", 2)))
  if (length(sample_list) > 8) {
    p2 <- ggplot2::ggplot(dfASpg, ggplot2::aes(x = type, y = ASpg, fill = type)) + ggplot2::geom_dotplot(binaxis='y', stackdir='center') + ggplot2::theme_bw()+
      ggplot2::scale_fill_manual(values=c("brown", "chartreuse4")) + ggplot2::xlab("Group") +
      ggplot2::ylab(paste0("Mean ", exon_type, " \n per gene")) + ggplot2::geom_violin(fill = NA)
  } else {
    p2 <- ggplot2::ggplot(dfASpg, ggplot2::aes(x = type, y = ASpg, fill = type)) + ggplot2::geom_dotplot(binaxis='y', stackdir='center') + ggplot2::theme_bw()+
      ggplot2::scale_fill_manual(values=c("brown", "chartreuse4")) + ggplot2::xlab("Group") +
      ggplot2::ylab(paste0("Mean ", exon_type, " \n per gene"))
  }
  ## PSI distribution across phenotype
  psiVals <- lapply(data_list, function(x) x[,grep("PSI", colnames(x))])
  dfECDF <- data.frame(val = unlist(psiVals),
                       type = unlist(lapply(1:length(psiVals), function(x) rep(unlist(lapply(sample_list, "[[", 2))[x], length(psiVals[[x]])))))

  dfECDF <- dfECDF[dfECDF$val < 1 & dfECDF$val > 0,]
  p3 <- ggplot2::ggplot(dfECDF, ggplot2::aes(x = val, colour = type, fill = type)) +
    ggplot2::stat_ecdf(geom = "step") + ggplot2::theme_bw() + ggplot2::scale_color_manual(breaks=c("control","test"),
                                                                                 values=c("brown", "chartreuse4")) + ggplot2::xlab("PSI") + ggplot2::ylab(paste0(exon_type, " PSI ECDF"))
  comb_plot <- ggpubr::ggarrange(p1, p2, p3, labels = c("A", "B", "C"))

  return(list(p1 = p1, p2 = p2, p3 = p3, comb_plot = comb_plot))
}
