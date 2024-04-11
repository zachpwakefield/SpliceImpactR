#' Get comparison between phenotypes: count of AS event, count of AS event per gene, ECDF compare and save to Foreground subdir
#'
#' @param control_names paths of control samples
#' @param test_names paths of test samples
#' @param exon_type type of AS being investigated
#' @param output_location location where everything is being saved
#' @return 4 individual plots and 1 combined plot.
#' @examples
#' getOverviewComparison(c("path_control1", "path_control2"), c("path_test1", "path_test2"), "AFE", "path_to_output")
getOverviewComparison <- function(control_names, test_names, exon_type, output_location, plot = T) {
  sample_types <- list()

  # Categorize each sample name as 'test' or 'control'
  for (i in test_names) {
    sample_types <- c(sample_types, list(c(i, 'test')))
  }

  for (i in control_names) {
    sample_types <- c(sample_types, list(c(i, 'control')))
  }
  # Sort samples by type (control then test)
  sample_list <- c(sample_types[which(unlist(lapply(sample_types, "[[", 2)) == "control")], sample_types[which(unlist(lapply(sample_types, "[[", 2)) == "test")])

  # Load PSI values for each sample and splicing event type
  data_list <- lapply(sample_list, function(x) read.table(paste0(x[1], paste0(".", exon_type, "PSI")), header = T, sep = '\t'))

  ## Number of AS across phenotype
  if (exon_type %in% c("AFE", "ALE", "HFE", "HLE")) {
    colNameInc <- "PSI"
  } else {
    colNameInc <- "IncLevel1"
  }
  nASE <- lapply(data_list, function(x) {
    vals <- x[,grep(colNameInc, colnames(x), fixed = T)]
    vals <- vals[!is.na(vals)]
    length(vals[vals > 0 & vals < 1])
  })

  dfCount <- data.frame(AScount = unlist(nASE),
                        type = unlist(lapply(sample_list, "[[", 2)))
  if (length(sample_list) > 8) {
    p1 <- ggplot2::ggplot(dfCount, ggplot2::aes(x = type, y = AScount, fill = type)) +
      ggplot2::geom_dotplot(binaxis='y', stackdir='center') + ggplot2::theme_bw()+
      ggplot2::scale_fill_manual(values=c("brown", "chartreuse4"))+ ggplot2::xlab("Group") +
      ggplot2::ylab(paste0("Count of ", exon_type)) +ggplot2::geom_violin(fill = NA)
  } else {
    p1 <- ggplot2::ggplot(dfCount, ggplot2::aes(x = type, y = AScount, fill = type)) +
      ggplot2::geom_dotplot(binaxis='y', stackdir='center') + ggplot2::theme_bw()+
      ggplot2::scale_fill_manual(values=c("brown", "chartreuse4"))+ ggplot2::xlab("Group") +
      ggplot2::ylab(paste0("Count of ", exon_type))
  }

  ## Number of AS per gene across phenotype
  nASpg <- lapply(data_list, function(x) {
    if (exon_type %in% c("AFE", "HFE")) {
      mean(as.numeric(table(x$gene[x$AFEPSI > 0])))
    } else if (exon_type %in% c("ALE", "HLE")) {
      mean(as.numeric(table(x$gene[x$ALEPSI > 0])))
    } else {
      mean(as.numeric(table(x$geneSymbol[x$IncLevel1 > 0])))
    }
  })
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

  ## Distribution of Isoforms
  ## Distribution of Isoforms
  if (exon_type %in% c("AFE", "ALE", "HFE", "HLE")) {
    dASpg <- lapply(data_list, function(x) as.numeric(table(x$gene)))

  } else {
    dASpg <- lapply(data_list, function(x) {
      as.numeric(table(x$geneSymbol[x$IncLevel1 > 0 & x$IncLevel1 < 1]))

    })
  }
  dfdASpg <- data.frame(ASpg = unlist(dASpg),
                        type = unlist(lapply(1:length(dASpg), function(x) rep(unlist(lapply(sample_list, "[[", 2))[x], length(dASpg[[x]])))))

  p3 <- ggplot2::ggplot(dfdASpg, ggplot2::aes(y = .data$ASpg, fill = .data$type)) +
    ggplot2::geom_histogram(binwidth = 1) + ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(breaks=seq(1,max(dfdASpg$ASpg), floor(max(dfdASpg$ASpg)/5)))+
    ggplot2::scale_x_continuous(breaks=seq(0,
                                           max(c(as.integer(table(dfdASpg$ASpg[dfdASpg$type == "control"])),
                                                 as.integer(table(dfdASpg$ASpg[dfdASpg$type == "test"])))),
                                           floor(max(c(as.integer(table(dfdASpg$ASpg[dfdASpg$type == "control"])),
                                                       as.integer(table(dfdASpg$ASpg[dfdASpg$type == "test"]))))/2))) +
    ggplot2::facet_wrap(ggplot2::vars(.data$type), strip.position = "bottom") +
    ggplot2::scale_fill_manual(values=c("brown", "chartreuse4")) +
    ggplot2::xlab(paste0("Gene count")) +
    ggplot2::ylab(paste0(exon_type, " count per gene")) +
    theme(strip.background=ggplot2::element_rect(colour="black",
                                                 fill="white"),
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))

  ## PSI distribution across phenotype
  psiVals <- lapply(data_list, function(x) x[,grep(colNameInc, colnames(x))])
  dfECDF <- data.frame(val = unlist(psiVals),
                       type = unlist(lapply(1:length(psiVals), function(x) rep(unlist(lapply(sample_list, "[[", 2))[x], length(psiVals[[x]])))))

  dfECDF <- dfECDF[dfECDF$val < 1 & dfECDF$val > 0,]
  p4 <- ggplot2::ggplot(dfECDF, ggplot2::aes(x = val, colour = type, fill = type)) +
    ggplot2::stat_ecdf(geom = "step") + ggplot2::theme_bw() + ggplot2::scale_color_manual(breaks=c("control","test"),
                                                                                          values=c("brown", "chartreuse4")) + ggplot2::xlab("PSI") + ggplot2::ylab(paste0(exon_type, " PSI eCDF"))
  comb_plot <- ggpubr::ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"),
                                 common.legend = TRUE, legend="bottom")

  pdf(paste0(output_location, "Foreground/", "comparison_plots.pdf"))
  print(comb_plot)
  dev.off()



  return(list(p1 = p1, p2 = p2, p3 = p3, comb_plot = comb_plot))
}
