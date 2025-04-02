#' Get comparison between phenotypes: count of AS event, count of AS event per gene, ECDF compare and save to Foreground subdir
#'
#' @param exon_type type of AS being investigated
#' @param output_location location where everything is being saved, NULL if not desired
#' @param data_df data_df made in specific format, used elsewhere
#' @param minReads to use while filtering
#' @importFrom data.table fread
#' @importFrom ggpubr stat_compare_means
#' @importFrom stats ks.test
#' @return 4 individual plots and 1 combined plot.
#'
#' @examples
#'
#'  pdir <- system.file("extdata", package="SpliceImpactR")
#'  dataDirectory <- paste0(pdir, "/rawData/")
#'  test_group <- paste0(dataDirectory, c("test1","test2", "test3"))
#'  control_group <- paste0(dataDirectory, c("control1", "control2", "control3"))
#'  data_df <- data.frame(sample_names = c(control_group, test_group),
#'                        phenotype_names = c(rep("control", length(control_group)),
#'                        rep("test", length(test_group))))
#'
#'  data_df$utc <- "control"
#'  data_df$utc[data_df$phenotype_names == unique(data_df$phenotype_names)[2]] <- "test"
#'  overview <- getOverviewComparison(data_df, "AFE", output_location = NULL, minReads = 10)
#'
#' @export
getOverviewComparison <- function(data_df, exon_type, output_location, minReads = 10) {
  sample_types <- list()

  for (i in seq_len(nrow(data_df))) {
    sample_types <- c(sample_types, list(c(data_df$sample_names[i],
                                           data_df$utc[i], data_df$phenotype_names[i])))
  }

  # Sort samples by type (control then test)
  sample_list <- c(sample_types[which(unlist(lapply(sample_types, "[[", 2)) == "control")],
                   sample_types[which(unlist(lapply(sample_types, "[[", 2)) == "test")])

  # Get normalization vals from exon files
  if (file.exists(paste0(sample_list[[1]][1], paste0(".exon")))) {
    depth <- unlist(lapply(sample_list, function(x) {
      df <- data.table::fread(paste0(x[1], paste0(".exon")), select = c('nUP', 'nDOWN'))
      depth <- mean(df$nUP + df$nDOWN)
      depth
    }))
  } else {
    depth <- rep(1, length(sample_list))
  }

  # Load PSI values for each sample and splicing event type
  data_list <- lapply(sample_list, function(x) {
    df <- data.frame(fread(paste0(x[1], paste0(".", exon_type, "PSI"))))
    if (exon_type %in% c("AFE", "ALE", "HFE", "HLE")) {
      df[((df$nUP + df$nDOWN) >= minReads),]
    } else {
      df[(df$IJC_SAMPLE_1  >= minReads),]
    }

  })

  ## Number of AS across phenotype
  if (exon_type %in% c("AFE", "ALE", "HFE", "HLE")) {
    colNameInc <- "PSI"
    nASE <- lapply(data_list, function(x) {
      vals <- x[,grep(colNameInc, colnames(x), fixed = TRUE)]
      vals <- vals[!is.na(vals)]
      length(vals)
    })
  } else {
    colNameInc <- "IncLevel1"
    nASE <- lapply(data_list, function(x) {
      vals <- x[,grep(colNameInc, colnames(x), fixed = TRUE)]
      vals <- vals[!is.na(vals)]
      length(vals[vals > 0 & vals < 1])
    })
  }


  dfCount <- data.frame(AScount = unlist(nASE),
                        type = unlist(lapply(sample_list, "[[", 3)))
  dfCount$normAScount <- dfCount$AScount/depth

  if (length(sample_list) > 10) {
    p1 <- ggplot2::ggplot(dfCount, ggplot2::aes(x = type, y = normAScount, fill = type)) +
      ggplot2::geom_violin() + ggplot2::theme_bw()+
      ggplot2::scale_fill_manual(values=c("brown", "chartreuse4"))+ ggplot2::xlab("Group") +
      ggplot2::ylab(paste0("Count of ", exon_type)) +ggplot2::geom_violin(fill = NA) +
      ggpubr::stat_compare_means(method = "wilcox")
  } else {
    p1 <- ggplot2::ggplot(dfCount, ggplot2::aes(x = type, y = normAScount, fill = type)) +
      ggplot2::geom_violin() + ggplot2::theme_bw()+
      ggplot2::scale_fill_manual(values=c("brown", "chartreuse4"))+ ggplot2::xlab("Group") +
      ggplot2::ylab(paste0("Count of ", exon_type)) +
      ggpubr::stat_compare_means(method = "wilcox")

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
                       type = unlist(lapply(sample_list, "[[", 3)))
  dfASpg$normASpg <- dfASpg$ASpg/depth
  if (length(sample_list) > 10) {
    p2 <- ggplot2::ggplot(dfASpg, ggplot2::aes(x = type, y = normASpg, fill = type)) +
      ggplot2::geom_violin() +
      ggplot2::theme_bw() +
      ggplot2::scale_fill_manual(values=c("brown", "chartreuse4")) +
      ggplot2::xlab("Group") +
      ggplot2::ylab(paste0("Mean ", exon_type, " \n per gene")) +
      ggplot2::geom_violin(fill = NA) +
      ggpubr::stat_compare_means(method = 'wilcox')
  } else {
    p2 <- ggplot2::ggplot(dfASpg, ggplot2::aes(x = type, y = normASpg, fill = type)) +
      ggplot2::geom_dotplot(binaxis='y', stackdir='center') +
      ggplot2::theme_bw() +
      ggplot2::scale_fill_manual(values=c("brown", "chartreuse4")) +
      ggplot2::xlab("Group") +
      ggplot2::ylab(paste0("Mean ", exon_type, " \n per gene")) +
      ggpubr::stat_compare_means(method = 'wilcox')
  }

  if (exon_type %in% c("AFE", "ALE", "HFE", "HLE")) {
    if (exon_type %in% c("AFE", "HFE")) {
      idForDist <- lapply(data_list, function(x) {
        y <- x[!is.na(x$AFEPSI) & x$AFEPSI > 0 & ((x$nUP + x$nDOWN) > minReads),]
        paste0(y$gene, "#", y$exon)
      })
    } else if (exon_type %in% c("ALE", "HLE")) {
      idForDist <- lapply(data_list, function(x) {
        y <- x[!is.na(x$ALEPSI) & x$ALEPSI > 0 & ((x$nUP + x$nDOWN) > minReads),]
        paste0(y$gene, "#", y$exon)
      })
    }

  } else {

    getRMATSid <- function(temp, et) {
      temp <- data.table::data.table(temp)
      if (et == "SE") {
        temp <- temp %>% dplyr::select('GeneID', "chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "IncLevel1", "IJC_SAMPLE_1", "SJC_SAMPLE_1")
        temp$id <- paste0(temp$GeneID, "#", temp$chr, ":", temp$exonStart_0base, "-", temp$exonEnd, "#",
                          temp$strand, ";", temp$upstreamES, "-", temp$upstreamEE, ";",
                          temp$downstreamES, "-", temp$downstreamEE)

      } else if (et == "A3SS" | et == "A5SS") {
        temp <- temp %>% dplyr::select('GeneID', "chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "SJC_SAMPLE_1")
        temp$id <- paste0(temp$GeneID, "#", temp$chr, ":", temp$longExonStart_0base, "-", temp$longExonEnd, "#",
                          temp$strand, ";", temp$shortES, "-", temp$shortEE, ";",
                          temp$flankingES, "-", temp$flankingEE)

      } else if (et == "MXE") {
        colnames(temp)[colnames(temp) %in% c("1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd")] <- paste0("X", colnames(temp)[colnames(temp) %in% c("1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd")])
        temp <- temp %>% dplyr::select('GeneID', "chr", "strand", "X1stExonStart_0base", "X1stExonEnd", "X2ndExonStart_0base", "X2ndExonEnd", "upstreamES", "upstreamEE",
                                       "downstreamES", "downstreamEE", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "SJC_SAMPLE_1")

        temp[, X1stExonStart_0base := ifelse(strand == "+", X1stExonStart_0base, X2ndExonStart_0base)]
        temp[, X1stExonEnd := ifelse(temp$strand == "+", temp$X1stExonEnd, temp$X2ndExonEnd)]
        temp[, X2ndExonStart_0base := ifelse(temp$strand == "+", temp$X2ndExonStart_0base, temp$X1stExonStart_0base)]
        temp[, X2ndExonEnd := ifelse(temp$strand == "+", temp$X2ndExonEnd, temp$X1stExonEnd)]

        temp$id <- paste0(temp$GeneID, "#", temp$chr, ":", temp$X1stExonStart_0base, "-", temp$X1stExonEnd, "#",
                          temp$strand, ";", temp$X2ndExonStart_0base, "-", temp$X2ndExonEnd, ";",
                          temp$upstreamES, "-", temp$upstreamEE, ";",
                          temp$downstreamES, "-", temp$downstreamEE)

      } else if (et == "RI") {
        temp <- temp %>% dplyr::select('GeneID', "chr", "strand", "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "IncLevel1",
                                       "IncLevel2", "IJC_SAMPLE_1", "SJC_SAMPLE_1")
        temp$id <- paste0(temp$GeneID, "#", temp$chr, ":", temp$riExonStart_0base, "-", temp$riExonEnd, "#",
                          temp$strand, ";", temp$upstreamES, "-", temp$upstreamEE, ";",
                          temp$downstreamES, "-", temp$downstreamEE)

      }
      return(data.frame(temp))
    }


    idForDist <- lapply(data_list, function(x) {
      y <- x[!is.na(x$IncLevel1) & x$IncLevel1 > 0.05 & x$IncLevel1 < .95 & ((x$IJC_SAMPLE_1 + x$SJC_SAMPLE_1) > minReads),]
      rid <- getRMATSid(temp = y, et = exon_type)$id
      unique(unlist(lapply(strsplit(rid, "#[+-]"), "[[", 1)))
    })
  }

  distSummTest <- as.numeric(table(unlist(lapply(strsplit(
    unique(unlist(idForDist[unlist(lapply(sample_list, "[[", 2)) == "test"])), split = "[.]"), "[[", 1))))
  distSummControl <- as.numeric(table(unlist(lapply(strsplit(
    unique(unlist(idForDist[unlist(lapply(sample_list, "[[", 2)) == "control"])), split = "[.]"), "[[", 1))))

  dfdASpg <- data.frame(ASpg = c(distSummControl, distSummTest),
                        type = c(rep("control", length(distSummControl)), rep("test", length(distSummTest))),
                        type_label = c(rep(unique(data_df$phenotype_names[data_df$utc == "control"]), length(distSummControl)),
                                       rep(unique(data_df$phenotype_names[data_df$utc == "test"]), length(distSummTest))))


  p3 <- ggplot2::ggplot(dfdASpg, ggplot2::aes(y = .data$ASpg, fill = .data$type_label)) +
    ggplot2::geom_histogram(binwidth = 1) + ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(breaks=seq(1,max(dfdASpg$ASpg), max(floor(max(dfdASpg$ASpg)/5), 1)))+
    ggplot2::scale_x_continuous(breaks=seq(0,
                                           max(c(as.integer(table(dfdASpg$ASpg[dfdASpg$type == "control"])),
                                                 as.integer(table(dfdASpg$ASpg[dfdASpg$type == "test"])))),
                                           floor(max(c(as.integer(table(dfdASpg$ASpg[dfdASpg$type == "control"])),
                                                       as.integer(table(dfdASpg$ASpg[dfdASpg$type == "test"]))))/2))) +
    ggplot2::facet_wrap(ggplot2::vars(.data$type_label), strip.position = "bottom") +
    ggplot2::scale_fill_manual(values=c("brown", "chartreuse4")) +
    ggplot2::xlab(paste0("Gene count")) +
    ggplot2::ylab(paste0(exon_type, " count per gene ")) +
    ggplot2::theme(strip.background=ggplot2::element_rect(colour="black",
                                                 fill="white"),
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))

  ## PSI distribution across phenotype
  psiVals <- lapply(data_list, function(x) x[,grep(colNameInc, colnames(x))])
  dfECDF <- data.frame(val = unlist(psiVals),
                       type = unlist(lapply(seq_along(psiVals), function(x) rep(unlist(lapply(sample_list, "[[", 2))[x], length(psiVals[[x]])))))

  dfECDF <- dfECDF[dfECDF$val < 1 & dfECDF$val > 0 & !is.na(dfECDF$val),]
  dfECDF <- dfECDF[!is.na(dfECDF$val),]
  p4 <- ggplot2::ggplot(dfECDF, ggplot2::aes(x = val, colour = type, fill = type)) +
    ggplot2::stat_ecdf(geom = "step") +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(breaks=c("control","test"), values=c("brown", "chartreuse4")) +
    ggplot2::xlab("PSI") +
    ggplot2::ylab(paste0(exon_type, " PSI eCDF")) +
    ggplot2::annotate(
      "text", x = 0.7, y = 0.3,
      label = paste0("K-S p-value = ", signif(stats::ks.test(dfECDF$val[dfECDF$type == "control"],
                                                dfECDF$val[dfECDF$type == "test"])$p.value, 3))
    )

  comb_plot <- ggpubr::ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"),
                                 common.legend = TRUE, legend="bottom")

  if (!(is.null(output_location))) {
    pdf(paste0(output_location, "Foreground/", "comparison_plots.pdf"))
    print(comb_plot)
    dev.off()
  }




  return(list(p1 = p1, p2 = p2, p3 = p3, comb_plot = comb_plot))
}
