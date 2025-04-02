#' Perform various integrative analyses on multiple aRNAp results
#'
#' @param fg_list a list of the fglfc.csv output from getForeground, names are the event types
#' @param pfg_list a list of the paired_output from getPaired, names are the event types
#' @param domain_list a list of the data output from getDomainData, names are the event types
#' @return a dataframe with summary stats for each event type supplied along with various plots
#'
#' @importFrom dplyr select filter group_by mutate distinct
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom ggplot2 ggplot aes geom_bar geom_histogram theme_classic theme_minimal theme_bw
#' @importFrom ggplot2 geom_vline scale_fill_manual scale_y_continuous labs coord_flip facet_wrap
#' @importFrom ComplexUpset upset upset_set_size
#' @importFrom stats median
#'
#' @examples
#'
#' pdir <- system.file("extdata", package="SpliceImpactR")
#' dataDirectory <- paste0(pdir, "/")
#' test_group <- paste0(dataDirectory, "rawData/", c("test1","test2", "test3"))
#' control_group <- paste0(dataDirectory, "rawData/", c("control1", "control2", "control3"))
#' data_df <- data.frame(
#'     sample_names = c(control_group, test_group),
#'     phenotype_names = c(
#'       rep("control", length(control_group)),
#'       rep("test", length(test_group))
#'      ),
#'    stringsAsFactors = FALSE
#'   )
#'
#' transDF <- readr::read_csv(paste0(dataDirectory, "transcripts_limited_transDF.csv"))
#' c_trans <- readr::read_lines(paste0(dataDirectory, "transcripts_limited_c_trans.csv"))
#'
#' transcripts_sample <- list(transDF = transDF,
#'                            c_trans = c_trans)
#'
#' gtf_sample <- list(gtf = readr::read_csv(paste0(dataDirectory, "gtf_limited.csv")),
#'             transcript_gtf = readr::read_csv(paste0(dataDirectory, "transcript_gtf_limited.csv")),
#'             tgp_biomart = readr::read_csv(paste0(dataDirectory, "tgp_biomart_limited"))
#'             )
#' translations_sample <- readr::read_lines(paste0(dataDirectory, "translations_limited.csv"))
#'
#' ip <- readr::read_csv(paste0(dataDirectory, "biomart_ip.csv"))
#' code_regions <- readr::read_csv(paste0(dataDirectory, "biomart_code_regions.csv"))
#' pfam_exon_level <- readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon_level.csv"))
#' fsd_exon_data <- readr::read_csv(paste0(dataDirectory, "biomart_data_sample.csv"))
#' pfam_data = readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon.csv"))
#' biomart_data_sample <- list(ip = ip,
#'                      code_regions = code_regions,
#'                      fsd_exon_data = fsd_exon_data,
#'                      pfam_exon_level = pfam_exon_level,
#'                      pfam_data = pfam_data)
#'
#' initTTI <- init_ddi(pdir = dataDirectory,
#'                     output_location = NULL,
#'                     ppidm_class = c("Gold_Standard", "Gold", "Silver", "Bronze")[1],
#'                     removeDups = TRUE,
#'                     cores = 1,
#'                     pfam_data = biomart_data_sample$pfam_data)
#' twoASfullRun <- fullASoutcome(as_types = c("AFE", "SE", "HIT"),
#'                               output_directory = NULL,
#'                               data_directory = dataDirectory,
#'                               data_df,
#'                               outlier_handle = "Inf",
#'                               cutoff = .1,
#'                               cores = 1,
#'                               bg_pre = NA,
#'                               tti_location = NULL,
#'                               initTTI = initTTI,
#'                               mOverlap = .05,
#'                               s_gtf = gtf_sample,
#'                               plotAlignments = FALSE,
#'                               transcripts = transcripts_sample,
#'                               translations = translations_sample,
#'                               biomart_data = biomart_data_sample,
#'                               max_zero_prop = 1,
#'                               min_prop_samples = 0,
#'                               chosen_method = 'qbGLM')
#'
#' fg_list <- list("AFE" = twoASfullRun$AFE$fg$proBed)
#' pfg_list <- list("AFE" = twoASfullRun$AFE$pfg$paired_proBed)
#' domain_list <- list("AFE" = twoASfullRun$AFE$gD$data)
#'
#' integrated <- getIntegratedResults(fg_list,
#'                                    pfg_list,
#'                                    domain_list)
#' @export
#'
#'

getIntegratedResults <- function(fg_list, pfg_list, domain_list) {
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Please install scales to use getIntegratedResults.")
  }
  relativeUseList <- lapply(seq_along(fg_list), function(i) {
    list(
      event = names(fg_list)[i],
      relativeUse = nrow(pfg_list[[i]]) / nrow(fg_list[[i]]),
      Match = sum(pfg_list[[i]]$alignType == "Match") / 2,
      onePC = sum(pfg_list[[i]]$alignType == "onePC"),
      PartialMatch = sum(pfg_list[[i]]$alignType == "PartialMatch"),
      FrameShift = sum(pfg_list[[i]]$alignType == "FrameShift"),
      Rescue = sum(!(pfg_list[[i]]$alignType %in% c("Match", "onePC", "PartialMatch", "FrameShift", "noRescue"))),
      PropRescue = sum(!(pfg_list[[i]]$alignType %in% c("Match", "onePC", "PartialMatch", "FrameShift", "noRescue"))) /
        sum(pfg_list[[i]]$alignType == "FrameShift"),
      MedianScore = stats::median(pfg_list[[i]]$pMatch[pfg_list[[i]]$pMatch > 0]),
      AlignmentScores = pfg_list[[i]]$pMatch[pfg_list[[i]]$pMatch > 0],
      Genes = pfg_list[[i]]$gene,
      Domains = if (nrow(domain_list[[i]])) {
        domain_list[[i]]$domain
      } else {
        NA
      },
      DomainCount = ifelse(nrow(domain_list[[i]]) > 0, sum(domain_list[[i]]$sample_successes), 0),
      RelativeDomainProp = ifelse(
        nrow(domain_list[[i]]) > 0,
        sum(domain_list[[i]]$sample_successes) / sum(pfg_list[[i]]$alignType == "PartialMatch"),
        0
      )
    )
  })

  relativeUse <- data.frame(
    event = unlist(lapply(relativeUseList, function(x) x$event)),
    relativeUse = unlist(lapply(relativeUseList, function(x) x$relativeUse)),
    onePC = unlist(lapply(relativeUseList, function(x) x$onePC)),
    Match = unlist(lapply(relativeUseList, function(x) x$Match)),
    PartialMatch = unlist(lapply(relativeUseList, function(x) x$PartialMatch)),
    FrameShift = unlist(lapply(relativeUseList, function(x) x$FrameShift)),
    Rescue = unlist(lapply(relativeUseList, function(x) x$Rescue)),
    PropRescue = unlist(lapply(relativeUseList, function(x) x$PropRescue)),
    MedianScore = unlist(lapply(relativeUseList, function(x) x$MedianScore)),
    DomainCount = unlist(lapply(relativeUseList, function(x) x$DomainCount)),
    RelativeDomainProp = unlist(lapply(relativeUseList, function(x) x$RelativeDomainProp))
  )

  relativeUse$event <- factor(relativeUse$event, levels = relativeUse$event)

  RelativeUsePlot <- ggplot2::ggplot(relativeUse, ggplot2::aes(x = event, y = relativeUse, fill = "pink")) +
    ggplot2::geom_bar(stat = "identity", fill = "pink") +
    ggplot2::theme_classic()

  MedianAlignScorePlot <- ggplot2::ggplot(relativeUse, ggplot2::aes(x = event, y = MedianScore, fill = "cadetblue4")) +
    ggplot2::geom_bar(stat = "identity", fill = "cadetblue4") +
    ggplot2::theme_classic()

  PropRescuePlot <- ggplot2::ggplot(
    relativeUse[!is.nan(relativeUse$PropRescue) & relativeUse$PropRescue != 0, ],
    ggplot2::aes(x = event, y = PropRescue, fill = "lightblue")
  ) +
    ggplot2::geom_bar(stat = "identity", fill = "lightblue") +
    ggplot2::theme_classic()

  DomainCountPlot <- ggplot2::ggplot(relativeUse, ggplot2::aes(x = event, y = DomainCount, fill = "forestgreen")) +
    ggplot2::geom_bar(stat = "identity", fill = "forestgreen") +
    ggplot2::theme_classic()

  RelativeDomainPlot <- ggplot2::ggplot(relativeUse, ggplot2::aes(x = event, y = RelativeDomainProp, fill = "deeppink4")) +
    ggplot2::geom_bar(stat = "identity", fill = "deeppink4") +
    ggplot2::theme_classic()

  relativeUse_long <- relativeUse %>%
    dplyr::select(event, PartialMatch, onePC, Match, FrameShift) %>%
    tidyr::pivot_longer(cols = -event, names_to = "Category", values_to = "Proportion")

  custom_colors <- c(
    "PartialMatch" = "#E69F00",
    "onePC" = "#56B4E9",
    "Match" = "#009E73",
    "FrameShift" = "#D55E00"
  )

  relativeUse_long$event <- factor(relativeUse_long$event, levels = rev(relativeUse$event))

  pairedEventTypePlot <- ggplot2::ggplot(
    relativeUse_long,
    ggplot2::aes(x = event, y = Proportion, fill = Category)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::labs(
      x = "Event",
      y = "Proportion",
      fill = "Category",
      title = "Proportion of Different Match Types by Event"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::coord_flip()

  alignScoreFrame <- do.call(rbind, lapply(relativeUseList, function(x) {
    data.frame(event = x$event, AlignmentScores = x$AlignmentScores)
  })) %>% dplyr::filter(AlignmentScores > 0 & AlignmentScores < 1)

  dataInt <- alignScoreFrame %>%
    dplyr::group_by(event) %>%
    dplyr::mutate(Median = stats::median(.data$AlignmentScores))

  alignmentScorePlot <- ggplot2::ggplot(
    alignScoreFrame,
    ggplot2::aes(x = AlignmentScores, fill = "deeppink4")
  ) +
    ggplot2::geom_histogram(binwidth = 0.1) +
    ggplot2::scale_fill_manual(values = c("deeppink4")) +
    ggplot2::facet_wrap(ggplot2::vars(event), scales = "free", ncol = 1) +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(
      data = dataInt,
      ggplot2::aes(xintercept = Median),
      linetype = "dashed",
      color = "cadetblue",
      size = 1
    ) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme_classic()

  geneEventDf <- do.call(rbind, lapply(relativeUseList, function(x) {
    data.frame(event = x$event, Genes = x$Genes)
  })) %>% dplyr::distinct()

  gene_wide <- geneEventDf %>%
    dplyr::mutate(present = TRUE) %>%
    tidyr::pivot_wider(names_from = event, values_from = present, values_fill = list(present = FALSE))

  if ('HFE' %in% colnames(gene_wide)) {
    gene_wide$HFE[gene_wide$AFE & gene_wide$HFE] <- FALSE
  }
  if ('HLE' %in% colnames(gene_wide)) {
    gene_wide$HLE[gene_wide$ALE & gene_wide$HLE] <- FALSE
  }



  # Define event categories for UpSet plot
  event_types <- colnames(gene_wide)[-1]  # Exclude Genes column

  if (length(event_types) > 1) {
    geneUpSet <- ComplexUpset::upset(
      gene_wide,
      intersect = event_types,
      name = "Gene Intersection",
      set_sizes = ComplexUpset::upset_set_size(),
      width_ratio = 0.15,
      n_intersections = 20
    )



  domainEventDf <- do.call(rbind, lapply(relativeUseList, function(x) {
    data.frame(event = x$event, Domains = x$Domains)
  })) %>% dplyr::distinct()

  domain_wide <- domainEventDf %>%
    dplyr::mutate(present = TRUE) %>%
    tidyr::pivot_wider(names_from = event, values_from = present, values_fill = list(present = FALSE))

  # Ensure logical filtering for mutually exclusive events
  domain_wide$HFE[domain_wide$AFE & domain_wide$HFE] <- FALSE
  domain_wide$HLE[domain_wide$ALE & domain_wide$HLE] <- FALSE

  # Define event categories for UpSet plot
  event_types <- colnames(domain_wide)[-1]  # Exclude Genes column

  # Create the ULTIMATE UpSet plot
  domainUpSet <- ComplexUpset::upset(
    domain_wide,
    intersect = event_types,
    name = "Gene Intersection",
    set_sizes = ComplexUpset::upset_set_size(),
    width_ratio = 0.15,
    n_intersections = 20
  )
  } else {
    geneUpSet <- "only one type of event given"
    domainUpSet <- "only one type of event given"
  }
  plots <- list(
    RelativeUsePlot = RelativeUsePlot,
    MedianAlignScorePlot = MedianAlignScorePlot,
    PropRescuePlot = PropRescuePlot,
    DomainCountPlot = DomainCountPlot,
    RelativeDomainPlot = RelativeDomainPlot,
    pairedEventTypePlot = pairedEventTypePlot,
    geneUpSet = geneUpSet,
    domainUpSet = domainUpSet,
    alignmentScorePlot = alignmentScorePlot
  )

  # note: summaryTable is not defined in this snippet; please define or remove if not used
  return(list(
    summaryTable = relativeUse,
    plots = plots
  ))
}

