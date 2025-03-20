

# branch <- '/projectnb2/evolution/zwakefield/GTEx_proteinImpact/results/brain_cerebellum_v_heart_left-ventricle/AFE/DomainEnrichment/(-)domainEnrichment.csv'
# events <- c('AFE', 'HFE', 'SE', 'A5SS', 'A3SS', 'MXE', 'RI', 'HLE', 'ALE')
#
# pfg_list <- lapply(events, function(x) {
#   readr::read_csv(paste0(branch, x, '/pairedOutput/paired_combined_rows.csv'))
# })
# names(pfg_list) <- events
#
#
# fg_list <- lapply(events, function(x) {
#   readr::read_csv(paste0(branch, x, '/Foreground/fglfc.csv'))
# })
# names(fg_list) <- events
#
# fg_list <- lapply(events, function(x) {
#   readr::read_csv(paste0(branch, x, '/DomainEnrichment/fglfc.csv'))
# })
# names(fg_list) <- events
#
#
# domain_list <- lapply(events, function(x) {
#   if (file.exists(paste0(branch, x, '/DomainEnrichment/(-)domainEnrichment.csv'))) {
#     down <- readr::read_csv(paste0(branch, x, '/DomainEnrichment/(-)domainEnrichment.csv'))
#     if (file.exists(paste0(branch, x, '/DomainEnrichment/(+)domainEnrichment.csv'))) {
#       out <- rbind(readr::read_csv(paste0(branch, x, '/DomainEnrichment/(+)domainEnrichment.csv')),
#                   down)
#       out$event <- x
#     } else {
#       out <- down
#       out$event <- x
#     }
#   } else {
#     if (file.exists(paste0(branch, x, '/DomainEnrichment/(+)domainEnrichment.csv'))) {
#       out <- readr::read_csv(paste0(branch, x, '/DomainEnrichment/(+)domainEnrichment.csv'))
#       out$event <- x
#     } else {
#      out <- data.frame()
#     }
#   }
#   out
# })
# names(domain_list) <- events

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
#' @importFrom scales percent
#' @importFrom ComplexUpset upset upset_set_size
#' @export
getIntegratedResults <- function(fg_list, pfg_list, domain_list) {

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
      MedianScore = median(pfg_list[[i]]$pMatch[pfg_list[[i]]$pMatch > 0]),
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
    dplyr::mutate(Median = median(.data$AlignmentScores))

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

  # Ensure logical filtering for mutually exclusive events
  gene_wide$HFE[gene_wide$AFE & gene_wide$HFE] <- FALSE
  gene_wide$HLE[gene_wide$ALE & gene_wide$HLE] <- FALSE

  # Define event categories for UpSet plot
  event_types <- colnames(gene_wide)[-1]  # Exclude Genes column

  # Create the ULTIMATE UpSet plot
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
