#' Extract differentially included AS events from HIT Index data
#'
#' @param test_names a vector of test_names
#' @param control_names a vector of control_names
#' @param et string of the exon_type
#' @param cores the number of cores requested
#' @param outlier_threshold the thresholding of the cooks distance, no outlier removal is "Inf"
#' @param min_prop_samples min prop of samples to require an event to be identified in
#' @return a dataframe with differential inclusion information
#' @importFrom data.table := data.table fread rbindlist fifelse
#' @importFrom dplyr arrange select
#' @export
differential_inclusion_rMATS <- function(control_names, test_names, et,
                                         outlier_threshold = c("4/n", "1", "Inf")[1],
                                         minReads = 10, min_prop_samples = .5) {


  sample_types <- data.table::data.table(sample_name = c(test_names, control_names),
                                         type = rep(c("test", "control"), c(length(test_names), length(control_names))))

  # Load PSI values for each sample and splicing event type
  psi_data_list <- lapply(1:nrow(sample_types), function(x) {
    fdf <- data.table::fread(paste0(sample_types$sample_name[x], ".", et, "PSI"))
    fdf[,c(-1)]
  })
  names(psi_data_list) <- sample_types$sample_name

  # Bind all data tables together
  psi_data <- data.table::rbindlist(psi_data_list, idcol = "sample_name")
  psi_data <- merge(psi_data, sample_types, by = "sample_name", all.x = TRUE)

  temp <- psi_data
  if (et == "SE") {
    temp <- temp %>% dplyr::select('sample_name', 'GeneID',
                                   "chr", "strand", "exonStart_0base",
                                   "exonEnd", "upstreamES", "upstreamEE",
                                   "downstreamES", "downstreamEE", "IncLevel1",
                                   "IJC_SAMPLE_1", "SJC_SAMPLE_1", 'IncFormLen', 'SkipFormLen', 'type')
    temp$id <- paste0(temp$GeneID, "#", temp$chr, ":", temp$exonStart_0base, "-", temp$exonEnd, "#",
                      temp$strand, ";", temp$upstreamES, "-", temp$upstreamEE, ";",
                      temp$downstreamES, "-", temp$downstreamEE)

  } else if (et == "A3SS" | et == "A5SS") {
    temp <- temp %>% dplyr::select('sample_name', 'GeneID', "chr", "strand",
                                   "longExonStart_0base", "longExonEnd",
                                   "shortES", "shortEE", "flankingES", "flankingEE",
                                   "IncLevel1", "IncLevel2", "IJC_SAMPLE_1",
                                   "SJC_SAMPLE_1", 'IncFormLen', 'SkipFormLen', 'type')
    temp$id <- paste0(temp$GeneID, "#", temp$chr, ":", temp$longExonStart_0base, "-", temp$longExonEnd, "#",
                      temp$strand, ";", temp$shortES, "-", temp$shortEE, ";",
                      temp$flankingES, "-", temp$flankingEE)

  } else if (et == "MXE") {
    colnames(temp)[colnames(temp) %in% c("1stExonStart_0base", "1stExonEnd",
                                         "2ndExonStart_0base", "2ndExonEnd")] <- paste0("X",
                                                                                        colnames(temp)[colnames(temp) %in%
                                                                                                         c("1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd")])
    temp <- temp %>% dplyr::select('sample_name', 'GeneID', "chr", "strand",
                                   "X1stExonStart_0base", "X1stExonEnd", "X2ndExonStart_0base",
                                   "X2ndExonEnd", "upstreamES", "upstreamEE",
                                   "downstreamES", "downstreamEE", "IncLevel1",
                                   "IncLevel2", "IJC_SAMPLE_1", "SJC_SAMPLE_1", 'IncFormLen', 'SkipFormLen', 'type')

    temp[, X1stExonStart_0base := ifelse(strand == "+", X1stExonStart_0base, X2ndExonStart_0base)]
    temp[, X1stExonEnd := ifelse(temp$strand == "+", temp$X1stExonEnd, temp$X2ndExonEnd)]
    temp[, X2ndExonStart_0base := ifelse(temp$strand == "+", temp$X2ndExonStart_0base, temp$X1stExonStart_0base)]
    temp[, X2ndExonEnd := ifelse(temp$strand == "+", temp$X2ndExonEnd, temp$X1stExonEnd)]

    temp$id <- paste0(temp$GeneID, "#", temp$chr, ":", temp$X1stExonStart_0base, "-", temp$X1stExonEnd, "#",
                      temp$strand, ";", temp$X2ndExonStart_0base, "-", temp$X2ndExonEnd, ";",
                      temp$upstreamES, "-", temp$upstreamEE, ";",
                      temp$downstreamES, "-", temp$downstreamEE)

  } else if (et == "RI") {
    temp <- temp %>% dplyr::select('sample_name', 'GeneID', "chr", "strand",
                                   "riExonStart_0base", "riExonEnd", "upstreamES",
                                   "upstreamEE", "downstreamES", "downstreamEE", "IncLevel1",
                                   "IncLevel2", "IJC_SAMPLE_1", "SJC_SAMPLE_1", 'IncFormLen', 'SkipFormLen', 'type')
    temp$id <- paste0(temp$GeneID, "#", temp$chr, ":", temp$riExonStart_0base, "-", temp$riExonEnd, "#",
                      temp$strand, ";", temp$upstreamES, "-", temp$upstreamEE, ";",
                      temp$downstreamES, "-", temp$downstreamEE)

  }

  if (!(length(temp$IncLevel1) == sum(is.na(temp$IncLevel1)))) {
    temp <- temp %>% dplyr::select(sample_name, id, IncLevel1, IJC_SAMPLE_1, SJC_SAMPLE_1, IncFormLen, SkipFormLen, type)
  } else {
    temp <- temp %>% dplyr::select(sample_name, id, IncLevel2, IJC_SAMPLE_2, SJC_SAMPLE_2, IncFormLen, SkipFormLen, type)
  }
  colnames(temp)[3:5] <- c("psi",  "IJC", "SJC") # Rename columns for consistency

  # Calculate additional columns needed for analysis
  temp[, `:=` (
    valid_reads = (IJC >= minReads),
    psi_adjusted = fifelse((IJC >= minReads), psi, 0)
  )]
  temp <- temp[valid_reads == TRUE]


  ## filter by min_prop_size to reduce compute time
  temp <- temp[!is.na(temp$psi),]

  if (et == "SE") {
    keep_ids <- union(
      names(table(temp$id[temp$type == 'control'])[table(temp$id[temp$type == 'control']) >= min_prop_samples*sum(sample_types$type == 'control')]),
      names(table(temp$id[temp$type == 'test'])[table(temp$id[temp$type == 'test']) >= min_prop_samples*sum(sample_types$type == 'test')]))

    temp <- temp[temp$id %in% keep_ids,]
  }

  all_gene_exons <- unique(temp[, .(id)])

  # Ensure each sample has all gene/exon combinations that appear in any sample
  expanded_data <- all_gene_exons[, .(sample_name = sample_types$sample_name), by = .(id)]
  expanded_data <- merge(expanded_data,
                         temp[,c('id', 'SkipFormLen', "IncFormLen")][!duplicated(temp[,c('id', 'SkipFormLen', "IncFormLen")])],
                         by = "id", all.x = T)

  expanded_data <- merge(expanded_data, temp[,-c('SkipFormLen', 'IncFormLen')], by = c("sample_name", "id"), all.x = TRUE)
  expanded_data <- merge(expanded_data[,-c('type')], sample_types, by = "sample_name", all.x = TRUE)

  expanded_data[, median_sum_IJC_SJC := round(median(IJC + SJC, na.rm = TRUE), digits = 0), by = sample_name]
  expanded_data$median_sum_IJC_SJC[is.na(expanded_data$median_sum_IJC_SJC)] <- median(expanded_data$median_sum_IJC_SJC, na.rm = T)

  expanded_data$SJC[is.na(expanded_data$SJC)] <- expanded_data$median_sum_IJC_SJC[is.na(expanded_data$SJC)]
  expanded_data[is.na(psi), psi := 0]  # Fill missing PSI values with 0
  expanded_data[is.na(IJC), IJC := 0]  # Fill missing PSI values with 0

  expanded_data[is.na(psi_adjusted), psi_adjusted := 0]  # Fill missing PSI values with 0
  expanded_data[is.na(valid_reads), valid_reads := TRUE]  # Fill missing PSI values with 0

  psi_data <- expanded_data

  psi_data[, c("psi_adjusted", "IJC", "SJC") := lapply(.SD, function(x) data.table::fifelse(is.na(x), 0, x)),
           .SDcols = c("psi_adjusted", "IJC", "SJC")]

  psi_data[, valid_group := .N > 1 && data.table::uniqueN(type) > 1, by = .(id)]

  psi_data[, inclusion := as.integer(IJC)]
  psi_data[, exclusion := as.integer(SJC)]
  psi_data <- psi_data[psi_data$valid_group,]

  threshold <- switch(outlier_threshold,
                      "4/n" = 4 / nrow(sample_types),
                      "1" = 1,
                      as.numeric(outlier_threshold))


  psi_data[, c("LR_stat", "p.val", "cooks_d") := {
    # Full and Null Models
    full_model <- glm(
      cbind(inclusion, exclusion) ~ type,
      family = binomial(link = "logit"),
      data = .SD
    )

    # Calculate Cook's Distance
    cooks_d <- as.numeric(cooks.distance(full_model))
    cooks_d[is.na(cooks_d)] <- 0
    # Threshold for Cook's Distance
    noninfluential_points <- which(cooks_d <= threshold)

    # Remove influential points and refit models
    cleaned_data <- .SD[noninfluential_points, ]
    if (length(unique(cleaned_data$type)) > 1) {
      full_model_cleaned <- glm(
        cbind(inclusion, exclusion) ~ type,
        family = binomial(link = "logit"),
        data = cleaned_data
      )
      null_model_cleaned <- glm(
        cbind(inclusion, exclusion) ~ 1,
        family = binomial(link = "logit"),
        data = cleaned_data
      )

      # Likelihood Ratio Test
      lrt <- anova(null_model_cleaned, full_model_cleaned, test = "Chisq")
      LR_statistic <- lrt$Deviance[2]
      p.val <- lrt$`Pr(>Chi)`[2]

      # Return results
      .(LR_statistic = as.numeric(LR_statistic), p.val = as.numeric(p.val), cooks_d)
    } else {
      .(LR_statistic = 0, p.val = 1, cooks_d)
    }

  }, by = .(id)]
  psi_data[psi_data$cooks_d <= threshold,]

  psi_data[, `:=` (
    delta.psi = mean(psi_adjusted[type == "test"]) - mean(psi_adjusted[type == "control"]),
    test_average_psi = mean(psi_adjusted[type == "test"]),
    control_average_psi = mean(psi_adjusted[type == "control"])
  ), by = .(id)]

  #Filter data based on zero counts
  psi_data[, zero_count := sum(psi_adjusted == 0), by = .(id)]

  # Filter data based on min proportion of samples per phenotype and return result
  final_data <- psi_data[psi_data$valid_group, .(
    id, p.val, delta.psi, test_average_psi, control_average_psi,
    count_test = sum(type == "test"), count_control = sum(type == "control"), zero_count

  ), by = .(id)]

  # Adjust p-values for multiple testing
  final_data$p.adj <- p.adjust(final_data$p.val, method = "fdr")
  final_data$p.adj[is.na(final_data$p.adj)] <- 1
  final_data$p.val[is.na(final_data$p.val)] <- 1
  final_data$type <- et
  fd2 <- final_data

  idSplit <- strsplit(fd2$id, split = "#")
  fd2[, gene := unlist(lapply(idSplit, "[[", 1))]
  fd2[, exon := unlist(lapply(idSplit, "[[", 2))]
  fd2[, add_inf := unlist(lapply(idSplit, "[[", 3))]

  stats_out <- data.frame(fd2[!duplicated(fd2),c(-1)] %>% dplyr::arrange(.data$gene))

  if (et == "A5SS" | et == "A3SS" | et == "MXE") {
    # extract paired results for naturally paired output
    stats_out <- paired_rMATS_helper(stats_out)
  } else if (et == "SE") {
    stats_out <- stats_out[rep(1:nrow(stats_out), each = 2),]
    stats_out$add_inf <- paste0(stats_out$add_inf, ";", c("seInclusion", "seExclusion"), ";", rep(1:(nrow(stats_out)/2), each = 2))
    stats_out$delta.psi <- stats_out$delta.psi * c(1, -1)
  } else if (et == "RI") {
    stats_out <- stats_out[rep(1:nrow(stats_out), each = 2),]
    stats_out$add_inf <- paste0(stats_out$add_inf, ";", c("riInclusion", "riExclusion"), ";", rep(1:(nrow(stats_out)/2), each = 2))
    stats_out$delta.psi <- stats_out$delta.psi * c(1, -1)
  }

  return(stats_out)

}


paired_rMATS_helper <- function(df, type) {
  swapped_exon <- paste0(unlist(lapply(strsplit(df$exon, split = ":"), "[[", 1)), ":",
                         unlist(lapply(strsplit(df$add_inf, split = ";"), "[[", 2)))
  mod_df <- do.call(rbind, lapply(1:nrow(df), function(x) {
    i_df <- data.frame(gene = c(df$gene[x], df$gene[x]),
                       exon = c(df$exon[x], swapped_exon[x]),
                       type = c(df$type[x], df$type[x]),
                       delta.psi = c(df$delta.psi[x], -1*(df$delta.psi[x])),
                       p.val = c(df$p.val[x], df$p.val[x]),
                       p.adj = c(df$p.adj[x], df$p.adj[x]),
                       control_average_psi = c(df$control_average_psi[x], 1-df$control_average_psi[x]),
                       test_average_psi = c(df$test_average_psi[x], 1-df$test_average_psi[x]),
                       count_test = c(df$count_test[x], df$count_test[x]),
                       count_control = c(df$count_control[x], df$count_control[x]),
                       zero_count = c(df$zero_count[x], df$zero_count[x]),
                       add_inf = c(paste0(swapped_exon[x], ";prim;", x), paste0(df$exon[x], ";sec;", x)))

    i_df
  }))

  return(mod_df)
}
