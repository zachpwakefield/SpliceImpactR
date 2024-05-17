differential_inclusion_rMATS <- function(control_names, test_names,
                                         stat_model_bool = T, outlier_bool = T, et, cores, outlier_threshold,
                                         min_prop_samples = .5, minReads = 10, max_zero_prop = .5) {


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
    temp <- temp %>% dplyr::select('sample_name', 'GeneID', "chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "IncLevel1", "IJC_SAMPLE_1", "SJC_SAMPLE_1", 'type')
    temp$id <- paste0(temp$GeneID, "#", temp$chr, ":", temp$exonStart_0base, "-", temp$exonEnd, "#",
                      temp$strand, ";", temp$upstreamES, "-", temp$upstreamEE, ";",
                      temp$downstreamES, "-", temp$downstreamEE)

  } else if (et == "A3SS" | et == "A5SS") {
    temp <- temp %>% dplyr::select('sample_name', 'GeneID', "chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "SJC_SAMPLE_1", 'type')
    temp$id <- paste0(temp$GeneID, "#", temp$chr, ":", temp$longExonStart_0base, "-", temp$longExonEnd, "#",
                      temp$strand, ";", temp$shortES, "-", temp$shortEE, ";",
                      temp$flankingES, "-", temp$flankingEE)

  } else if (et == "MXE") {
    colnames(temp)[colnames(temp) %in% c("1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd")] <- paste0("X", colnames(temp)[colnames(temp) %in% c("1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd")])
    temp <- temp %>% dplyr::select('sample_name', 'GeneID', "chr", "strand", "X1stExonStart_0base", "X1stExonEnd", "X2ndExonStart_0base", "X2ndExonEnd", "upstreamES", "upstreamEE",
                                   "downstreamES", "downstreamEE", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "SJC_SAMPLE_1", 'type')

    correct_temp <- do.call(rbind, lapply(1:nrow(temp), function(x)
      data.frame(sc_X1start = ifelse(temp$strand[x] == "+", temp$X1stExonStart_0base[x], temp$X2ndExonStart_0base[x]),
                 sc_X1end = ifelse(temp$strand[x] == "+", temp$X1stExonEnd[x], temp$X2ndExonEnd[x]),
                 sc_X2start = ifelse(temp$strand[x] == "+", temp$X2ndExonStart_0base[x], temp$X1stExonStart_0base[x]),
                 sc_X2end = ifelse(temp$strand[x] == "+", temp$X2ndExonEnd[x], temp$X1stExonEnd[x])
      )


    )
    )
    temp$X1stExonStart_0base <- correct_temp$sc_X1start
    temp$X1stExonEnd <- correct_temp$sc_X1end
    temp$X2ndExonStart_0base <- correct_temp$sc_X2start
    temp$X2ndExonEnd <- correct_temp$sc_X2end

    temp$id <- paste0(temp$GeneID, "#", temp$chr, ":", temp$X1stExonStart_0base, "-", temp$X1stExonEnd, "#",
                      temp$strand, ";", temp$X2ndExonStart_0base, "-", temp$X2ndExonEnd, ";",
                      temp$upstreamES, "-", temp$upstreamEE, ";",
                      temp$downstreamES, "-", temp$downstreamEE)

  } else if (et == "RI") {
    temp <- temp %>% dplyr::select('sample_name', 'GeneID', "chr", "strand", "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "IncLevel1",
                                   "IncLevel2", "IJC_SAMPLE_1", "SJC_SAMPLE_1", 'type')
    temp$id <- paste0(temp$GeneID, "#", temp$chr, ":", temp$riExonStart_0base, "-", temp$riExonEnd, "#",
                      temp$strand, ";", temp$upstreamES, "-", temp$upstreamEE, ";",
                      temp$downstreamES, "-", temp$downstreamEE)

  }

  if (!(length(temp$IncLevel1) == sum(is.na(temp$IncLevel1)))) {
    temp <- temp %>% dplyr::select(sample_name, id, IncLevel1, IJC_SAMPLE_1, SJC_SAMPLE_1, type)
  } else {
    temp <- temp %>% dplyr::select(sample_name, id, IncLevel2, IJC_SAMPLE_2, SJC_SAMPLE_2, type)
  }
  colnames(temp)[3:5] <- c("psi",  "IJC", "SJC") # Rename columns for consistency

  # Calculate additional columns needed for analysis
  temp[, `:=` (
    valid_reads = (IJC + SJC >= minReads),
    psi_adjusted = fifelse((IJC + SJC >= minReads), psi, 0)
  )]
  temp <- temp[valid_reads == TRUE]


  temp <- temp[!is.na(temp$psi),]
  all_gene_exons <- unique(temp[, .(id)])

  # Ensure each sample has all gene/exon combinations that appear in any sample
  # expanded_data <- all_gene_exons[, .(sample_name = sample_types$sample_name), by = .(id)]
  # expanded_data <- merge(expanded_data, temp, by = c("sample_name", "id"), all.x = TRUE)
  expanded_data <- merge(temp, sample_types, by = c("sample_name"), all.x = T)

  expanded_data$type <- expanded_data$type.y
  expanded_data$type.y <- NULL
  expanded_data$type.x <- NULL
  expanded_data[is.na(psi), psi := 0]  # Fill missing PSI values with 0
  expanded_data[is.na(IJC), IJC := 0]  # Fill missing PSI values with 0
  expanded_data[is.na(SJC), SJC := 0]  # Fill missing PSI values with 0
  expanded_data[is.na(psi_adjusted), psi_adjusted := 0]  # Fill missing PSI values with 0
  expanded_data[is.na(valid_reads), valid_reads := T]  # Fill missing PSI values with 0

  psi_data <- expanded_data

  # psi_data[, has_nonzero_psi := any(psi_adjusted > 0), by = .(sample_name, gene)]
  # psi_data <- psi_data[psi_data$has_nonzero_psi]

  psi_data[, c("psi_adjusted", "IJC", "SJC") := lapply(.SD, function(x) data.table::fifelse(is.na(x), 0, x)),
           .SDcols = c("psi_adjusted", "IJC", "SJC")]


  # Linear modeling and Cook's distance for outlier detection if enabled
  psi_data[, valid_group := .N > 1 && data.table::uniqueN(type) > 1, by = .(id)]
  if (outlier_bool) {
    psi_data[psi_data$valid_group, cooks_d := {
      model <- lm(psi_adjusted ~ type, data = .SD)
      cooks.distance(model)
    }, by = .(id)]

    psi_data[is.na(psi_data$cooks_d)] <- 0
    # Determine the threshold for outliers
    threshold <- switch(outlier_threshold,
                        "4/n" = 4 / nrow(sample_types),
                        "1" = 1,
                        as.numeric(outlier_threshold))
    psi_data <- psi_data[cooks_d <= threshold & valid_group]
  }
  psi_data[, valid_group := .N > 1 && uniqueN(type) > 1, by = .(id)]
  psi_data[psi_data$valid_group, c("LR_stat", "p.val") := {
    reduced_model <- lm(psi_adjusted ~ 1, data = .SD)
    full_model <- lm(psi_adjusted ~ type, data = .SD)
    reduced_ll <- logLik(reduced_model)
    full_ll <- logLik(full_model)
    LR_statistic <- -2 * (reduced_ll - full_ll)
    p.val <- pchisq(LR_statistic, df = 1, lower.tail = FALSE)
    .(LR_statistic, p.val)
  }, by = .(id)]

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
  fd2 <- final_data[final_data$count_control >= min_prop_samples * sum(sample_types$type == "control") &
                    final_data$count_test >= min_prop_samples * sum(sample_types$type == "test") &
                    final_data$zero_count <= max_zero_prop * nrow(sample_types),]

  idSplit <- strsplit(fd2$id, split = "#")
  fd2[, gene := unlist(lapply(idSplit, "[[", 1))]
  fd2[, exon := unlist(lapply(idSplit, "[[", 2))]
  fd2[, add_inf := unlist(lapply(idSplit, "[[", 3))]

  stats_out <- data.frame(fd2[!duplicated(fd2),c(-1)] %>% dplyr::arrange(.data$gene))

  if (et == "A5SS" | et == "A3SS" | et == "MXE") {
    # extract paired results for naturally paired output
    stats_out <- paired_rMATS_helper2(stats_out)
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


paired_rMATS_helper2 <- function(df, type) {
  swapped_exon <- paste0(unlist(lapply(strsplit(df$exon, split = ":"), "[[", 1)), ":", unlist(lapply(strsplit(df$add_inf, split = ";"), "[[", 2)))
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
