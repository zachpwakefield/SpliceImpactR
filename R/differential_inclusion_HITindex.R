differential_inclusion_HITindex <- function(test_names, control_names, et, cores = 2,
                                            outlier_threshold = c("4/n", "1", 1)[1],
                                            min_prop_samples = .5,
                                            minReads = 10, max_zero_prop = .5) {


  # Create sample type vector efficiently
  sample_types <- data.table::data.table(sample_name = c(test_names, control_names),
                                         type = rep(c("test", "control"), c(length(test_names), length(control_names))))

  # Load PSI values for each sample and splicing event type
  psi_data_list <- lapply(1:nrow(sample_types), function(x) {
    fdf <- data.table::fread(paste0(sample_types$sample_name[x], ".", et, "PSI"))
    fdf$type <- sample_types$type[x]
    fdf
  })
  psi_data_list <- lapply(paste0(sample_types$sample_name, ".", et, "PSI"), fread)
  names(psi_data_list) <- sample_types$sample_name

  # Bind all data tables together
  psi_data <- data.table::rbindlist(psi_data_list, idcol = "sample_name")
  psi_data <- merge(psi_data, sample_types, by = "sample_name", all.x = TRUE)
  colnames(psi_data)[grep("PSI", colnames(psi_data))] <- 'psi'

  all_gene_exons <- unique(psi_data[, .(gene, exon, strand)])

  # Ensure each sample has all gene/exon combinations that appear in any sample
  expanded_data <- all_gene_exons[, .(sample_name = sample_types$sample_name), by = .(gene, exon, strand)]
  expanded_data <- merge(expanded_data, psi_data, by = c("sample_name", "gene", "exon", "strand"), all.x = TRUE)
  expanded_data <- merge(expanded_data, sample_types, by = "sample_name", all.x = T)
  expanded_data$type <- expanded_data$type.y
  expanded_data$type.y <- NULL
  expanded_data$type.x <- NULL
  expanded_data[is.na(psi), psi := 0]  # Fill missing PSI values with 0
  expanded_data[is.na(nUP), nUP := 0]  # Fill missing PSI values with 0
  expanded_data[is.na(nDOWN), nDOWN := 0]  # Fill missing PSI values with 0
  psi_data <- expanded_data

  # Calculate additional columns needed for analysis
  psi_data[, `:=` (
    nDiff = if (et == "AFE") nDOWN - nUP else nUP - nDOWN,
    valid_reads = (nUP + nDOWN >= minReads),
    psi_adjusted = fifelse((nUP + nDOWN >= minReads), psi, 0)
  )]

  psi_data[, has_nonzero_psi := any(psi_adjusted > 0), by = .(sample_name, gene)]
  psi_data <- psi_data[psi_data$has_nonzero_psi]

  psi_data[, c("psi_adjusted", "nDiff", "nUP", "nDOWN", "HITindex") := lapply(.SD, function(x) data.table::fifelse(is.na(x), 0, x)),
           .SDcols = c("psi_adjusted", "nDiff", "nUP", "nDOWN", "HITindex")]


  # Linear modeling and Cook's distance for outlier detection if enabled
  psi_data[, valid_group := .N > 1 && data.table::uniqueN(type) > 1, by = .(gene, exon)]
  if (outlier_bool) {
    psi_data[psi_data$valid_group, cooks_d := {
      model <- lm(psi_adjusted ~ type, data = .SD)
      cooks.distance(model)
    }, by = .(gene, exon)]

    # Determine the threshold for outliers
    threshold <- switch(outlier_threshold,
                        "4/n" = 4 / nrow(sample_types),
                        "1" = 1,
                        as.numeric(outlier_threshold))
    psi_data <- psi_data[cooks_d <= threshold & valid_group]
  }
  psi_data[, valid_group := .N > 1 && uniqueN(type) > 1, by = .(gene, exon)]
  psi_data[psi_data$valid_group, c("LR_stat", "p.val") := {
    reduced_model <- lm(psi_adjusted ~ 1, data = .SD)
    full_model <- lm(psi_adjusted ~ type, data = .SD)
    reduced_ll <- logLik(reduced_model)
    full_ll <- logLik(full_model)
    LR_statistic <- -2 * (reduced_ll - full_ll)
    p.val <- pchisq(LR_statistic, df = 1, lower.tail = FALSE)
    .(LR_statistic, p.val)
  }, by = .(gene, exon)]

  psi_data[, `:=` (
    delta.psi = mean(psi_adjusted[type == "test"]) - mean(psi_adjusted[type == "control"]),
    test_average_psi = mean(psi_adjusted[type == "test"]),
    control_average_psi = mean(psi_adjusted[type == "control"])
  ), by = .(gene, exon)]


  psi_data[, zero_count := sum(psi_adjusted == 0), by = .(gene, exon)]

  # Filter data based on min proportion of samples per phenotype and return result
  final_data <- psi_data[psi_data$valid_group, .(
    gene, exon, p.val, delta.psi, test_average_psi, control_average_psi,
    count_test = sum(type == "test"), count_control = sum(type == "control"), zero_count

  ), by = .(gene, exon)]

  final_data <- final_data[!duplicated(final_data),-c(1:2)] %>% dplyr::arrange(.data$gene)

  # Adjust p-values for multiple testing
  final_data$p.adj <- p.adjust(final_data$p.val, method = "fdr")
  final_data$p.adj[is.na(final_data$p.adj)] <- 1
  final_data$p.val[is.na(final_data$p.val)] <- 1
  final_data$add_inf <- "none"
  final_data$type <- et
  final_data <- final_data[final_data$count_control >= min_prop_samples * sum(sample_types$type == "control") &
                final_data$count_test >= min_prop_samples * sum(sample_types$type == "test") &
                final_data$zero_count <= max_zero_prop * nrow(sample_types),]
  return(data.frame(final_data))

}
