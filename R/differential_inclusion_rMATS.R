differential_inclusion_rMATS <- function(control_names, test_names, et = "SE", cores, outlier_threshold, min_proportion_samples_per_phenotype = .333) {

  sample_types <- list()

  # Categorize each sample name as 'test' or 'control'
  for (i in test_names) {
    sample_types <- c(sample_types, list(c(i, 'test')))
  }

  for (i in control_names) {
    sample_types <- c(sample_types, list(c(i, 'control')))
  }

  # Sort samples by type (control then test)
  sample_types_sorted <- c(sample_types[which(unlist(lapply(sample_types, "[[", 2)) == "control")], sample_types[which(unlist(lapply(sample_types, "[[", 2)) == "test")])

  # Load PSI values for each sample and splicing event type
  load_output <- lapply(sample_types_sorted, function(x) read.table(paste0(x[1], paste0(".", et, "PSI")), header = T, sep = '\t'))

  # Extract and process rMATS data from loaded PSI values
  rMATS_df <- extract_rMATS(et = et, frM.list = load_output,
                            sample_ids = unlist(lapply(sample_types_sorted, "[[", 1)),
                            cores = cores)

  # Filter out rows with all 0 or all 1 PSI values
  rMATS_df <- rMATS_df[rowSums(rMATS_df[,grep("psi", colnames(rMATS_df))]) > 0 &
                         rowMeans(rMATS_df[,grep("psi", colnames(rMATS_df))]) != 1,]


  # Perform statistical analysis for each row in parallel
  stats_out <- do.call(rbind, parallel::mclapply(1:nrow(rMATS_df), mc.cores = cores, function(x) {
    cR <- rMATS_df[x,]

    # Extract control and test group PSI values
    cont.cR <- cR[,unlist(lapply(control_group, function(y) grep(y, colnames(cR)))), drop = T]
    test.cR <- cR[,unlist(lapply(test_group, function(y) grep(y, colnames(cR)))), drop = T]

    # Calculate Cook's distance to identify outliers
    cont.psi <- as.numeric(cont.cR[grep("psi", names(cont.cR))])
    test.psi <- as.numeric(test.cR[grep("psi", names(test.cR))])
    model <- lm(c(cont.psi, test.psi) ~ c(rep(0, length(cont.psi)), rep(1, length(test.psi))))
    influence <- as.numeric(cooks.distance(model))

    # Determine usable data points based on outlier threshold
    if (outlier_threshold == "4/n") {
      usable <- which(influence <= 4/length(sample_types_sorted))
    } else if (outlier_threshold == "4/mean") {
      usable <- which(influence <= 4/mean(influence))
    } else if (outlier_threshold == "1") {
      usable <- which(influence <= 1)
    } else {
      usable <- which(influence <= outlier_threshold)
    }

    # Apply outlier exclusion for control and test groups
    useIndices_cont <- usable[usable <= length(cont.psi)]
    useIndices_test <- usable[usable > length(test.psi)] - length(cont.psi)
    ol_init <- paste(c(control_group, test_group)[usable], collapse = "#")
    outliers <- ifelse(ol_init == "", "none", ol_init)

    # Recalculate PSI, SJC, IJC, mean values excluding outliers for both phenotypes
    cont.SJC.noOut <- as.numeric(cont.cR[grep("SJC", names(cont.cR))])[useIndices_cont]
    cont.IJC.noOut <- as.numeric(cont.cR[grep("IJC", names(cont.cR))])[useIndices_cont]
    cont.psi.noOut <- as.numeric(cont.cR[grep("psi", names(cont.cR))])[useIndices_cont]

    mean.cont.psi.noOut <- mean(cont.psi.noOut)
    mean.cont.IJC.noOut <- mean(cont.IJC.noOut)
    mean.cont.SJC.noOut <- mean(cont.SJC.noOut)

    test.SJC.noOut <- as.numeric(test.cR[grep("SJC", names(test.cR))])[useIndices_test]
    test.IJC.noOut <- as.numeric(test.cR[grep("IJC", names(test.cR))])[useIndices_test]
    test.psi.noOut <- as.numeric(test.cR[grep("psi", names(test.cR))])[useIndices_test]

    mean.test.psi.noOut <- mean(test.psi.noOut)
    mean.test.IJC.noOut <- mean(test.IJC.noOut)
    mean.test.SJC.noOut <- mean(test.SJC.noOut)

    # Continue with the calculation of delta PSI excluding outliers
    delta.psi <- mean.test.psi.noOut - mean.cont.psi.noOut

    # Perform statistical testing comparing control and test PSI values
    if ((sum(cont.psi != 0) >= min_proportion_samples_per_phenotype*length(control_names)) |
        (sum(test.psi != 0) >= min_proportion_samples_per_phenotype*length(control_names))) {

      data <- data.frame(
        counts = c(cont.IJC.noOut, cont.SJC.noOut, test.IJC.noOut, test.SJC.noOut),
        success = c(rep(1, length(cont.IJC.noOut)), rep(0, length(cont.SJC.noOut)),
                    rep(1, length(test.IJC.noOut)), rep(0, length(test.SJC.noOut))),  # 1 for inclusion, 0 for skipping
        condition = rep(c(rep("cont", length(cont.IJC.noOut)), rep("test", length(test.IJC.noOut))), each = 2)
      )
      data$counts <- as.factor(data$count)

      # Fit the null model (assuming same Psi for both conditions)
      null_model <- glm(counts ~ success, family = binomial(link = "logit"), data = data)

      # Fit the alternative model (allowing different Psi for each condition)
      alt_model <- glm(counts ~ success * condition, family = binomial(link = "logit"), data = data)

      # Compute the likelihood ratio test statistic
      lrt_stat <- -2 * (logLik(null_model) - logLik(alt_model))

      # Degrees of freedom for the test (difference in number of parameters)
      df <- length(coef(alt_model)) - length(coef(null_model))

      # Compute the p-value
      p_value <- pchisq(lrt_stat, df = df, lower.tail = FALSE)
    } else {
      p_value <- -1

    }
    print(x)

    strsplit(rMATS_df$id[x], split = "#")[[1]][1]
    strsplit(rMATS_df$id[x], split = "#")[[1]][2]
    strsplit(rMATS_df$id[x], split = "#")[[1]][3]
    stats_info <- data.frame(t(c(strsplit(rMATS_df$id[x], split = "#")[[1]][1],
                                 strsplit(rMATS_df$id[x], split = "#")[[1]][2],
                                 type = "SE"
                                 delta.psi, p_value,
                                 mean.cont.psi.noOut, mean.test.psi.noOut, outliers,
                                 mean.cont.IJC.noOut, mean.cont.SJC.noOut,
                                 mean.test.IJC.noOut, mean.test.SJC.noOut,
                                 strsplit(rMATS_df$id[x], split = "#")[[1]][3],
                                 influence, c(cont.psi, test.psi))))
    colnames(stats_info) <- c("gene", "exon", "type", "delta.psi", "p.val",
                              "control_average_psi", "test_average_psi", "outlier",
                              "control_average_IJC", "control_average_SJC",
                              "test_average_IJC", "test_average_SJC", "add_inf",
                              paste0(unlist(lapply(sample_types_sorted, "[[", 1)),
                                                                                       "_cooks_d"),
                              paste0(unlist(lapply(sample_types_sorted, "[[", 1)),
                                     "_psi"))


    # Return a data frame with statistical results for each row analyzed
    stats_info

  }
  ))

  # Convert columns to numeric, adjust p-values for multiple testing, and reorder columns
  stats_out <- stats_out %>% dplyr::mutate_at(colnames(stats_out)[c(-1, -2, -3, -c(ncol(stats_out)))], as.numeric)
  stats_out$p.adj <- p.adjust(stats_out$p.val, method = "fdr")
  stats_out <- stats_out %>% dplyr::relocate(p.adj, .after = p.val)
  stats_out$p.adj[stats_out$p.adj < 0] <- -1

  return(stats_out)
}


extract_rMATS <- function(et = "SE", frM.list, sample_ids, cores) {
  all.names <- c() # Initialize a vector to store all names, though it seems unused in this snippet

  # Process each rMATS output file in parallel, extracting relevant columns based on the event type
  rM.vals <- parallel::mclapply(1:length(frM.list), mc.cores = cores, function(i) {
    temp <- frM.list[[i]] # Extract the ith data frame from the list

    # Select columns based on the event type and create a unique identifier for each event
    if (et == "SE") {
      temp <- temp %>% dplyr::select('GeneID', "chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "IncLevel1", "IJC_SAMPLE_1", "SJC_SAMPLE_1")
      temp$id <- paste0(temp$GeneID, "#", temp$chr, ":", temp$exonStart_0base, "-", temp$exonEnd, "#",
                       temp$strand, ";", temp$upstreamES, "-", temp$upstreamEE, ";",
                       temp$downstreamES, "-", temp$downstreamEE)

    } else if (et == "A3SS" | et == "A5SS") {
      temp <- temp %>% dplyr::select('GeneID', "chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE", "IncLevel1", "IncLevel2")
      temp$id <- paste(temp$GeneID, temp$chr, temp$strand, temp$longExonStart_0base, temp$longExonEnd, temp$shortES,
                       temp$shortEE, temp$flankingES, temp$flankingEE,sep = ";")

    } else if (et == "MXE") {
      temp <- temp %>% dplyr::select('GeneID', "chr", "strand", "X1stExonStart_0base", "X1stExonEnd", "X2ndExonStart_0base", "X2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "IncLevel1", "IncLevel2")
      temp$id <- paste(temp$GeneID, temp$chr, temp$strand, temp$X1stExonStart_0base, temp$X1stExonEnd, temp$X2ndExonStart_0base, temp$X2ndExonEnd,temp$upstreamES,
                       temp$upstreamEE, temp$downstreamES, temp$downstreamEE,sep = ";")

    } else if (et == "RI") {
      temp <- temp %>% dplyr::select('GeneID', "chr", "strand", "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "IncLevel1", "IncLevel2")
      temp$id <- paste(temp$GeneID, temp$chr, temp$strand, temp$riExonStart_0base, temp$riExonEnd, temp$upstreamES,
                       temp$upstreamEE, temp$downstreamES, temp$downstreamEE,sep = ";")
    }

    # Select columns for IncLevel and sample inclusion and skipping counts, adjusting column names
    # The selection depends on whether IncLevel1 or IncLevel2 is not NA
    if (!(length(temp$IncLevel1) == sum(is.na(temp$IncLevel1)))) {
      temp <- temp %>% dplyr::select(id, IncLevel1, IJC_SAMPLE_1, SJC_SAMPLE_1)
    } else {
      temp <- temp %>% dplyr::select(id, IncLevel2, IJC_SAMPLE_2, SJC_SAMPLE_2)
    }
    colnames(temp)[2:4] <- c("psi",  "IJC", "SJC") # Rename columns for consistency
    temp <- temp[!(is.na(temp$psi)),] # Remove rows with NA psi values
    temp <- temp %>% dplyr::arrange(id) # Sort by unique identifier
    temp <- temp[!duplicated(temp$id),] # Remove duplicate events
    temp
  })

  # Create a master list of unique event identifiers across all samples
  ids <- unique(unlist(lapply(rM.vals, function(x) x$id)))
  # Initialize a dataframe to handle data
  inco <- data.frame(id = ids, psi = NA, IJC = NA, SJC = NA)

  # For each processed file, merge the event data with the master list to ensure consistent structure
  try2 <- parallel::mclapply(1:length(rM.vals), mc.cores = cores, function(i) {
    hm <- rM.vals[[i]]
    try <- rbind(hm, inco)
    try <- try[!(duplicated(try$id)),] %>% dplyr::arrange(id)
    colnames(try)[2:4] <- paste(sample_ids[i], "_", c("psi", "IJC", "SJC"), sep="")
    try
  })
  comb.df <- do.call(cbind, try2)

  # Remove redundant id columns from the combination
  rM.comb <- comb.df[-c(seq(5, length(colnames(comb.df)), by = 4))]
  rM.comb[is.na(rM.comb)] <- 0

  return(rM.comb) # Return the combined and processed dataframe
}



