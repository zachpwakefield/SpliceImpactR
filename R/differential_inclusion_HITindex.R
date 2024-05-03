## This function conducts differential inclusion analysis for alternative splicing events, specifically AFE (Alternative First Exon) and ALE (Alternative Last Exon).
## It compares splicing event inclusion levels between test and control groups and identifies significant differences using statistical tests.


differential_inclusion_HITindex <- function(test_names, control_names, et, cores = 2, stat_model_bool = T, outlier_bool = T,
                                            outlier_threshold = c("4/n", "4/mean", "1", 1)[1], minReads = 10, min_proportion_samples_per_phenotype = .2) {
  if ((length(test_names) + length(control_names)) < 6 & (stat_model_bool | outlier_bool) ) {
    warning("The stat model and outlier filtering will not have much power with low sample number (eg: <= 6)")
  }

  sample_types <- list()

  # Categorize each sample name as 'test' or 'control'
  for (i in test_names) {
    sample_types <- c(sample_types, list(c(i, 'test')))
  }

  for (i in control_names) {
    sample_types <- c(sample_types, list(c(i, 'control')))
  }

# Process each splicing event type (AFE, ALE) separately

  # Sort samples by type (control then test)
  sample_types_sorted <- c(sample_types[which(unlist(lapply(sample_types, "[[", 2)) == "control")], sample_types[which(unlist(lapply(sample_types, "[[", 2)) == "test")])

  # Load PSI values for each sample and splicing event type
  load_output <- lapply(sample_types_sorted, function(x) read.table(paste0(x[1], paste0(".", et, "PSI")), header = T, sep = '\t'))

  # Extract unique genes and exons
  genes <- unique(unlist(lapply(load_output, function(x) x$gene)))
  exons <- lapply(genes, function(x) unique(unlist(lapply(load_output, function(y) y$exon[y$gene == x]))))
  dic <- exons
  names(dic) <- genes

  # Prepare a dataframe with gene-exon pairs
  df <- data.frame(gene = unlist(lapply(1:length(dic), function(x) rep(names(dic)[x], length(dic[[x]])))),
                   exon = unlist(exons))

  # Identify indices for control and test samples
  ctrl_index <- which(unlist(lapply(sample_types_sorted, "[[", 2)) == "control")
  test_index <- which(unlist(lapply(sample_types_sorted, "[[", 2)) == "test")


  # Extract values for each gene-exon pair and perform linear modeling
  vals <- do.call(rbind, mclapply(1:length(df$gene), mc.cores = cores, function(i) {
    val_extract <- lapply(load_output, function(x) {
      # Extract PSI values and calculate nDiff depending on event type
      if (et == "AFE") {
        psi_holder <- x$AFEPSI[x$gene == df$gene[i] & x$exon == df$exon[i]]
        nDiff_i <- x$nDOWN[x$gene == df$gene[i] & x$exon == df$exon[i]]-x$nUP[x$gene == df$gene[i] & x$exon == df$exon[i]]
      } else if (et == "ALE") {
        psi_holder <- x$ALEPSI[x$gene == df$gene[i] & x$exon == df$exon[i]]
        nDiff_i <- x$nUP[x$gene == df$gene[i] & x$exon == df$exon[i]]-x$nDOWN[x$gene == df$gene[i] & x$exon == df$exon[i]]
      }
      # Return values if available, else return default values
      if (length(psi_holder) > 0) {
        c(psi_holder,
          nDiff_i,
          x$nUP[x$gene == df$gene[i] & x$exon == df$exon[i]],
          x$nDOWN[x$gene == df$gene[i] & x$exon == df$exon[i]],
          x$HITindex[x$gene == df$gene[i] & x$exon == df$exon[i]])
      } else {
        c(0, 0, 0, 0, NA)
      }
    })

    # Aggregate extracted values
    psi <- unlist(lapply(val_extract, "[[", 1))
    nDiff <- unlist(lapply(val_extract, "[[", 2))
    nUP <- unlist(lapply(val_extract, "[[", 3))
    nDOWN <- unlist(lapply(val_extract, "[[", 4))
    HITindex <- unlist(lapply(val_extract, "[[", 5))

    psi[is.na(psi)] <- 0
    nDiff[is.na(nDiff)] <- 0
    nUP[is.na(nUP)] <- 0
    nDOWN[is.na(nDOWN)] <- 0
    HITindex[is.na(HITindex)] <- 0

    # Perform linear regression and calculate Cook's distance for outlier detection
    model <- lm(psi ~ unlist(lapply(sample_types_sorted, "[[", 2)))
    influence <- as.numeric(cooks.distance(model))
    influence[is.na(influence)] <- 0

    # Determine outliers based on the specified threshold
    if (!(outlier_bool)) {
      usable <- which(influence <= max(influence))
    } else if (outlier_threshold == "4/n") {
      usable <- which(influence <= 4/length(sample_types_sorted))
    } else if (outlier_threshold == "4/mean") {
      usable <- which(influence <= 4/mean(influence))
    } else if (outlier_threshold == "1") {
      usable <- which(influence <= 1)
    } else {
      usable <- which(influence <= outlier_threshold)
    }
    condition <- intersect(usable, test_index)
    log_nDiff <- log(nDiff+1)

    # Prepare dataframes with and without outliers for further analysis
    with_outliers_df <- data.frame(gene = df$gene[i],
                                   exon = df$exon[i],
                                   nDiff = nDiff,
                                   log_nDiff = log_nDiff,
                                   psi = psi,
                                   sample = 1:length(sample_types_sorted),
                                   sample_name = unlist(lapply(sample_types_sorted, "[[", 1)),
                                   condition = 0,
                                   cooks_d = influence,
                                   type = et,
                                   nUP = nUP,
                                   nDOWN = nDOWN,
                                   HITindex = HITindex)
    with_outliers_df$condition[condition] <- rep(1, length(condition))
    remove_outliers_df <- with_outliers_df[usable,]
    gdf <- remove_outliers_df
    gdf

  }))

  vals <- vals[vals$nUP + vals$nDOWN >= minReads,]

  # Compute p-values for each exon using a likelihood ratio test
  p_vals <- mclapply(unique(vals$exon), mc.cores = cores, function(x) {
    # Only include exons present in at least 1/3 of each of the phenotypes
    if (sum(vals$exon == x & vals$condition == 1) > min_proportion_samples_per_phenotype*length(test_names) |
        sum(vals$exon == x & vals$condition == 0) > min_proportion_samples_per_phenotype*length(control_names)) {

      # Subset data for the current exon
      gene_exon <- unique(vals$gene[vals$exon == x])
      df <- vals[!is.na(match(vals$exon, x)),]

      # Mark test samples for the current exon
      df$cond2 <- 0
      df$cond2[df$exon == x & df$condition == 1] <- 1
      # df$var <- rep(c(var(df$psi[df$condition == 0]), var(df$psi[df$condition == 1])), each = 2)
      x1 <- df[,c('psi', 'condition')]

      full_model <- lm(psi ~ condition, data = x1)
      full_ll <- logLik(full_model)

      ## reduced model
      reduced_model <- lm(psi ~ 1, data = x1)
      reduced_ll <- logLik(reduced_model)

      # Calculate likelihood ratio statistic and p-value
      LR_statistic <- -2 * (as.numeric(reduced_ll) - as.numeric(full_ll))
      p_val <- pchisq(LR_statistic, df = 1, lower.tail = FALSE)
      ifelse(!is.na(p_val), p_val, -1)} else {-1}

  })

  # Adjust p-values for multiple testing
  p.adj <- p.adjust(p_vals, method = 'fdr')
  p.adj[p.adj < 0] <- -1
  if (!(stat_model_bool)) {
    p.adj <- rep(0.01, length(p.adj))
  }

  # Prepare the final output data
  init_output <- vals[!duplicated(vals[,c('gene', 'exon')]),c('gene', 'exon')]

  # Assemble final data including gene, exon, type, delta.psi, adjusted p-values, and other statistics
  data <- do.call(rbind, mclapply(1:length(unique(vals$exon)), mc.cores = cores, function(h) {
    ex <- init_output[h,2]
    gene_ex <- init_output[h,1]

    vals_min <- vals[vals$gene == gene_ex & vals$exon == ex,]


    dpsi <- mean(vals_min$psi[vals_min$exon == ex & vals_min$condition == 1]) - mean(vals_min$psi[vals_min$exon == ex & vals_min$condition == 0])
    outlier <- paste(setdiff(unlist(lapply(sample_types_sorted, "[[", 1)), vals_min$sample_name[vals_min$exon == ex]), collapse = "#", sep = "#")
    init_df <- data.frame(gene = unique(vals_min$gene[vals_min$exon == ex]),
                          exon = ex,
                          type = unique(vals_min$type[vals_min$exon == ex]),
                          delta.psi = dpsi,
                          p.val = p_vals[h][[1]],
                          p.adj = p.adj[h],
                          control_average_psi = mean(vals_min$psi[vals_min$exon == ex & vals_min$condition == 0]),
                          test_average_psi = mean(vals_min$psi[vals_min$exon == ex & vals_min$condition == 1]),
                          control_average_nDiff = mean(vals_min$nDiff[vals_min$exon == ex & vals_min$condition == 0]),
                          test_average_nDiff = mean(vals_min$nDiff[vals_min$exon == ex & vals_min$condition == 1]),
                          outlier = ifelse(outlier != "" , outlier, "none"),
                          add_inf = "none"
                          , check.names = F)
    add_numerics <- data.frame(t(unlist(lapply(unlist(lapply(sample_types_sorted, "[[", 1)), function(s) {
      ifelse(length(intersect(vals_min$sample_name[vals_min$exon == ex], s)) >=1,
             out <- c(vals_min$psi[vals_min$exon == ex & vals_min$sample_name == s],
                      vals_min$nDiff[vals_min$exon == ex & vals_min$sample_name == s],
                      vals_min$cooks_d[vals_min$exon == ex & vals_min$sample_name == s]),
             out <- c(0, 0, 0))
      out
    }))))
    colnames(add_numerics) <- paste0(rep(unlist(lapply(sample_types_sorted, "[[", 1)), each = 3), c("_psi", "_nDiff", "_cooks_d"))
    ex_df <- cbind(init_df, add_numerics)
  }))

  # Return the assembled data for this splicing event type

  # Return the combined data for all analyzed splicing event types
  return(data)
}

