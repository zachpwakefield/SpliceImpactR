#' Extract significance stats for each aRNAp
#'
#' @param chosen_method stats method selected
#' @param psi_data_sf a df with exons for signifcance calc
#' @return a dataframe with psi values, LR, cooks d
#' @importFrom data.table :=
#' @export
getSignificance <- function(psi_data_sf, chosen_method) {
  if (chosen_method == "nbGLM") {
    psi_data_sf[, c("LR_stat", "p.val", "cooks_d") := {

      # All data for the current gene
      gene_data <- .SD

      # Unique exons within the gene
      exon_ids <- unique(gene_data$exon)

      # Initialize lists to store results for each exon
      LR_stats <- numeric(nrow(gene_data))
      p_vals <- numeric(nrow(gene_data))
      cooks_ds <- numeric(nrow(gene_data))

      # Iterate over each exon
      for (exon_id in exon_ids) {

        # Apply the nbGLM_model function
        model_results <- suppressWarnings(
          nbGLM_model(
            exon_data = gene_data[exon == exon_id, ],
            gene_data = gene_data,
            exon_of_interest = exon_id,
            threshold = Inf
          )
        )

        # Identify the rows corresponding to the current exon
        exon_rows <- gene_data$exon == exon_id

        # Assign the results to the corresponding rows
        LR_stats[exon_rows] <- model_results$LR_statistic
        p_vals[exon_rows] <- model_results$p.val
        cooks_ds[exon_rows] <- model_results$cooks_d
      }

      # Return a list of vectors to assign
      list(LR_stat = LR_stats, p.val = p_vals, cooks_d = cooks_ds)

    }, by = gene]

  } else if (chosen_method == "qbGLM") {
    psi_data_sf[, c("LR_stat", "p.val", "cooks_d") := qbGLM_model(.SD, threshold = Inf),
                by = .(gene, exon)]

  } else if (chosen_method == "wilcox") {
    psi_data_sf[, c("LR_stat", "p.val", "cooks_d") := wilcox_model(.SD),
                by = .(gene, exon)]
  }
  return(psi_data_sf)
}



#' Extract size factors for each sample
#'
#' @param df dataframe with key info for size factor calc
#' @return a dataframe with sizeFactors added
#' @importFrom data.table :=
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr group_by summarize left_join
#' @export
getSizeFactors <- function(df) {
  df$reads <- df$nUP + df$nDOWN
  counts_wide <- df %>%
    dplyr::group_by(sample_name, exon) %>%
    dplyr::summarize(reads = sum(reads), .groups = "drop") %>%  # aggregate duplicates
    tidyr::pivot_wider(
      names_from  = sample_name,   # columns become sample names
      values_from = reads,
      values_fill = 0              # fill missing combos with 0
    )
  exon_ids <- counts_wide$exon

  # Drop the 'exon' column
  counts_mat <- as.matrix(counts_wide[, -1])

  rownames(counts_mat) <- exon_ids

  geoMean <- exp(rowMeans(log(counts_mat+1)))

  sizeFactors_list <- apply(counts_mat, 2, function(col_j) {
    median((col_j+1)/geoMean)
  })

  df <- dplyr::left_join(df, data.frame(sample_name = names(sizeFactors_list),
                                        sizeFactor = as.numeric(sizeFactors_list)))
  return(df)
}

#' Extract significance stats for each aRNAp
#'
#' @param data dataframe with key info for diff inc calc
#' @param threshold threshold for cooks.d
#' @return p vals, LR vals, cooks d vals
#' @importFrom stats cooks.distance glm anova quasibinomial
#' @export
qbGLM_model <- function(data, threshold = Inf) {
  # Ensure necessary columns are present
  required_cols <- c("psi", "type", "total")
  if (!all(required_cols %in% names(data))) {
    stop("Data must contain the following columns: psi, type, total")
  }

  # Initialize default results
  result <- list(
    LR_statistic = 0,
    p.val        = 1,
    cooks_d      = rep(1, nrow(data))
  )

  # Fit the full model
  tryCatch({
    full_model <- stats::glm(
      psi ~ type,
      weights = total,
      family = stats::quasibinomial,
      data = data,
      control = stats::glm.control(maxit = 1000)
    )

    # Calculate Cook's Distance
    cooks_d <- as.numeric(stats::cooks.distance(full_model))

    # Identify non-influential points
    noninfluential_points <- which(cooks_d <= threshold)

    # Clean the data by removing influential points
    cleaned_data <- data[noninfluential_points, ]

    # Check if there is sufficient variation in 'type'
    if (length(unique(cleaned_data$type)) > 1) {
      # Refit the full and null models on cleaned data
      full_model_cleaned <- stats::glm(
        psi ~ type,
        weights = total,
        family = stats::quasibinomial,
        data = cleaned_data,
        control = stats::glm.control(maxit = 1000)
      )

      null_model_cleaned <- stats::glm(
        psi ~ 1,
        weights = total,
        family = stats::quasibinomial,
        data = cleaned_data,
        control = stats::glm.control(maxit = 1000)
      )

      # Perform Likelihood Ratio Test
      lrt <- stats::anova(null_model_cleaned, full_model_cleaned, test = "Chisq")
      LR_statistic <- lrt$Deviance[2]
      p_val <- lrt$`Pr(>Chi)`[2]

      # Assign results
      result <- list(
        LR_statistic = as.numeric(LR_statistic),
        p.val        = as.numeric(p_val),
        cooks_d      = cooks_d
      )
    }

  }, error = function(e) {
    message(paste("qbGLM_model failed:", e$message))
    # Defaults are already set
  })

  return(result)
}

#' Extract significance stats for each aRNAp
#'
#' @param data dataframe with key info for diff inc calc
#' @return p vals for each event
#' @importFrom stats wilcox.test
#' @export
wilcox_model <- function(data) {
  # Ensure necessary columns are present
  required_cols <- c("psi", "type")
  if (!all(required_cols %in% names(data))) {
    stop("Data must contain the following columns: psi, type")
  }

  # Initialize default results
  result <- list(
    LR_statistic = 0,
    p.val        = 1,
    cooks_d      = 1
  )

  # Check for sufficient variation and difference in means
  if (length(unique(data$type)) > 1 &&
      mean(data$psi[data$type == 'test'], na.rm = TRUE) !=
      mean(data$psi[data$type == 'control'], na.rm = TRUE)) {

    # Perform Wilcoxon test
    wilcox <- tryCatch({
      stats::wilcox.test(psi ~ type, exact = FALSE, data = data)
    }, error = function(e) {
      message(paste("Wilcoxon test failed:", e$message))
      return(NULL)
    })

    if (!is.null(wilcox)) {
      result$p.val <- wilcox$p.value
    }
  }

  return(result)
}

#' Extract significance stats for each aRNAp
#'
#' @param data dataframe with key info for diff inc calc
#' @param threshold threshold for cooks.d
#' @return p vals, LR vals, cooks d vals
#' @importFrom stats cooks.distance glm anova glm.control
#' @importFrom MASS glm.nb
#' @export
nbGLM_model <- function(exon_data, gene_data, exon_of_interest, threshold = Inf) {

  # Check if multiple exons exist and if each exon has
  if (length(unique(gene_data$exon)) > 1 &
      length(unique(exon_data$type)) > 1) {
    # Ensure factors are properly set
    gene_data$sample_name <- factor(gene_data$sample_name)
    gene_data$type        <- factor(gene_data$type)
    gene_data$exon_i      <- as.factor(ifelse(gene_data$exon == exon_of_interest, 1, 0))

    # Initialize result list with default values
    result <- list(
      LR_statistic = rep(3, sum(gene_data$exon == exon_of_interest)),
      p.val        = rep(3, sum(gene_data$exon == exon_of_interest)),
      cooks_d      = rep(3, sum(gene_data$exon == exon_of_interest))
    )

    # Try fitting the full model
    tryCatch({
      full_model <- MASS::glm.nb(
        formula = nDiff ~ sample_name + exon_i + type:exon_i + offset(log(sizeFactor)),
        data    = gene_data,
        control = stats::glm.control(maxit = 1000)
      )

      # Calculate Cook's Distance
      cooks_d <- as.numeric(stats::cooks.distance(full_model))

      # Threshold for Cook's Distance
      noninfluential_points <- which(!is.na(cooks_d) & cooks_d <= threshold)

      # Remove influential points and refit models
      cleaned_data <- gene_data[noninfluential_points, ]

      # Check if cleaned_data has sufficient variation
      if (length(unique(cleaned_data$type)) > 1) {
        full_model_cleaned <-  MASS::glm.nb(
          formula = nDiff ~ sample_name + exon_i + type:exon_i + offset(log(sizeFactor)),
          data    = cleaned_data,
          control = stats::glm.control(maxit = 1000)
        )



        null_model_cleaned <- MASS::glm.nb(
          formula = nDiff ~ sample_name + exon_i + offset(log(sizeFactor)),
          data    = cleaned_data,
          control = stats::glm.control(maxit = 1000)
        )

        # Perform Likelihood Ratio Test
        lrt <- stats::anova(null_model_cleaned, full_model_cleaned, test = "Chisq")
        LR_statistic <- lrt$`LR stat.`[2]
        dexseq_pvals <- lrt$`Pr(Chi)`[2]

        # Extract Cook's Distance for the exon of interest
        cooks_d_length_adjusted <- cooks_d[gene_data$exon == exon_of_interest]

        # Populate the result list
        result <- list(
          LR_statistic = rep(as.numeric(LR_statistic), length(cooks_d_length_adjusted)),
          p.val        = rep(as.numeric(dexseq_pvals), length(cooks_d_length_adjusted)),
          cooks_d      = cooks_d_length_adjusted
        )
      }

    }, error = function(e) {
      # On error, return default result with -1
      message(paste("Model fitting failed for exon:", exon_of_interest, "\nError:", e$message))
      # The default result is already set
    })

  } else {
    # If only one exon, return default results
    repLength <- sum(gene_data$exon == exon_of_interest)
    result <- list(
      LR_statistic = rep(2, repLength),
      p.val        = rep(2, repLength),
      cooks_d      = rep(2, repLength)
    )
  }

  return(result)
}



