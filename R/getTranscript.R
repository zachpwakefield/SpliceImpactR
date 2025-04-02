#' gets the transcripts to match to each as event or exon
#'
#' @param gtf gtf from get_gtf from setup
#' @param redExon input of simplified exon coords internally
#' @param exon_type type of exon being queried
#' @param minOverlap minimum overlap to classify as matched to annotation
#' @param cores number of requested cores
#' @param transcript_gtf output dataframe frmo setup_gtf
#' @return figures and dataframes with paired data
#' @keywords internal
getTranscriptForeground <- function(gtf, redExon, ex_type, minOverlap = .05, cores = 1, transcript_gtf) {

  # Parallel computation for each exon in redExon using multiple cores
  results <- matcher(ex_type, FALSE, cores, redExon,  minOverlap, gtf, transcript_gtf)

  # Double for convenience if not AFE/ALE
  compliment_redExon <- redExon

  # Filter out 0s and ensure results are unique to minimize redundant operations
  valid_results <- results[results != 0]

  # Directly subset the gtf dataframe without parallelization if not necessary
  filtered_gtf <- gtf[gtf$rownum %in% valid_results, ]

  # Ensure gtf is a data.frame or a tibble for this operation
  out_matched <- gtf[match(valid_results, gtf$rownum), ]

  # Add additional information to the matched data for further analysis
  out_matched$input_id <- paste(compliment_redExon$geneR, ";",
                                compliment_redExon$chr, ":",
                                compliment_redExon$start, "-",
                                compliment_redExon$stop, sep = "")[results != 0]
  matched <- out_matched %>% dplyr::relocate(input_id)

  # Add additional information to the matched data for further foreground analysis
  matched$delta.psi <- compliment_redExon$delta.psi[results != 0]
  matched$p.adj <- compliment_redExon$p.adj[results != 0]
  matched$add_inf <- compliment_redExon$add_inf[results != 0]


  return(matched = matched) # Return the matched data
}

#' gets the transcripts to match to each as event or exon
#'
#' @param gtf gtf from get_gtf from setup
#' @param redExon input of simplified exon coords internally
#' @param exon_type type of exon being queried
#' @param minOverlap minimum overlap to classify as matched to annotation
#' @param cores number of requested cores
#' @param transcript_gtf output dataframe frmo setup_gtf
#' @return figures and dataframes with paired data
#' @keywords internal
getTranscriptBackground <- function(gtf, redExon, ex_type, minOverlap = .05, cores = 1, transcript_gtf) {

    # Parallel computation for each exon in redExon using multiple cores
    results <- matcher(ex_type, TRUE, cores, redExon, minOverlap, gtf, transcript_gtf)
    compliment_redExon <- redExon[rep(seq_len(nrow(redExon)), each = 1),]

    # Filter out 0s
    valid_results <- results[results != 0]

    # subset the gtf dataframe
    filtered_gtf <- gtf[gtf$rownum %in% valid_results, ]

    # maintaining order, match valid_results to gtf rownum
    out_matched <- gtf[match(valid_results, gtf$rownum), ]

    # Add additional information to the matched data for further analysis
    out_matched$input_id <- paste(compliment_redExon$geneR, ";",
                                  compliment_redExon$chr, ":",
                                  compliment_redExon$start, "-",
                                  compliment_redExon$stop, sep = "")[results != 0]
    matched <- out_matched %>% dplyr::relocate(input_id)

    # Add additional information to the matched data for further foreground analysis
    matched$add_inf <- compliment_redExon$add_inf[results != 0]


  return(matched = matched) # Return the matched data
}

