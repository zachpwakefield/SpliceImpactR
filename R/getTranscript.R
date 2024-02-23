## This function searches for specific exon types within a given GTF file, comparing them to a set of exons provided (redExon) based on the exon type specified.
## It uses the Jaccard index to determine the overlap between the exons in the GTF file and the provided exons, filtering based on a minimum overlap threshold.


getTranscriptForeground <- function(gtf = gtf, redExon = redExon, ex_type = exon_type, minOverlap = .5, cores) {

  # Print a message indicating the start of the search for the specified exon type
  print(paste("searching for ", ex_type, "...", sep = ""))

  # Parallel computation for each exon in redExon using multiple cores
  results <- matcher(ex_type = ex_type, cores = cores, redExon = redExon,  minOverlap=minOverlap)

  # Double for convenience if not AFE/ALE
  compliment_redExon <- redExon
  # if (ex_type != "SE") {
  #   compliment_redExon <- redExon
  # } else {
  #   compliment_redExon <- redExon[rep(1:nrow(redExon), each = 2),]
  #   add_inf_extension <- rep(1:nrow(redExon), each = 2)
  #
  #   compliment_redExon$add_inf <- paste0(compliment_redExon$add_inf, rep(c("inclusion;", ";exclusion"), nrow(redExon)))
  # }

  # Filter out 0s and ensure results are unique to minimize redundant operations
  valid_results <- results[results != 0]

  # Directly subset the gtf dataframe without parallelization if not necessary
  filtered_gtf <- gtf[gtf$rownum %in% valid_results, ]

  # Ensure gtf is a data.frame or a tibble for this operation
  out_matched <- gtf[match(valid_results, gtf$rownum), ]

  # Add additional information to the matched data for further analysis
  out_matched$input_id <- paste(compliment_redExon$geneR, ";", compliment_redExon$chr, ":", compliment_redExon$start, "-", compliment_redExon$stop, sep = "")[results != 0]
  matched <- out_matched %>% dplyr::relocate(input_id)

  # Add additional information to the matched data for further foreground analysis
  matched$delta.psi <- compliment_redExon$delta.psi[results != 0]
  matched$p.adj <- compliment_redExon$p.adj[results != 0]
  matched$add_inf <- compliment_redExon$add_inf[results != 0]


  return(matched = matched) # Return the matched data
}

## This function identifies transcripts from a background dataset that correspond to a specified set of exons,
## based on exon type and minimum overlap criteria. It computes the Jaccard index to determine the best match.

getTranscriptBackground <- function(gtf = gtf, redExon = redExon, ex_type = exon_type, minOverlap = .5, cores) {

  # Print a message indicating the start of the search for the specified exon type
  print(paste("searching for ", ex_type, "...", sep = ""))

  # Parallel computation for each exon in redExon using multiple cores
  results <- matcher(ex_type = ex_type, background = T, cores = cores, redExon = redExon, minOverlap = minOverlap)

  # Double for convenience if not AFE/ALE
  compliment_redExon <- redExon[rep(1:nrow(redExon), each = 1),]

  # Filter out 0s and ensure results are unique to minimize redundant operations
  valid_results <- results[results != 0]

  # Directly subset the gtf dataframe without parallelization if not necessary
  filtered_gtf <- gtf[gtf$rownum %in% valid_results, ]

  # Ensure gtf is a data.frame or a tibble for this operation
  out_matched <- gtf[match(valid_results, gtf$rownum), ]

  # Add additional information to the matched data for further analysis
  out_matched$input_id <- paste(compliment_redExon$geneR, ";", compliment_redExon$chr, ":", compliment_redExon$start, "-", compliment_redExon$stop, sep = "")[results != 0]
  matched <- out_matched %>% dplyr::relocate(input_id)

  # Add additional information to the matched data for further foreground analysis
  matched$add_inf <- compliment_redExon$add_inf[results != 0]


  return(matched = matched) # Return the matched data
}

