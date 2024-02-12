## This function searches for specific exon types within a given GTF file, comparing them to a set of exons provided (redExon) based on the exon type specified.
## It uses the Jaccard index to determine the overlap between the exons in the GTF file and the provided exons, filtering based on a minimum overlap threshold.


getTranscriptForeground <- function(gtf = gtf, redExon = redExon, ex_type = exon_type, minOverlap = .5, cores) {

  # Print start message indicating the search for the specified exon type
  print(paste("searching for ", ex_type, "...", sep = ""))

  # Parallel computation for each exon in redExon using multiple cores
  rowOuts <- unlist(parallel::mclapply(1:length(redExon$geneR), mc.cores = cores, function(i) {
    # rc_out <- lapply(1:length(redExon$geneR), function(i) {
    hyb_stat <- "no" # Initialize hybrid status as "no"


    # Determine the limit based on exon type for filtering the GTF data
    lim <- switch(ex_type,
                  "AFE" = "first",
                  "ALE" = "last",
                  "SE" = "internal",
                  "A5SS" = "internal",
                  "A3SS" = "internal",
                  "MXE" = "internal",
                  "RI" = "internal")

    # Filter the GTF data for the current gene and exon type
    gtf_min <- gtf[gtf$geneID == redExon$geneR[i] & gtf$type == "exon" & gtf$classification %in% lim,]

    # Calculate Jaccard index for overlap between GTF exons and the current exon

    if (dim(gtf_min)[1] != 0) {
      inEx <- seq(redExon$start[i], redExon$stop[i])
      gtfEx <- lapply(1:length(gtf_min$geneID), function(x) seq(gtf_min$start[x], gtf_min$stop[x]))
      un <- unlist(lapply(gtfEx, function(x) length(union(inEx, unlist(x)))))
      ins <- unlist(lapply(gtfEx, function(x) length(intersect(inEx, unlist(x)))))
      gtf_min$length_jacc <- ins/lengths(gtfEx)
      gtf_min$jaccard <- ins/un
      gtf_min <- gtf_min[order(-gtf_min[, "jaccard"]),]
    }

    # Determine the best match based on the Jaccard index and additional criteria
    # Proceed through various checks to make sure matches are both identified and valid

    ## If no matches or all bad matches
    if (dim(gtf_min)[1] == 0 | max(gtf_min$jaccard) < minOverlap) {
      rowOuts <- 0
    }

    ## If only one match or one best match
    else if (nrow(gtf_min) == 1 | gtf_min$jaccard[1] > gtf_min$jaccard[2]) {
      rowOuts <- gtf_min$rownum[1]
    }

    ## If multiple best, use length of matched exon
    else if (gtf_min$jaccard[1] == gtf_min$jaccard[2]) {
      c_gtf <- gtf_min[gtf_min$jaccard == max(gtf_min$jaccard),] %>% dplyr::arrange(desc(length_jacc))
      rowOuts <- c_gtf$rownum[1]
    }

    ## Other
    else {rowOuts <- 0}

    rowOuts # Return selected row number
  }))

  # Extract the matched rows from the GTF data based on the identified row numbers
  out_matched <- do.call(rbind, parallel::mclapply(unlist(rowOuts)[unlist(rowOuts) != 0], mc.cores = cores, function(x) gtf[gtf$rownum == x, ]))

  # Add additional information to the matched data for further analysis
  out_matched$input_id <- paste(redExon$geneR, ";", redExon$chr, ":", redExon$start, "-", redExon$stop, sep = "")[unlist(rowOuts) != 0]
  out_matched$delta.psi <- redExon$delta.psi[unlist(rowOuts) != 0]
  out_matched$p.adj <- redExon$p.adj[unlist(rowOuts) != 0]

  matched <- out_matched %>% dplyr::relocate(input_id)

  return(matched = matched) # Return the matched data
}

## This function identifies transcripts from a background dataset that correspond to a specified set of exons,
## based on exon type and minimum overlap criteria. It computes the Jaccard index to determine the best match.

getTranscriptBackground <- function(gtf = gtf, redExon = redExon, ex_type = exon_type, minOverlap = .5, cores) {

  # Print a message indicating the start of the search for the specified exon type
  print(paste("searching for ", ex_type, "...", sep = ""))

  # Filter the GTF data for faster computation based on the specified exon type
  lim <- switch(ex_type,
                "AFE" = "first",
                "ALE" = "last",
                "SE" = "internal",
                "A5SS" = "internal",
                "A3SS" = "internal",
                "MXE" = "internal",
                "RI" = "internal")

  gtf_filtered <- gtf[gtf$type == "exon" & gtf$classification %in% lim,]


  # Parallel computation for each exon in redExon using multiple cores
  results <- unlist(parallel::mclapply(1:nrow(redExon), mc.cores = cores, function(i) {

    # Initiate geneR, start, stop
    geneR <- redExon$geneR[i]
    start <- redExon$start[i]
    stop <- redExon$stop[i]

    # Further filter for computational efficiency to current gene
    gtf_min <- gtf_filtered[gtf_filtered$geneID == geneR,]

    if (nrow(gtf_min) == 0) {
      return(0) # Skip if no relevant gtf entries found
    }

    # Calculate Jaccard index for each gtf entry
    jacc <- calculate_jaccard(start, stop, gtf_min$start, gtf_min$stop)
    gtf_min$length_jacc <- jacc[[2]]
    gtf_min$jaccard <- jacc[[1]]

    # Filter based on minOverlap
    gtf_min <- gtf_min[gtf_min$jaccard >= minOverlap,]

    if (nrow(gtf_min) == 0) {
      return(0) # Skip again if no relevant gtf entries found with jaccard index over the minOverlap
    }

    # Order by Jaccard index to find best match
    gtf_min <- gtf_min[order(-gtf_min[, "jaccard"]),]

    ## If only one match or one best match
    if (nrow(gtf_min) == 1 | gtf_min$jaccard[1] > gtf_min$jaccard[2]) {
      return(gtf_min$rownum[1])
    }
    ## If multiple best, use length of matched exon
    else if (gtf_min$jaccard[1] == gtf_min$jaccard[2]) {
      c_gtf <- gtf_min[gtf_min$jaccard == max(gtf_min$jaccard),] %>% dplyr::arrange(desc(length_jacc))
      return(c_gtf$rownum[1])
    } else {return(0)}

  }))

  # Extract the matched rows from the GTF data based on the identified row numbers
  out_matched <- do.call(rbind, parallel::mclapply(results[results != 0], mc.cores = cores, function(x) gtf[gtf$rownum == x, ]))

  # Add additional information to the matched data for further analysis
  out_matched$input_id <- paste(redExon$geneR, ";", redExon$chr, ":", redExon$start, "-", redExon$stop, sep = "")[results != 0]
  matched <- out_matched %>% dplyr::relocate(input_id)

  return(matched = matched)
}
