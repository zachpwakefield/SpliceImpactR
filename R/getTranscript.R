## This function searches for specific exon types within a given GTF file, comparing them to a set of exons provided (redExon) based on the exon type specified.
## It uses the Jaccard index to determine the overlap between the exons in the GTF file and the provided exons, filtering based on a minimum overlap threshold.


getTranscriptForeground <- function(gtf = gtf, redExon = redExon, ex_type = exon_type, minOverlap = .5, cores = 8) {

  # Print start message indicating the search for the specified exon type
  print(paste("searching for ", ex_type, "...", sep = ""))
  rowOuts <- list()  # Initialize list to store output rows


  # Parallel computation for each exon in redExon using multiple cores
  rc_out <- parallel::mclapply(1:length(redExon$geneR), mc.cores = cores, function(i) {
    # rc_out <- lapply(1:length(redExon$geneR), function(i) {
    hyb_stat <- "no" # Initialize hybrid status as "no"


    # Determine the limit based on exon type for filtering the GTF data
    lim <- switch(ex_type,
                  "AFE" = "first",
                  "ALE" = "last",
                  "SE" = "internal")

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

    ## If no matches
    if (dim(gtf_min)[1] == 0) {
      rowOuts <- 0
    }
    ## If all bad matches
    else if (max(gtf_min$jaccard) < minOverlap) {
      rowOuts <- 0
    }

    ## If only one match
    else if (nrow(gtf_min) == 1) {
      rowOuts <- gtf_min$rownum[1]
    }

    ## If one best match
    else if (gtf_min$jaccard[1] > gtf_min$jaccard[2]) {
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
  })

  # Extract the row numbers from the parallel computation results
  rowOuts <- unlist(rc_out)

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

getTranscriptBackground <- function(gtf = gtf, redExon = redExon, ex_type = exon_type, minOverlap = .5, cores = 8) {

  # Print a message indicating the start of the search for the specified exon type
  print(paste("searching for ", ex_type, "...", sep = ""))

  rowOuts <- list()  # Initialize list to store output rows

  # Parallel computation for each exon in redExon using multiple cores
  rc_out <- parallel::mclapply(1:length(redExon$geneR), mc.cores = cores, function(i) {
    # rc_out <- lapply(1:length(redExon$geneR), function(i) {
    hyb_stat <- "no" # Initialize hybrid status as "no"


    # Filter the GTF data for faster computation based on the specified exon type
    lim <- switch(ex_type,
                  "AFE" = "first",
                  "ALE" = "last",
                  "SE" = "internal")

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

    ## If no matches
    if (dim(gtf_min)[1] == 0) {
      rowOuts <- 0
    }
    ## If all bad matches
    else if (max(gtf_min$jaccard) < minOverlap) {
      rowOuts <- 0
    }

    ## If only one match
    else if (dim(gtf_min)[1] == 1) {
      rowOuts <- gtf_min$rownum[1]
    }

    ## If one best match
    else if (gtf_min$jaccard[1] > gtf_min$jaccard[2]) {
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
  })

  # Extract the row numbers from the parallel computation results
  rowOuts <- unlist(rc_out)

  # Extract the matched rows from the GTF data based on the identified row numbers
  out_matched <- do.call(rbind, parallel::mclapply(unlist(rowOuts)[unlist(rowOuts) != 0], mc.cores = cores, function(x) gtf[gtf$rownum == x, ]))

  # Add additional information to the matched data for further analysis
  out_matched$input_id <- paste(redExon$geneR, ";", redExon$chr, ":", redExon$start, "-", redExon$stop, sep = "")[unlist(rowOuts) != 0]
  matched <- out_matched %>% dplyr::relocate(input_id)

  return(matched = matched)
}
