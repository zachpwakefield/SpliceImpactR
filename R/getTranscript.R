getTranscript <- function(gtf = gtf, redExon = redExon, ex_type = exon_type, minOverlap = .5, swaps = F, cores = 8) {
  ## Start up Message
  print(paste("searching for ", ex_type, "...", sep = ""))
  rowOuts <- list()
  checks <- c()
  rc_out <- parallel::mclapply(1:length(redExon$geneR), mc.cores = cores, function(i) {
  # rc_out <- lapply(1:length(redExon$geneR), function(i) {
    hyb_stat <- "no" #HFE" "HLE" "no"


    ## Progress Messagge
    if ((i %% max(1, round(length(redExon[,1])/10, 0))) == 0) {
      print(paste(i, " exons or ", i/length(redExon[,1]), " completed", sep=""))
    }
    ## Reduce gtf for faster computation

    if (ex_type == "AFE") {
      lim <- c("first")
    } else if (ex_type == "ALE") {
      lim <- c("last")
    } else if (ex_type == "SE") {
      lim <- c("internal")
    }
    gtf_min <- gtf[gtf$geneID == redExon$geneR[i] & gtf$type == "exon" & gtf$classification %in% lim,]

    ## Calculate Jaccard
    if (dim(gtf_min)[1] != 0) {
      inEx <- seq(redExon$start[i], redExon$stop[i])
      gtfEx <- lapply(1:length(gtf_min$geneID), function(x) seq(gtf_min$start[x], gtf_min$stop[x]))
      un <- unlist(lapply(gtfEx, function(x) length(union(inEx, unlist(x)))))
      ins <- unlist(lapply(gtfEx, function(x) length(intersect(inEx, unlist(x)))))
      gtf_min$length_jacc <- ins/lengths(gtfEx)
      gtf_min$jaccard <- ins/un
      gtf_min <- gtf_min %>% dplyr::arrange(desc(jaccard))
    }
    ## Proceed through various checks to make sure matches are both identified and valid

    ## If no matches
    if (dim(gtf_min)[1] == 0) {
      checks <- 1
      rowOuts <- 0
    }
    ## If all bad matches
    else if (max(gtf_min$jaccard) < minOverlap) {
      checks <- 2
      rowOuts <- 0
    }

    ## If only one match
    else if (dim(gtf_min)[1] == 1) {
      checks <- 2
      rowOuts <- gtf_min$rownum[1]
    }

    ## If one best match
    else if (gtf_min$jaccard[1] > gtf_min$jaccard[2]) {
      checks <- 3
      rowOuts <- gtf_min$rownum[1]
    }

    ## If multiple best, use length of matched exon
    else if (gtf_min$jaccard[1] == gtf_min$jaccard[2]) {
      checks <- 4
      c_gtf <- gtf_min[gtf_min$jaccard == max(gtf_min$jaccard),] %>% dplyr::arrange(desc(length_jacc))
      rowOuts <- c_gtf$rownum[1]
    }

    ## Other
    else {checks <- 5
    rowOuts <- 0}

    c(rowOuts, checks)
  })
  rowOuts <- unlist(lapply(rc_out, "[[", 1))
  checks <- unlist(lapply(rc_out, "[[", 2))
  ## Remove erroneous matches
  out_matched <- do.call(rbind, parallel::mclapply(unlist(rowOuts)[unlist(rowOuts) != 0], mc.cores = cores, function(x) gtf[gtf$rownum == x, ]))
  ## Data arrangement
  out_matched$input_id <- paste(redExon$geneR, ";", redExon$chr, ":", redExon$start, "-", redExon$stop, sep = "")[unlist(rowOuts) != 0]
  out_matched <- out_matched %>% dplyr::relocate(input_id)
  tot_matched <- out_matched

  ## If looking for swaps, extract paired transcripts
  if (swaps == T) {
    out_matched <- out_matched %>% dplyr::arrange(geneID) %>% dplyr::filter(geneID %in% names(table(out_matched$geneID))[table(out_matched$geneID) == 2])
  }
  return(list(out_matched = out_matched,
              typeBreakdown = checks,
              tot_matched = tot_matched
  ))
}
