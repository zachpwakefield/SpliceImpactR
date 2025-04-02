#' matches the given locations to annotations and transcripts
#'
#' @param background whether making background or not
#' @param redExon input of simplified exon coords internally
#' @param ex_type type of exon being queried
#' @param minOverlap minimum overlap to classify as matched to annotation
#' @param cores number of requested cores
#' @param gtf dataframe from setup_gtf
#' @param gtf_transcript dataframe from setup_gtf
#' @return figures and dataframes with paired data
#' @importFrom parallel mclapply
#' @keywords internal
matcher <- function(ex_type, background = FALSE, cores = 1, redExon, minOverlap=.05, gtf, gtf_transcripts) {

  gtf_exons <- gtf

  # Filter protein_coding_transcripts once
  protein_coding_transcripts <- unique(gtf$transcriptID[gtf$tpc == "protein_coding"
                                                        & !is.na(gtf$tpc)])

  # Pre-compute a lookup for start positions of all transcripts in gtf
  transcript_starts <- setNames(gtf_transcripts$start,
                                gtf_transcripts$transcriptID)

  # Explicit evals
  redExon_globe <- redExon
  transcript_starts_globe <- transcript_starts
  gtf_exons_globe <- gtf_exons
  gtf_transcripts_globe <- gtf_transcripts
  protein_coding_transcripts_globe <- protein_coding_transcripts
  minOverlap_globe <- minOverlap

  if (background) {
    gtf_filtered <- gtf[gtf$classification %in% c("first", "internal", "last"),]
    gtf_filtered_globe <- gtf_filtered  # Explicitly evaluate `gtf_filtered`
    results <- unlist(parallel::mclapply(seq_len(nrow(redExon_globe)), mc.cores = cores, function(i) {
      HITmatcher(i, redExon = redExon_globe,
                 gtf_filtered=gtf_filtered_globe,
                 minOverlap = minOverlap_globe,
                 protein_coding_transcripts = protein_coding_transcripts_globe)
    }))
  } else if (ex_type %in% c("AFE", "ALE", "HFE", "HLE")) {
    gtf_filtered_globe <- gtf[gtf$classification %in% c("first", "internal", "last"),]
    if (ex_type == "HFE" | ex_type == "AFE") {
      lim <- "first"
      gtf_filtered_globe <- gtf[gtf$classification == lim,]
    } else if (ex_type == "HLE" | ex_type == "ALE") {
      lim <- "last"
      gtf_filtered_globe <- gtf[gtf$classification == lim,]
    }
    results <- unlist(parallel::mclapply(seq_len(nrow(redExon_globe)), mc.cores = cores, function(i) {
      HITmatcher(i, redExon = redExon_globe,
                 gtf_filtered=gtf_filtered_globe,
                 minOverlap = minOverlap_globe,
                 protein_coding_transcripts = protein_coding_transcripts_globe)
    }))
  } else if (ex_type %in% c("A5SS", "A3SS")) {
    results <- unlist(parallel::mclapply(seq(1, nrow(redExon_globe), by = 2), mc.cores = cores, function(i) {
      ASmatcher(i, redExon = redExon_globe, minOverlap = minOverlap_globe,
                 gtf_transcripts = gtf_transcripts_globe,
                 gtf_exons = gtf_exons_globe,
                 protein_coding_transcripts = protein_coding_transcripts_globe,
                 transcript_starts = transcript_starts_globe, whichType = ex_type)
    }))
  } else if (ex_type %in% c("MXE")) {
    results <- unlist(parallel::mclapply(seq(1, nrow(redExon_globe), by = 2), mc.cores = cores, function(i) {
      MXmatcher(i, redExon = redExon_globe, minOverlap = minOverlap_globe,
                gtf_transcripts = gtf_transcripts_globe,
                gtf_exons = gtf_exons_globe,
                protein_coding_transcripts = protein_coding_transcripts_globe,
                transcript_starts = transcript_starts_globe)
    }))
  } else if (ex_type %in% c("SE")) {
    results <- unlist(parallel::mclapply(seq(1, nrow(redExon_globe), by = 2), mc.cores = cores, function(i) {
      SEmatcher(i, redExon = redExon_globe, minOverlap = minOverlap_globe,
                gtf_transcripts = gtf_transcripts_globe,
                gtf_exons = gtf_exons_globe,
                protein_coding_transcripts = protein_coding_transcripts_globe,
                transcript_starts = transcript_starts_globe)
    }))
  } else if (ex_type %in% c("RI")) {
    results <- unlist(parallel::mclapply(seq(1, nrow(redExon_globe), by = 2), mc.cores = cores, function(i) {
      RImatcher(i, redExon = redExon_globe, minOverlap = minOverlap_globe,
                gtf_transcripts = gtf_transcripts_globe,
                gtf_exons = gtf_exons_globe,
                protein_coding_transcripts = protein_coding_transcripts_globe,
                transcript_starts = transcript_starts_globe)
    }))
  }
  return(results)
}

#' specific matcher for HIT
#' @param i index of redExon
#' @param redExon dataframe generated initially in getForeground
#' @param gtf_filtered gtf filtered for gene
#' @param minOverlap minimum overlap to count as same exon
#' @param protein_coding_transcripts protein coding transcripts
#' @return matched transcript rownumber
#' @keywords internal
HITmatcher <- function(i, redExon, gtf_filtered, minOverlap,
                       protein_coding_transcripts) {

  # Initiate geneR, start, stop
  geneR <- redExon$geneR[i]
  start <- redExon$start[i]
  stop <- redExon$stop[i]

  # Further filter for computational efficiency to current gene
  gtf_min <- gtf_filtered[gtf_filtered$geneID == geneR,]

  # Skip if no relevant gtf entries found
  if (nrow(gtf_min) == 0) {
    return(0)
  }

  # Calculate Jaccard index for each gtf entry
  jacc <- calculate_jaccard(start, stop, gtf_min$start, gtf_min$stop)
  gtf_min$length_jacc <- jacc[[2]]
  gtf_min$jaccard <- jacc[[1]]

  # Filter based on minOverlap
  gtf_min <- gtf_min[gtf_min$jaccard >= minOverlap,]

  # Skip again if no relevant gtf entries found with jaccard index over the minOverlap
  if (nrow(gtf_min) == 0) {
    return(0)
  }

  # Further refine to only include those that are also in the protein coding transcripts, if necessary
  if (sum(gtf_min$transcriptID %in% protein_coding_transcripts) > 0) {
    gtf_min <- gtf_min[gtf_min$transcriptID %in% protein_coding_transcripts,]
  }

  # Order by Jaccard index to find best match
  gtf_min <- gtf_min %>% dplyr::arrange(dplyr::desc(jaccard))

  ## If only one match or one best match
  if (nrow(gtf_min) == 1 | gtf_min$jaccard[1] > gtf_min$jaccard[2]) {
    return(gtf_min$rownum[1])
  }
  ## If multiple best, use length of matched exon
  else if (gtf_min$jaccard[1] == gtf_min$jaccard[2]) {

    c_gtf <- gtf_min[gtf_min$jaccard == max(gtf_min$jaccard),] %>% dplyr::arrange(desc(.data$length_jacc))

    return(c_gtf$rownum[1])
  } else {return(0)}
}

#' specific matcher for SE
#' @param i index of redExon
#' @param redExon dataframe generated initially in getForeground
#' @param gtf_filtered gtf filtered for gene
#' @param gtf_transcripts transcripts gtf dataframe
#' @param gtf_exons exons gtf dataframe
#' @param minOverlap minimum overlap to count as same exon
#' @param protein_coding_transcripts protein coding transcripts
#' @param below_thresh threshold to make skipped form match less than
#' @return matched transcript rownumber
#' @keywords internal
SEmatcher <- function(i, below_thresh = 0, redExon, minOverlap = .05,
                      gtf_transcripts,
                      gtf_exons,
                      protein_coding_transcripts,
                      transcript_starts) {

  # Initiate geneR, start, stop
  geneR <- redExon$geneR[i]
  start <- redExon$start[i]
  stop <- redExon$stop[i]

  # Pre-process redExon$add_inf to avoid repeated calculations
  split_add_inf <- strsplit(redExon$add_inf[i], split = ";")[[1]][2:3]
  up_down <- as.numeric(unlist(strsplit(split_add_inf, split = "-")))

  # Assign values directly without unnecessary lapply and unlist
  up_down_start <- up_down[c(TRUE, FALSE)]
  up_down_stop <- up_down[c(FALSE, TRUE)]

  # Filter for computational efficiency at the beginning to minimize dataset size
  gtf_filtered <- subset(gtf_exons, geneID == geneR & chr == redExon$chr[i])
  if (nrow(gtf_filtered) <= 1) return(c(0, 0)) # Skip if less than 2 relevant gtf entries found

  # Calculate Jaccard index for each gtf entry once instead of three times
  jaccard_indices <- lapply(list(c(start, stop), c(up_down_start[1], up_down_stop[1]), c(up_down_start[2], up_down_stop[2])),
                            function(range) calculate_jaccard(range[1], range[2], gtf_filtered$start, gtf_filtered$stop))


  # Compute intersections more efficiently
  possible_inclusion_transcripts <- Reduce(intersect, lapply(jaccard_indices, function(jacc)
    unique(gtf_filtered$transcriptID[jacc$jaccard_index > minOverlap])))

  # Function to determine if all Jaccard indices for a transcript are below a threshold
  is_below_threshold <- function(transcript_id, threshold = below_thresh) {
    jacc_indices <- jaccard_indices[[1]]$jaccard_index[gtf_filtered$transcriptID == transcript_id]
    all(jacc_indices <= threshold)
  }
  # Compute intersections more efficiently
  possible_exclusion_transcripts <- Reduce(intersect, c(lapply(jaccard_indices[2:3], function(jacc) unique(gtf_filtered$transcriptID[jacc$jaccard_index > minOverlap])),
                                                        list(gtf_filtered$transcriptID[unlist(lapply(gtf_filtered$transcriptID, function(x)
                                                          is_below_threshold(x, threshold = below_thresh)))])))


  # Further refine to only include those that are also in the protein coding transcripts, if necessary
  pc_exclusion <- possible_exclusion_transcripts[possible_exclusion_transcripts %in% protein_coding_transcripts]

  # Apply protein coding filter again
  pc_inclusion <- possible_inclusion_transcripts[possible_inclusion_transcripts %in% protein_coding_transcripts]

  inclusion <- ifelse(length(pc_inclusion) == 0, list(possible_inclusion_transcripts), list(pc_inclusion))[[1]]
  exclusion <- ifelse(length(pc_exclusion) == 0, list(possible_exclusion_transcripts), list(pc_exclusion))[[1]]
  exclusion <- exclusion[!(exclusion %in% inclusion)]
  exclusion_lengths <- unlist(lapply(exclusion, function(tid) {
    abs(gtf_transcripts$start[gtf_transcripts$transcriptID == tid]-
          gtf_transcripts$stop[gtf_transcripts$transcriptID == tid])
  }))
  exclusion_rownum <- gtf_transcripts$rownum[gtf_transcripts$transcriptID == exclusion[which.max(exclusion_lengths)] & gtf_transcripts$chr == redExon$chr[i]]
  # Skip if no exclusion
  if (length(exclusion) == 0) {return(c(0, 0))}

  # Apply jaccard index and classification filter directly
  gtf_filtered$jaccard <- jaccard_indices[[1]]$jaccard_index
  gtf_filtered$length_jacc <- jaccard_indices[[1]]$length_jacc

  gtf_filtered <- subset(gtf_filtered, jaccard > minOverlap & classification == "internal" & transcriptID %in% inclusion)

  if (nrow(gtf_filtered) == 0) return(c(0, 0)) # Skip if no entries found

  # Order and find the best match more efficiently
  gtf_filtered <- gtf_filtered[order(-gtf_filtered$jaccard),]


  # Directly return results based on condition
  if (nrow(gtf_filtered) == 1 || gtf_filtered$jaccard[1] > gtf_filtered$jaccard[2]) {
    return(c(gtf_filtered$rownum[1], exclusion_rownum)) # Complete the logic for finding the correct rownum
  } else if (gtf_filtered$jaccard[1] == gtf_filtered$jaccard[2]) {

    gtf_min <- gtf_filtered[gtf_filtered$jaccard == max(gtf_filtered$jaccard),]

    # Pre-compute the start positions for transcripts in exclusion
    exclusion_start <- gtf_transcripts$start[gtf_transcripts$rownum == exclusion_rownum]


    # Vectorize the calculation of start distances
    start_distances <- abs(transcript_starts[gtf_min$transcriptID] - exclusion_start)

    # Ensure start_distances is a vector if it's not due to subsetting a single element
    start_distances <- as.numeric(start_distances)

    c_gtf <- gtf_min[start_distances == min(start_distances),] %>% dplyr::arrange(desc(.data$length_jacc))

    inclusion_rownum <- c_gtf$rownum[1]

    return(c(inclusion_rownum, exclusion_rownum))
  } else {
    return(c(0, 0))
  }
}

#' specific matcher for MXE
#' @return matched transcript rownumber
#' @param i index of redExon
#' @param redExon dataframe generated initially in getForeground
#' @param gtf_transcripts transcripts gtf dataframe
#' @param gtf_exons exons gtf dataframe
#' @param minOverlap minimum overlap to count as same exon
#' @param protein_coding_transcripts protein coding transcripts
#' @param below_thresh threshold to make skipped form match less than
#' @param transcript_starts start location of transcript
#' @keywords internal
MXmatcher <- function(i, below_thresh = 0, redExon, minOverlap = .05,
                      gtf_transcripts,
                      gtf_exons,
                      protein_coding_transcripts,
                      transcript_starts) {

  # Initiate geneR, start, stop
  geneR <- redExon$geneR[i]
  start <- redExon$start[i]
  stop <- redExon$stop[i]

  # Pre-process redExon$add_inf to avoid repeated calculations
  split_add_inf <- strsplit(strsplit(redExon$add_inf[i], split = ";")[[1]][1], split = ":")[[1]][2]
  up_down <- as.numeric(unlist(strsplit(split_add_inf, split = "-")))

  # Assign values directly without unnecessary lapply and unlist
  excl_start <- up_down[1]
  excl_stop <- up_down[2]

  # Filter for computational efficiency at the beginning to minimize dataset size
  gtf_filtered <- subset(gtf_exons, classification == "internal" & geneID == geneR & chr == redExon$chr[i])
  if (nrow(gtf_filtered) <= 1) return(c(0, 0)) # Skip if less than 2 relevant gtf entries found

  # Calculate Jaccard index for inclusion and finds intersect with exclusion of second exon
  inclusion_indices <- calculate_jaccard(start, stop, gtf_filtered$start, gtf_filtered$stop)

  exclusion_indices <- calculate_jaccard(excl_start, excl_stop, gtf_filtered$start, gtf_filtered$stop)

  # Function to determine if all Jaccard indices for a transcript are below a threshold
  is_below_threshold <- function(transcript_id, jacc, threshold =  below_thresh) {
    jacc_indices <- jacc$jaccard_index[gtf_filtered$transcriptID == transcript_id]
    all(jacc_indices <= threshold)
  }

  possible_inclusion_transcripts <- intersect(unique(gtf_filtered$transcriptID[inclusion_indices$jaccard_index > minOverlap]),
                                              gtf_filtered$transcriptID[unlist(lapply(gtf_filtered$transcriptID, function(tid)
                                                is_below_threshold(tid, exclusion_indices)))])

  possible_exclusion_transcripts <- intersect(unique(gtf_filtered$transcriptID[exclusion_indices$jaccard_index > minOverlap]),
                                              gtf_filtered$transcriptID[unlist(lapply(gtf_filtered$transcriptID, function(tid)
                                                is_below_threshold(tid, inclusion_indices)))])


  # Further refine to only include those that are also in the protein coding transcripts, if necessary
  pc_exclusion <- possible_exclusion_transcripts[possible_exclusion_transcripts %in% protein_coding_transcripts]

  # Apply protein coding filter again
  pc_inclusion <- possible_inclusion_transcripts[possible_inclusion_transcripts %in% protein_coding_transcripts]

  inclusion <- ifelse(length(pc_inclusion) == 0, list(possible_inclusion_transcripts), list(pc_inclusion))[[1]]
  exclusion <- ifelse(length(pc_exclusion) == 0, list(possible_exclusion_transcripts), list(pc_exclusion))[[1]]
  exclusion <- exclusion[!(exclusion %in% inclusion)]
  exclusion_lengths <- unlist(lapply(exclusion, function(tid) {
    abs(gtf_transcripts$start[gtf_transcripts$transcriptID == tid]-
          gtf_transcripts$stop[gtf_transcripts$transcriptID == tid])
  }))
  exclusion_rownum <- gtf_transcripts$rownum[gtf_transcripts$transcriptID == exclusion[which.max(exclusion_lengths)] & gtf_transcripts$chr == redExon$chr[i]]
  # Skip if no exclusion
  if (length(exclusion) == 0 | length(inclusion) == 0) {return(c(0, 0))}

  exclusion_rownum <- getMXE_internal(gtf_filtered, exclusion_indices, exclusion, "excl", minOverlap)
  inclusion_rownum <- getMXE_internal(gtf_filtered, inclusion_indices, inclusion, "incl", minOverlap)
  if (sum(exclusion_rownum) == 0 | sum(inclusion_rownum) == 0) {
    return(c(0, 0))
  }
  return(c(inclusion_rownum, exclusion_rownum))
}


#' helper for MXEmatcher
#' @return rownums for the clusion
#' @keywords internal
getMXE_internal <- function(g, indices, clusion, strVar, minOverlap) {
  g$jaccard <- indices$jaccard_index
  g$length_jacc <- indices$length_jacc
  g <- subset(g, jaccard > minOverlap & classification == "internal" & transcriptID %in% clusion)
  if (nrow(g) == 0) return(c(0, 0)) # Skip if no entries found
  g <- g[order(-g$jaccard),]


  if (nrow(g) == 1 || g$jaccard[1] > g$jaccard[2]) {
    return(c(g$rownum[1])) # Complete the logic for finding the correct rownum
  } else if (g$jaccard[1] == g$jaccard[2]) {

    gtf_min <- g[g$jaccard == max(g$jaccard),]
    clusion <- clusion[g$jaccard == max(g$jaccard)]
    # Pre-compute the start positions for transcripts in exclusion
    if (strVar == "inc") {
      exclusion_start <- gtf_transcripts$start[gtf_transcripts$rownum == exclusion_rownum]


      # Vectorize the calculation of start distances
      start_distances <- abs(transcript_starts[clusion] - exclusion_start)

      # Ensure start_distances is a vector if it's not due to subsetting a single element
      start_distances <- as.numeric(start_distances)
      c_gtf <- gtf_min[start_distances == min(start_distances),] %>% dplyr::arrange(desc(.data$length_jacc))

    } else {
      c_gtf <- gtf_min %>% dplyr::arrange(desc(.data$length_jacc))

    }

    clusion_rownum <- c_gtf$rownum[1]

    return(clusion_rownum)
  } else {
    return(c(0, 0))
  }
}

#' specific matcher for A5SS/A3SS
#' @return matched transcript rownumber
#' @keywords internal
ASmatcher <- function(i, below_thresh = .2, redExon, minOverlap = .05,
                      gtf_transcripts,
                      gtf_exons,
                      protein_coding_transcripts,
                      transcript_starts, whichType) {

  # Function to calculate Jaccard-like index more efficiently
  calculate_jaccard_like <- function(start1, stop1, start2, stop2) {
    intersection_length <- pmax(0, pmin(stop1, stop2) - pmax(start1, start2) + 1)
    union_length <- (stop1 - start1 + 1) + (stop2 - start2 + 1) - intersection_length
    jaccard_index <- intersection_length / (abs(start1-stop1)+1)
    length_jacc <- intersection_length/(abs(start2-stop2)+1)
    return(list(jaccard_index = jaccard_index,
                length_jacc = length_jacc))
  }

  # Initiate geneR, start, stop
  geneR <- redExon$geneR[i]
  start <- redExon$start[i]
  stop <- redExon$stop[i]

  # Pre-process redExon$add_inf to avoid repeated calculations
  split_add_inf <- strsplit(strsplit(redExon$add_inf[i], split = ";")[[1]][1], split = ":")[[1]][2]
  up_down <- as.numeric(unlist(strsplit(split_add_inf, split = "-")))

  # Extract set difference between two sets of coords
  idSS <- function(start1, stop1, start2, stop2) {
    if (stop1-stop2 == 0) {
      return(sort(c(start1, start2), decreasing=FALSE)-c(0, 1))
    } else {
      return(sort(c(stop1, stop2), decreasing =FALSE)+c(1, 0))}}


  # Assign values directly without unnecessary lapply and unlist
  exclusion_location <- idSS(start, stop, up_down[1], up_down[2])

  # Filter for computational efficiency at the beginning to minimize dataset size
  gtf_filtered <- subset(gtf_exons, classification %in% c("first", "last", "internal") & geneID == geneR & chr == redExon$chr[i])
  if (nrow(gtf_filtered) <= 1) return(c(0, 0))
  if (whichType == "A5SS" & unique(gtf_filtered$strand) == "+") {
    gtf_filtered <- gtf_filtered[gtf_filtered$classification %in% c("first", "internal"),]
  } else if (whichType == "A5SS" & unique(gtf_filtered$strand) == "-") {
    gtf_filtered <- gtf_filtered[gtf_filtered$classification %in% c("last", "internal"),]
  } else if (whichType == "A3SS" & unique(gtf_filtered$strand) == "+") {
    gtf_filtered <- gtf_filtered[gtf_filtered$classification %in% c("last", "internal"),]
  } else if (whichType == "A3SS" & unique(gtf_filtered$strand) == "-") {
    gtf_filtered <- gtf_filtered[gtf_filtered$classification %in% c("first", "internal"),]
  }

  if (nrow(gtf_filtered) <= 1) return(c(0, 0)) # Skip if less than 2 relevant gtf entries found

  # Calculate Jaccard index for inclusion and finds intersect with exclusion of second exon
  inclusion_indices <- calculate_jaccard(start, stop, gtf_filtered$start, gtf_filtered$stop)

  exclusion_extension_indices <- calculate_jaccard_like(exclusion_location[1], exclusion_location[2], gtf_filtered$start, gtf_filtered$stop)

  exclusion_overlap_indices <- calculate_jaccard(up_down[1], up_down[2], gtf_filtered$start, gtf_filtered$stop)

  # Function to determine if all Jaccard indices for a transcript are below a threshold
  is_below_threshold <- function(transcript_id, jacc, threshold =  below_thresh) {
    jacc_indices <- jacc$jaccard_index[gtf_filtered$transcriptID == transcript_id]
    all(jacc_indices <= threshold)
  }

  possible_inclusion_transcripts <- unique(gtf_filtered$transcriptID[inclusion_indices$jaccard_index > minOverlap &
                                                                       exclusion_extension_indices$jaccard_index > minOverlap])


  possible_exclusion_transcripts <- unique(gtf_filtered$transcriptID[exclusion_overlap_indices$jaccard_index > minOverlap &
                                                                       exclusion_extension_indices$jaccard_index < minOverlap])



  # Further refine to only include those that are also in the protein coding transcripts, if necessary
  pc_exclusion <- possible_exclusion_transcripts[possible_exclusion_transcripts %in% protein_coding_transcripts]

  # Apply protein coding filter again
  pc_inclusion <- possible_inclusion_transcripts[possible_inclusion_transcripts %in% protein_coding_transcripts]

  inclusion <- ifelse(length(pc_inclusion) == 0, list(possible_inclusion_transcripts), list(pc_inclusion))[[1]]
  exclusion <- ifelse(length(pc_exclusion) == 0, list(possible_exclusion_transcripts), list(pc_exclusion))[[1]]
  exclusion <- exclusion[!(exclusion %in% inclusion)]
  exclusion_lengths <- unlist(lapply(exclusion, function(tid) {
    abs(gtf_transcripts$start[gtf_transcripts$transcriptID == tid]-
          gtf_transcripts$stop[gtf_transcripts$transcriptID == tid])
  }))

  # Skip if no exclusion
  if (length(exclusion) == 0 | length(inclusion) == 0) {return(c(0, 0))}

  exclusion_rownum <- getAS_internal(gtf_filtered, exclusion_overlap_indices, exclusion, "excl", minOverlap)
  inclusion_rownum <- getAS_internal(gtf_filtered, inclusion_indices, inclusion, "incl", minOverlap)

  if (sum(exclusion_rownum) == 0 | sum(inclusion_rownum) == 0) {
    return(c(0, 0))
  }
  return(c(inclusion_rownum, exclusion_rownum))
}


#' helper for MXEmatcher
#' @return rownums for the clusion
#' @keywords internal
getAS_internal <- function(g, indices, clusion, strVar, minOverlap) {
  g$jaccard <- indices$jaccard_index
  g$length_jacc <- indices$length_jacc
  g <- subset(g, jaccard > minOverlap & transcriptID %in% clusion)
  if (nrow(g) == 0) return(c(0, 0)) # Skip if no entries found
  g <- g[order(-g$jaccard),]


  if (nrow(g) == 1 || g$jaccard[1] > g$jaccard[2]) {
    return(c(g$rownum[1])) # Complete the logic for finding the correct rownum
  } else if (g$jaccard[1] == g$jaccard[2]) {

    gtf_min <- g[g$jaccard == max(g$jaccard),]
    clusion <- clusion[g$jaccard == max(g$jaccard)]
    # Pre-compute the start positions for transcripts in exclusion
    if (strVar == "inc") {
      exclusion_start <- gtf_transcripts$start[gtf_transcripts$rownum == exclusion_rownum]


      # Vectorize the calculation of start distances
      start_distances <- abs(transcript_starts[clusion] - exclusion_start)

      # Ensure start_distances is a vector if it's not due to subsetting a single element
      start_distances <- as.numeric(start_distances)
      c_gtf <- gtf_min[start_distances == min(start_distances),] %>% dplyr::arrange(desc(.data$length_jacc))

    } else {
      c_gtf <- gtf_min %>% dplyr::arrange(desc(.data$length_jacc))

    }

    clusion_rownum <- c_gtf$rownum[1]

    return(clusion_rownum)
  } else {
    return(c(0, 0))
  }
}

#' specific matcher for RI
#' @return matched transcript rownumber
#' @keywords internal
RImatcher <- function(i, below_thresh = .2, redExon, minOverlap = .05,
                      gtf_transcripts,
                      gtf_exons,
                      protein_coding_transcripts,
                      transcript_starts) {

  # Function to calculate Jaccard-like index more efficiently
  calculate_jaccard_like <- function(start1, stop1, start2, stop2) {
    intersection_length <- pmax(0, pmin(stop1, stop2) - pmax(start1, start2) + 1)
    union_length <- (stop1 - start1 + 1) + (stop2 - start2 + 1) - intersection_length
    jaccard_index <- intersection_length / (abs(start1-stop1)+1)
    length_jacc <- intersection_length/(abs(start2-stop2)+1)
    return(list(jaccard_index = jaccard_index,
                length_jacc = length_jacc))
  }

  # Initiate geneR, start, stop
  geneR <- redExon$geneR[i]
  start <- redExon$start[i]
  stop <- redExon$stop[i]

  # Pre-process redExon$add_inf to avoid repeated calculations
  split_add_inf <- strsplit(strsplit(redExon$add_inf[i], split = ";")[[1]][2:3], split = "-")
  up_down <- list(c(split_add_inf[[1]][1], split_add_inf[[2]][1]),
                  c(split_add_inf[[1]][2], split_add_inf[[2]][2]))

  # Assign values directly without unnecessary lapply and unlist
  up_down_start <- as.numeric(up_down[c(TRUE, FALSE)][[1]])
  up_down_stop <- as.numeric(up_down[c(FALSE, TRUE)][[1]])

  # intron retention location
  ir_exclusion_location <- c(up_down_stop[which.min(up_down_stop)]+1, up_down_start[which.max(up_down_start)]-1)


  # Filter for computational efficiency at the beginning to minimize dataset size
  gtf_filtered <- subset(gtf_exons, geneID == geneR & chr == redExon$chr[i])
  if (nrow(gtf_filtered) <= 1) return(c(0, 0)) # Skip if less than 2 relevant gtf entries found


  # Calculate Jaccard index for each gtf entry once instead of three times
  jaccard_indices <- lapply(list(c(start, stop), c(up_down_start[1], up_down_stop[1]),
                                 c(up_down_start[2], up_down_stop[2]), c(ir_exclusion_location[1], ir_exclusion_location[2])),
                            function(range) calculate_jaccard(range[1], range[2], gtf_filtered$start, gtf_filtered$stop))


  # Compute intersections more efficiently
  possible_inclusion_transcripts <- unique(gtf_filtered$transcriptID[jaccard_indices[[1]]$jaccard_index > minOverlap])

  # Function to determine if all Jaccard indices for a transcript are below a threshold
  is_below_threshold_ri <- function(transcript_id, threshold = below_thresh) {
    jacc_indices <- jaccard_indices[[4]]$jaccard_index[gtf_filtered$transcriptID == transcript_id]
    all(jacc_indices < threshold)
  }

  # Compute intersections more efficiently
  possible_exclusion_exons <- lapply(jaccard_indices[2:3], function(jacc) {
    ji <- jacc$jaccard_index[jacc$jaccard_index > minOverlap]
    gff <- gtf_filtered[jacc$jaccard_index > minOverlap,] %>% arrange(desc(ji))
    if (nrow(gff) == 0) {return(c(0, 0))}
    gff$exonID[1]
  })
  possible_exclusion_exons[[1]] <- unique(possible_exclusion_exons[[1]])
  possible_exclusion_exons[[2]] <- unique(possible_exclusion_exons[[2]])[!(unique(possible_exclusion_exons[[2]]) %in%
                                                                             possible_exclusion_exons[[1]])]
  if (sum(lengths(possible_exclusion_exons) > 0) < 2) {return(c(0, 0))}
  possible_exclusion_transcripts <- unique(gtf_filtered$transcriptID[unlist(lapply(gtf_filtered$transcriptID,
                                                                            function(x) is_below_threshold_ri(x, threshold = below_thresh)))])

  pet_split <- split(gtf_filtered[gtf_filtered$transcriptID %in% possible_exclusion_transcripts,],
                     gtf_filtered$transcriptID[gtf_filtered$transcriptID %in% possible_exclusion_transcripts])

  possible_exclusion_transcripts <- names(pet_split)[unlist(lapply(pet_split, function(x) {
    sum(sum(x$exonID %in% possible_exclusion_exons[[1]]),
        sum(x$exonID %in% possible_exclusion_exons[[2]])) == 2
  }))]


  # Further refine to only include those that are also in the protein coding transcripts, if necessary
  pc_exclusion <- possible_exclusion_transcripts[possible_exclusion_transcripts %in% protein_coding_transcripts]

  # Apply protein coding filter again

  inclusion <- possible_inclusion_transcripts
  exclusion <- ifelse(length(pc_exclusion) == 0, list(possible_exclusion_transcripts), list(pc_exclusion))[[1]]
  exclusion <- exclusion[!(exclusion %in% inclusion)]
  exclusion_lengths <- unlist(lapply(exclusion, function(tid) {
    abs(gtf_transcripts$start[gtf_transcripts$transcriptID == tid]-
          gtf_transcripts$stop[gtf_transcripts$transcriptID == tid])
  }))

  # Skip if no exclusion
  if (length(exclusion) == 0 | length(inclusion) == 0) {return(c(0, 0))}

  inclusion_rownum <- getRI_internal(gtf_filtered, jaccard_indices[[1]], inclusion, "incl", minOverlap)

  if (unique(gtf_filtered$strand) == "+") {
    downstream_exclusion_exon <- unlist(possible_exclusion_exons)[which.max(gtf_filtered$start[match(unlist(possible_exclusion_exons),
                                                                                                     gtf_filtered$exonID) ])]
    exclusion_rownum <- gtf_filtered$rownum[gtf_filtered$transcriptID == exclusion[which.max(exclusion_lengths)] &
                                              gtf_filtered$chr == redExon$chr[i] &
                                              gtf_filtered$exonID == downstream_exclusion_exon]
  } else {
    downstream_exclusion_exon <- unlist(possible_exclusion_exons)[which.min(gtf_filtered$start[match(unlist(possible_exclusion_exons),
                                                                                                     gtf_filtered$exonID) ])]
    exclusion_rownum <- gtf_filtered$rownum[gtf_filtered$transcriptID == exclusion[which.max(exclusion_lengths)] &
                                              gtf_filtered$chr == redExon$chr[i] &
                                              gtf_filtered$exonID == downstream_exclusion_exon]
  }
  if (sum(exclusion_rownum) == 0 | sum(inclusion_rownum) == 0) {
    return(c(0, 0))
  }
  return(c(inclusion_rownum, exclusion_rownum))
}

#' helper for RImatcher
#' @return rownums for the clusion
#' @keywords internal
getRI_internal <- function(g, indices, clusion, strVar, minOverlap) {
  g$jaccard <- indices$jaccard_index
  g$length_jacc <- indices$length_jacc
  g <- subset(g, jaccard > minOverlap & transcriptID %in% clusion)
  if (nrow(g) == 0) return(c(0, 0)) # Skip if no entries found
  g <- g[order(-g$jaccard),]


  if (nrow(g) == 1 || g$jaccard[1] > g$jaccard[2]) {
    return(c(g$rownum[1])) # Complete the logic for finding the correct rownum
  } else if (g$jaccard[1] == g$jaccard[2]) {

    gtf_min <- g[g$jaccard == max(g$jaccard),]
    clusion <- clusion[g$jaccard == max(g$jaccard)]
    # Pre-compute the start positions for transcripts in exclusion
    if (strVar == "inc") {
      exclusion_start <- gtf_transcripts$start[gtf_transcripts$rownum == exclusion_rownum]


      # Vectorize the calculation of start distances
      start_distances <- abs(transcript_starts[clusion] - exclusion_start)

      # Ensure start_distances is a vector if it's not due to subsetting a single element
      start_distances <- as.numeric(start_distances)
      c_gtf <- gtf_min[start_distances == min(start_distances),] %>% dplyr::arrange(desc(.data$length_jacc))

    } else {
      c_gtf <- gtf_min %>% dplyr::arrange(desc(.data$length_jacc))

    }

    clusion_rownum <- c_gtf$rownum[1]

    return(clusion_rownum)
  } else {
    return(c(0, 0))
  }
}
