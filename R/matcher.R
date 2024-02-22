matcher <- function(ex_type, background = F, cores) {

  gtf_transcripts <- gtf[gtf$classification == 'transcript',]
  gtf_exons <- gtf[!(gtf$classification %in% c('gene', 'transcript')),]

  # Filter protein_coding_transcripts once and for all
  protein_coding_transcripts <- unique(gtf$transcriptID[gtf$tpc == "protein_coding" & !is.na(gtf$tpc)])

  # Pre-compute a lookup for start positions of all transcripts in gtf
  transcript_starts <- setNames(gtf$start[gtf$classification == 'transcript'], gtf$transcriptID[gtf$classification == 'transcript'])

  if (ex_type %in% c("AFE", "ALE") | background) {
    results <- unlist(parallel::mclapply(1:nrow(redExon), mc.cores = cores, function(i) {
      HITmatcher(i)
    }))
  } else if (ex_type %in% c("A5SS", "A3SS")) {
    results <- unlist(parallel::mclapply(1:nrow(redExon), mc.cores = cores, function(i) {
      ASmatcher(i)
    }))
  } else if (ex_type %in% c("MXE")) {
    results <- unlist(parallel::mclapply(1:nrow(redExon), mc.cores = cores, function(i) {
      MXmatcher(i)
    }))
  } else if (ex_type %in% c("SE")) {
    results <- unlist(parallel::mclapply(1:nrow(redExon), mc.cores = cores, function(i) {
      SEmatcher(i)
    }))
  }
  return(results)
}
HITmatcher <- function(i) {
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

    # Further refine to only include those that are also in the protein coding transcripts, if necessary
    pc_gtf_min <- gtf_min[gtf_min$transcriptID %in% protein_coding_transcripts,]
    if (nrow(pc_gtf_min) == 0) {
      c_gtf <- gtf_min[gtf_min$jaccard == max(gtf_min$jaccard),] %>% dplyr::arrange(desc(length_jacc))
    } else {
      c_gtf <- pc_gtf_min[pc_gtf_min$jaccard == max(pc_gtf_min$jaccard),] %>% dplyr::arrange(desc(length_jacc))
    }
    return(c_gtf$rownum[1])
  } else {return(0)}
}

SEmatcher <- function(i, below_thresh = .2) {

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
    all(jacc_indices < threshold)
  }
  # Compute intersections more efficiently
  possible_exclusion_transcripts <- Reduce(intersect, c(lapply(jaccard_indices[2:3], function(jacc) unique(gtf_filtered$transcriptID[jacc$jaccard_index > minOverlap])), list(gtf_filtered$transcriptID[sapply(gtf_filtered$transcriptID, function(x) is_below_threshold(x, threshold = below_thresh))])))


  # Further refine to only include those that are also in the protein coding transcripts, if necessary
  pc_exclusion <- possible_exclusion_transcripts[possible_exclusion_transcripts %in% protein_coding_transcripts]

  # Apply protein coding filter again
  pc_inclusion <- possible_inclusion_transcripts[possible_inclusion_transcripts %in% protein_coding_transcripts]

  inclusion <- ifelse(length(pc_inclusion) == 0, list(possible_inclusion_transcripts), list(pc_inclusion))[[1]]
  exclusion <- ifelse(length(pc_exclusion) == 0, list(possible_exclusion_transcripts), list(pc_exclusion))[[1]]
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

    c_gtf <- gtf_min[start_distances == min(start_distances),] %>% dplyr::arrange(desc(length_jacc))

    inclusion_rownum <- c_gtf$rownum[1]

    return(c(inclusion_rownum, exclusion_rownum))
  } else {
    return(c(0, 0))
  }
}

MXmatcher <- function(i, below_thresh = .2) {

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
    all(jacc_indices < threshold)
  }

  possible_inclusion_transcripts <- intersect(unique(gtf_filtered$transcriptID[inclusion_indices$jaccard_index > minOverlap]),
                                              gtf_filtered$transcriptID[sapply(gtf_filtered$transcriptID, function(tid)
                                                is_below_threshold(tid, exclusion_indices))])

  possible_exclusion_transcripts <- intersect(unique(gtf_filtered$transcriptID[exclusion_indices$jaccard_index > minOverlap]),
                                              gtf_filtered$transcriptID[sapply(gtf_filtered$transcriptID, function(tid)
                                                is_below_threshold(tid, inclusion_indices))])


  # Further refine to only include those that are also in the protein coding transcripts, if necessary
  pc_exclusion <- possible_exclusion_transcripts[possible_exclusion_transcripts %in% protein_coding_transcripts]

  # Apply protein coding filter again
  pc_inclusion <- possible_inclusion_transcripts[possible_inclusion_transcripts %in% protein_coding_transcripts]

  inclusion <- ifelse(length(pc_inclusion) == 0, list(possible_inclusion_transcripts), list(pc_inclusion))[[1]]
  exclusion <- ifelse(length(pc_exclusion) == 0, list(possible_exclusion_transcripts), list(pc_exclusion))[[1]]
  exclusion_lengths <- unlist(lapply(exclusion, function(tid) {
    abs(gtf_transcripts$start[gtf_transcripts$transcriptID == tid]-
          gtf_transcripts$stop[gtf_transcripts$transcriptID == tid])
  }))
  exclusion_rownum <- gtf_transcripts$rownum[gtf_transcripts$transcriptID == exclusion[which.max(exclusion_lengths)] & gtf_transcripts$chr == redExon$chr[i]]
  # Skip if no exclusion
  if (length(exclusion) == 0) {return(c(0, 0))}

  # Apply jaccard index and classification filter directly
  gtf_filtered$jaccard <- inclusion_indices$jaccard_index
  gtf_filtered$length_jacc <- inclusion_indices$length_jacc

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
    start_distances <- abs(transcript_starts[inclusion] - exclusion_start)

    # Ensure start_distances is a vector if it's not due to subsetting a single element
    start_distances <- as.numeric(start_distances)

    c_gtf <- gtf_min[start_distances == min(start_distances),] %>% dplyr::arrange(desc(length_jacc))

    inclusion_rownum <- c_gtf$rownum[1]

    return(c(inclusion_rownum, exclusion_rownum))
  } else {
    return(c(0, 0))
  }
}

ASmatcher <- function(i, below_thresh = .2) {

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
      return(sort(c(start1, start2), decreasing=F)-c(0, 1))
    } else {
      return(sort(c(stop1, stop2), decreasing = F)+c(1, 0))}}


  # Assign values directly without unnecessary lapply and unlist
  exclusion_location <- idSS(start, stop, up_down[1], up_down[2])

  # Filter for computational efficiency at the beginning to minimize dataset size
  gtf_filtered <- subset(gtf_exons, classification %in% c("first", "last", "internal") & geneID == geneR & chr == redExon$chr[i])
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
  exclusion_lengths <- unlist(lapply(exclusion, function(tid) {
    abs(gtf_transcripts$start[gtf_transcripts$transcriptID == tid]-
          gtf_transcripts$stop[gtf_transcripts$transcriptID == tid])
  }))
  exclusion_rownum <- gtf_transcripts$rownum[gtf_transcripts$transcriptID == exclusion[which.max(exclusion_lengths)] & gtf_transcripts$chr == redExon$chr[i]]
  # Skip if no exclusion
  if (length(exclusion) == 0) {return(c(0, 0))}

  # Apply jaccard index and classification filter directly
  gtf_filtered$jaccard <- inclusion_indices$jaccard_index
  gtf_filtered$length_jacc <- inclusion_indices$length_jacc

  gtf_filtered <- subset(gtf_filtered, jaccard > minOverlap & transcriptID %in% inclusion)

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

    c_gtf <- gtf_min[start_distances == min(start_distances),] %>% dplyr::arrange(desc(length_jacc))

    inclusion_rownum <- c_gtf$rownum[1]

    return(c(inclusion_rownum, exclusion_rownum))
  } else {
    return(c(0, 0))
  }
}
