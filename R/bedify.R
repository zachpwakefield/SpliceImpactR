bedifyForeground <- function(matched, outname, cores) {

  # Function to process each transcriptID and construct BED rows
  processTranscriptID_fg <- function(i) {
    # Efficiently filter for relevant rows
    bed <- gtf[gtf$transcriptID == matched$transcriptID[i] & gtf$type == "exon", ]
    bed <- bed %>%
      dplyr::select(chr, start, stop, transcriptID, geneID, strand) %>%
      dplyr::mutate(
        eiID = matched$input_id[i],
        score = 0,
        name = paste(transcriptID, "#", matched$input_id[i], sep = ""),
        delta.psi = matched$delta.psi[i],
        p.adj = matched$p.adj[i],
      )

    # Adjust chromStart based on strand to adhere to BED format (0-based) conventions
    bed <- bed %>%
      dplyr::mutate(chromStart = ifelse(strand == "+", as.integer(start) - 1, as.integer(stop) + 1)) %>%
      dplyr::select(chrom = chr, chromStart, chromEnd = stop, name, score, strand, delta.psi, p.adj)

    return(bed)
  }

  toBed <- parallel::mclapply(1:nrow(matched), processTranscriptID_fg, mc.cores = cores)

  # Combine all BED rows into one dataframe
  toBed <- do.call(rbind, toBed)

  return(toBed)
}


bedifyBackground <- function(matched, outname, cores) {

  # Function to process each transcriptID and construct BED rows
  processTranscriptID_bg <- function(i) {
    # Efficiently filter for relevant rows
    bed <- gtf[gtf$transcriptID == matched$transcriptID[i] & gtf$type == "exon", ]
    bed <- bed %>%
      dplyr::select(chr, start, stop, transcriptID, geneID, strand) %>%
      dplyr::mutate(
        eiID = matched$input_id[i],
        score = 0,
        name = paste(transcriptID, "#", matched$input_id[i], sep = "")
      )

    # Adjust chromStart based on strand to adhere to BED format (0-based) conventions
    bed <- bed %>%
      dplyr::mutate(chromStart = ifelse(strand == "+", as.integer(start) - 1, as.integer(stop) + 1)) %>%
      dplyr::select(chrom = chr, chromStart, chromEnd = stop, name, score, strand)

    return(bed)
  }

  toBed <- parallel::mclapply(1:nrow(matched), processTranscriptID_bg, mc.cores = cores)

  # Combine all BED rows into one dataframe
  toBed <- do.call(rbind, toBed)

  return(toBed)
}
