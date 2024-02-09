bedifyForeground <- function(matched = matched, outname = outname, cores = 8) {

  # Extract necessary columns from the 'matched' dataframe
  tID <- matched$transcriptID
  eiID <- matched$input_id
  p.adj <- matched$p.adj
  delta.psi <- matched$delta.psi

  toBed <- list() # Initialize an empty list to store BED format data

  toBed <- parallel::mclapply(1:length(tID), mc.cores = cores, function(i) {
    bed <- gtf[gtf$transcriptID == tID[i] & gtf$type == "exon",] %>% dplyr::select(chr, start, stop, transcriptID, geneID, strand)

    # Add extra columns to the BED data: input ID, score (set to 0), constructed 'name' field combining transcriptID and input ID, delta.psi, and adjusted p-value
    bed$eiID <- rep(eiID[i], length(gtf$geneID[gtf$transcriptID == tID[i] & gtf$type == "exon"]))
    bed$score <- 0
    bed$squish <- paste(bed$transcriptID, "#", bed$eiID, sep = "")
    bed$delta.psi <- rep(delta.psi[i], length(gtf$geneID[gtf$transcriptID == tID[i] & gtf$type == "exon"]))
    bed$p.adj <- rep(p.adj[i], length(gtf$geneID[gtf$transcriptID == tID[i] & gtf$type == "exon"]))

    # Rename columns to match BED format specifications and select only necessary columns
    colnames(bed) <- c('chrom','chromStart', 'chromEnd', 'transcriptID', "geneID", "strand", "eiID", "score", "name", "delta.psi", "p.adj")
    bed <- bed %>% dplyr::select(chrom, chromStart, chromEnd, name, score, strand, delta.psi, p.adj)

    # Adjust chromStart based on strand to adhere to BED format (0-based) conventions
    if (unique(bed$strand) == "+") {
      bed$chromStart <- as.integer(bed$chromStart) - 1
    } else {
      bed$chromStart <- as.integer(bed$chromEnd) + 1
    }
    bed # Return the modified BED data for this transcript
  })
  toBed <- do.call(rbind, toBed) # Combine all BED data rows into a single dataframe
  return(toBed) # Return the final BED formatted data
}


bedifyBackground <- function(matched = matched, outname = outname, cores) {

  # Extract necessary columns from the 'matched' dataframe
  tID <- matched$transcriptID
  eiID <- matched$input_id

  toBed <- list() # Initialize an empty list to store BED format data

  toBed <- parallel::mclapply(1:length(tID), mc.cores = cores, function(i) {
    bed <- gtf[gtf$transcriptID == tID[i] & gtf$type == "exon",] %>% dplyr::select(chr, start, stop, transcriptID, geneID, strand)

    # Add extra columns to the BED data: input ID, score (set to 0), constructed 'name' field combining transcriptID and input ID, delta.psi, and adjusted p-value
    bed$eiID <- rep(eiID[i], length(gtf$geneID[gtf$transcriptID == tID[i] & gtf$type == "exon"]))
    bed$score <- 0
    bed$squish <- paste(bed$transcriptID, "#", bed$eiID, sep = "")

    # Rename columns to match BED format specifications and select only necessary columns
    colnames(bed) <- c('chrom','chromStart', 'chromEnd', 'transcriptID', "geneID", "strand", "eiID", "score", "name")
    bed <- bed %>% dplyr::select(chrom, chromStart, chromEnd, name, score, strand)

    # Adjust chromStart based on strand to adhere to BED format (0-based) conventions
    if (unique(bed$strand) == "+") {
      bed$chromStart <- as.integer(bed$chromStart) - 1
    } else {
      bed$chromStart <- as.integer(bed$chromEnd) + 1
    }
    bed # Return the modified BED data for this transcript
  })
  toBed <- do.call(rbind, toBed)
  return(toBed) # Return the final BED formatted data
}
