bedifyForeground <- function(matched, outname, cores) {

  # Ensure gtf is filtered for "exon" type entries and select relevant columns upfront
  gtf_exons_limited <- gtf %>%
    dplyr::filter(.data$type == "exon") %>%
    dplyr::select(.data$start, .data$stop, .data$transcriptID)

  # Creat a carbon copy of matched with changed column names to avoid ".x" and ".y"
  matched_carbon <- matched %>%
    dplyr::rename(.data$matchedStart = start) %>%
    dplyr::rename(.data$matchedStop = stop)

  # Left join matched data with gtf_exons based on transcriptID
  bed <- dplyr::left_join(matched_carbon, gtf_exons_limited, by = "transcriptID")

  # Manipulate into BED format
  bed <- bed %>%
    dplyr::mutate(
      score = 0,
      name = paste(.data$transcriptID, "#", .data$input_id, sep = ""),
      chromStart = ifelse(.data$strand == "+", bed$start - 1, bed$start),
      chrom = bed$chr,
      chromEnd =ifelse(.data$strand == "+", bed$stop, bed$stop + 1)
    ) %>%
    dplyr::select(
      .data$chrom, .data$chromStart, .data$chromEnd, .data$name, .data$score, .data$strand,
      .data$delta.psi, .data$p.adj, add_inf = .data$add_inf
    )

  return(bed)
}


bedifyBackground <- function(matched, outname, cores) {

  # Ensure gtf is filtered for "exon" type entries and select relevant columns upfront
  gtf_exons_limited <- gtf %>%
    dplyr::filter(.data$type == "exon") %>%
    dplyr::select(.data$start, .data$stop, .data$transcriptID)

  # Creat a carbon copy of matched with changed column names to avoid ".x" and ".y"
  matched_carbon <- matched %>%
    dplyr::rename(matchedStart = .data$start) %>%
    dplyr::rename(matchedStop = .data$stop)

  # Left join matched data with gtf_exons based on transcriptID
  bed <- dplyr::left_join(matched_carbon, gtf_exons_limited, by = "transcriptID")

  # Manipulate into BED format
  bed <- bed %>%
    dplyr::mutate(
      score = 0,
      name = paste(.data$transcriptID, "#", .data$input_id, sep = ""),
      chromStart = ifelse(.data$strand == "+", bed$start - 1, bed$start),
      chrom = bed$chr,
      chromEnd =ifelse(.data$strand == "+", bed$stop, bed$stop + 1)
    ) %>%
    dplyr::select(
      .data$chrom, .data$chromStart, .data$chromEnd, .data$name, .data$score, .data$strand
    )

  return(bed)
}
