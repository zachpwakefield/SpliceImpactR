bedifyForeground <- function(matched, outname, cores) {

  # Ensure gtf is filtered for "exon" type entries and select relevant columns upfront
  gtf_exons_limited <- gtf %>%
    dplyr::filter(type == "exon") %>%
    dplyr::select(start, stop, transcriptID)

  # Creat a carbon copy of matched with changed column names to avoid ".x" and ".y"
  matched_carbon <- matched %>%
    dplyr::rename(matchedStart = start) %>%
    dplyr::rename(matchedStop = stop)

  # Left join matched data with gtf_exons based on transcriptID
  bed <- dplyr::left_join(matched_carbon, gtf_exons_limited, by = "transcriptID")

  # Manipulate into BED format
  bed <- bed %>%
    dplyr::mutate(
      score = 0,
      name = paste(transcriptID, "#", input_id, sep = ""),
      chromStart = ifelse(strand == "+", bed$start - 1, bed$start),
      chrom = bed$chr,
      chromEnd =ifelse(strand == "+", bed$stop, bed$stop + 1)
    ) %>%
    dplyr::select(
      chrom, chromStart, chromEnd, name, score, strand,
      delta.psi, p.adj, add_inf = add_inf
    )

  return(bed)
}


bedifyBackground <- function(matched, outname, cores) {

  # Ensure gtf is filtered for "exon" type entries and select relevant columns upfront
  gtf_exons_limited <- gtf %>%
    dplyr::filter(type == "exon") %>%
    dplyr::select(start, stop, transcriptID)

  # Creat a carbon copy of matched with changed column names to avoid ".x" and ".y"
  matched_carbon <- matched %>%
    dplyr::rename(matchedStart = start) %>%
    dplyr::rename(matchedStop = stop)

  # Left join matched data with gtf_exons based on transcriptID
  bed <- dplyr::left_join(matched_carbon, gtf_exons_limited, by = "transcriptID")

  # Manipulate into BED format
  bed <- bed %>%
    dplyr::mutate(
      score = 0,
      name = paste(transcriptID, "#", input_id, sep = ""),
      chromStart = ifelse(strand == "+", bed$start - 1, bed$start),
      chrom = bed$chr,
      chromEnd =ifelse(strand == "+", bed$stop, bed$stop + 1)
    ) %>%
    dplyr::select(
      chrom, chromStart, chromEnd, name, score, strand,
      add_inf = add_inf
    )

  return(bed)
}
