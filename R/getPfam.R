## This function associates transcripts from background and foreground datasets with Pfam protein domains,
## using a reference file that maps transcripts to protein codes. It outputs two files containing Pfam domain
## information for the foreground and background datasets.


getPfam <- function(background, foreground, pdir, output_location, cores) {
  # Load the Pfam reference data, which maps transcripts to protein codes
  pfam_hg38 <- readr::read_tsv(paste0(pdir, '/protein_code_from_gencodev43_headerFix.txt.tsv'), col_names = F)

  # Extract transcript IDs from the reference data
  pfam_hg38$transcriptID <- unlist(lapply(strsplit(pfam_hg38$X1, split = "#"),
                                          "[[", 1))

  # Reorder the dataframe to have transcriptID as the first column
  pfam_hg38 <- pfam_hg38 %>% dplyr::relocate(transcriptID)

  # Extract gene IDs from the reference data
  pfam_hg38$geneID <- unlist(lapply(strsplit(unlist(lapply(strsplit(pfam_hg38$X1,
                                                                    split = "#"), "[[", 2)), split = ";"), "[[", 1))
  # Reorder the dataframe to have geneID right after transcriptID
  pfam_hg38 <- pfam_hg38 %>% dplyr::relocate(geneID, .after = transcriptID)

  # Process foreground data to match Pfam domains with transcripts
  fg_out <- do.call(rbind, parallel::mclapply(1:length(foreground$proBed$transcript), mc.cores = cores, function(i) {
    # Check if the transcript is in the Pfam reference
    if (foreground$proBed$transcript[i] %in% pfam_hg38$transcriptID) {
      # Subset the reference data for the matching transcript
      df <- pfam_hg38[pfam_hg38$transcriptID == foreground$proBed$transcript[i],]
      # Construct an identifier combining various elements of the transcript
      id <- paste0(foreground$proBed$transcript[i], "#",
                   foreground$proBed$gene[i], ";",
                   foreground$proBed$chr[i], ":",
                   foreground$proBed$start[i], "-",
                   foreground$proBed$stop[i], ";",
                   foreground$proBed$strand[i])
      # Replace the first column of the dataframe with the new identifier
      df$X1 <- rep(id, length(df$X1))
      df
    }
  }))

  # Process background data similarly to the foreground
  bg_out <- do.call(rbind, parallel::mclapply(1:length(background$proBed$transcript), mc.cores = cores, function(i) {
    if (background$proBed$transcript[i] %in% pfam_hg38$transcriptID) {
      df <- pfam_hg38[pfam_hg38$transcriptID == background$proBed$transcript[i],]
      id <- paste0(background$proBed$transcript[i], "#",
                   background$proBed$gene[i], ";",
                   background$proBed$chr[i], ":",
                   background$proBed$start[i], "-",
                   background$proBed$stop[i], ";",
                   background$proBed$strand[i])
      df$X1 <- rep(id, length(df$X1))
      df
    }
  }))

  # Write the processed foreground and background data to TSV files
  write_tsv(fg_out, paste0(output_location, "Foreground/", "fgoutFast.fa.tsv"))
  write_tsv(bg_out, paste0(output_location, "Background/", "bgoutFast.fa.tsv"))
  return(list(fg_out = fg_out,
              bg_out = bg_out))
}
