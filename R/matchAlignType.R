matchAlignType <- function(proBed, protCode, nucleotides) {

  df <- dplyr::left_join(proBed, nucleotides$transDF, by = c("transcript" = "transcriptID")) # transcripts

  system(paste0("mkdir ", output_location, "pairedAlignments"))
  setwd(paste0(output_location, "pairedAlignments/"))
  # Iterate through protein codes in pairs
  alignmentTypes <- lapply(seq(1, nrow(df), by = 2), function(i)  {
    if (df$prot[i] == "none" & df$prot[i+1] == "none") {
      alignType <- c("noPC")
      pMatch <- 0
    } else if (df$prot[i] == "none" | df$prot[i+1] == "none") {
      alignType <- c("onePC")
      pMatch <- 0
    } else if (df$prot[i] == df$prot[i+1]) {
      pMatch <- 1.04
      alignType <- c("Match")
    } else {
      fs_check <- frameShiftDetectorSum(df, i)
      pMatch <- fs_check[[2]]
      if (alignType[[1]] == alignType[[3]]) {
        alignType <- fs_check[[1]]
      } else {
        alignType <- c("PartialMatch")
      }

      try(msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE,
                         alFile = paste(output_location, "pairedAlignments/", fs_check, "_", proBed$transcript[i], "_", proBed$transcript[i+1], "_pm_Alignment.fasta", sep = ""),
                         file = paste(output_location, "pairedAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_pm_Alignment.pdf", sep = ""), output = "pdf"))
    }



    c(pMatch, alignType)

  })
  pMatch <- unlist(lapply(alignmentTypes, "[[", 1))
  alignType <- unlist(lapply(alignmentTypes, "[[", 2))


  # Print a table of protein categories for diagnostic purposes
  print(table(alignType))

  # Add matching information to the 'proBed' dataframe
  proBed$prop <- rep(pMatch, each = 2)
  proBed$matchType <- rep(alignType, each = 2)

  # Return a list containing the modified 'proBed' dataframe, match percentages, and alignment types
  return(list(proBed,
              pMatch,
              alignType))
}
