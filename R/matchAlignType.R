#' getPaired helper function to identify align classifications
#'
#' @param proBed proBed from getForeground in getPaired
#' @param protCode protein code for each protein / transcript / exon associated
#' @param nucleotides get_c_nucs from setup
#' @param saveAlignments whether to save alignment output pdf
#' @param output_location output location
#'
#' @return altered probed and alignment info
#' @importFrom Biostrings AAStringSet
#' @importFrom dplyr left_join
#' @importFrom msa msaPrettyPrint msa msaClustalW
#' @importFrom stats median
#' @keywords internal
matchAlignType <- function(proBed, protCode, nucleotides, output_location = NULL, saveAlignments = FALSE, exon_type, newgtf, exon_data) {

  df <- dplyr::left_join(proBed, nucleotides$transDF, by = c("transcript" = "transcriptID")) # transcripts

  if (!is.null(output_location)) {
  system(paste0("mkdir ", output_location, "pairedAlignments"))
  setwd(paste0(output_location, "pairedAlignments/"))
  }
  # Iterate through protein codes in pairs

  alignmentScores <- alignmentScorer(exon_type, df)
  alignmentScore <- alignmentScores$alignScore
  alignmentScore[alignmentScore == -1] <- rep(stats::median(alignmentScore[alignmentScore != 0]), sum(alignmentScore == -1))

  alignmentTypesIntermediary <- getFrameShift(df, et = exon_type, newgtf, exon_data)
  alignmentTypes <- alignmentTypesIntermediary[seq(1, length(alignmentTypesIntermediary), by = 2)]
  alignmentTypes[alignmentScore == 1] <- "Match"
  if (saveAlignments & !is.null(output_location)) {
    lapply(seq(1, length(alignmentTypes), by = 2), function(i) {
      if (alignmentTypes[i] %in% c("Match", "PartialMatch", "FrameShift")) {
        try(msa::msaPrettyPrint(msa::msa(Biostrings::AAStringSet(c(df$prot[i], df$prot[i+1])), verbose = FALSE, method = "ClustalW"), askForOverwrite=FALSE,
                                file = paste(output_location, "pairedAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1],
                                             "_", alignmentTypes[i], "_Alignment.pdf", sep = ""), output = "pdf"))
      }

    })
  }

  # Add matching information to the 'proBed' dataframe
  df$prop <- alignmentScore
  df$matchType <- alignmentTypes
  df$rescue <- alignmentTypesIntermediary[seq(2, length(alignmentTypesIntermediary), by = 2)]

  # Return a list containing the modified 'proBed' dataframe, match percentages, and alignment types
  return(list(df,
              alignmentScore,
              alignmentTypes))
}
