matchAlignType <- function(proBed, protCode) {

  # Initialize variables for storing alignment and matching information
  protAlign <- list()  # Stores alignment results
  protC <- c()         # Stores protein category (PC, nonPC, Same, Different)
  pMatch <- c()        # Stores match percentage
  alignType <- c()     # Stores type of alignment (onePC, Match, PartialMatch, FrameShift)


  system(paste0("mkdir ", output_location, "pairedAlignments"))
  setwd(paste0(output_location, "pairedAlignments/"))
  # Iterate through protein codes in pairs
  alignmentTypes <- lapply(seq(1, nrow(df), by = 2), function(i)  {
    pMatch <- 0
    alignN <- 0
    maxPc <- max(nchar(df$prot[i]), nchar(df$prot[i+1]))
    if (df$prot[i] == "none" & df$prot[i+1] == "none") {
      alignType <- c("noPC")
    } else if (df$prot[i] == "none" | df$prot[i+1] == "none") {
      alignType <- c("onePC")
    } else {
      alignN <- sum(strsplit(msa::msaConsensusSequence(msa::msa(Biostrings::AAStringSet(c(df$prot[i], df$prot[i+1])))), "")[[1]] != "?")
      pMatch <- alignN/maxPc
      if (alignN == 1.0) {
        pMatch <- 1.04
        alignType <- c("Match")
      } else if (alignN > .3*maxPc)  {
        alignType <- c("PartialMatch")
        try(msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE,
                           alFile = paste(output_location, "pairedAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_pm_Alignment.fasta", sep = ""),
                           file = paste(output_location, "pairedAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_pm_Alignment.pdf", sep = ""), output = "pdf"))
      } else {
        alignType <- c("FrameShift")
        try(msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE,
                           alFile = paste(output_location, "pairedAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_pm_Alignment.fasta", sep = ""),
                           file = paste(output_location, "pairedAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_fs_Alignment.pdf", sep = ""), output = "pdf"))
      }
    }
    c(alignN, maxPc, pMatch, alignType)

  })
  pMatch <- unlist(lapply(alignmentTypes, "[[", 3))
  alignType <- unlist(lapply(alignmentTypes, "[[", 4))


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
