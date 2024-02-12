matchAlignType <- function(proBed, protCode) {
  setwd(output_location)
  # Initialize variables for storing alignment and matching information
  protAlign <- list()  # Stores alignment results
  protC <- c()         # Stores protein category (PC, nonPC, Same, Different)
  pMatch <- c()        # Stores match percentage
  alignType <- c()     # Stores type of alignment (onePC, Match, PartialMatch, FrameShift)


  system(paste0("mkdir ", output_location, "pairedAlignments"))
  # Iterate through protein codes in pairs
  for (i in seq(from=1,to=(length(protCode)-1), by=2)) {

    # Handle cases where one or both protein codes are "none"
    if (protCode[i] == "none" | protCode[i+1] == "none") {

      # Case where the first protein code is "none" and the second is not
      if (protCode[i] == "none" & protCode[i+1] != "none") {
        protC <- c(protC, "nonPC", "PC")
        protAlign[[i]] <- "none"
        protAlign[[i]] <- "onePC"
        alignType <- c(alignType, "onePC")

        # Case where the first protein code is not "none" and the second is
      } else if (protCode[i] != "none" & protCode[i+1] == "none") {
        protC <- c(protC, "PC", "nonPC")
        protAlign[[i]] <- "onePC"
        alignType <- c(alignType, "onePC")

        # Case where both protein codes are "none"
      } else {
        protC <- c(protC, "nonPC", "nonPC")
        protAlign[[i]] <- "none"
        alignType <- c(alignType, "noPC")
      }
      pMatch <- c(pMatch, 0) # Set match percentage to 0 for "none" cases

      # Handle cases where protein codes are identical
    } else if (protCode[i] == protCode[i+1]) {
      protC <- c(protC, "Same", "Same")

      # Perform multiple sequence alignment for identical codes
      protAlign[[i]] <- msa::msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])))
      pMatch <- c(pMatch, 1.04)
      alignType <- c(alignType, "Match")

      # Handle cases where protein codes are different
    } else {
      protC <- c(protC, "Different", "Different")

      # Perform multiple sequence alignment for different codes
      protAlign[[i]] <- msa::msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE)

      # Calculate match percentage based on consensus sequence
      minPc <- min(nchar(protCode[i]), nchar(protCode[i+1]))
      pMatch <- c(pMatch, table(unlist(lapply(strsplit(msa::msaConsensusSequence(protAlign[[i]]), split = ""), function(x) x == "?")))[1]/min(nchar(protCode[i]), nchar(protCode[i+1])))

      # Determine alignment type based on length criteria
      if (nchar(paste(strsplit(msa::msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]][nchar(strsplit(msa::msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]]) > (.1*minPc)], collapse = "")) > .2*minPc)  {
        alignType <- c(alignType, "PartialMatch")
        msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE,
        alFile = paste(output_location, "pairedAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_pm_Alignment.fasta", sep = ""),
        file = paste(output_location, "pairedAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_pm_Alignment.pdf", sep = ""), output = "pdf")
      } else {
        alignType <- c(alignType, "FrameShift")
        msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE,
        alFile = paste(output_location, "pairedAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_pm_Alignment.fasta", sep = ""),
        file = paste(output_location, "pairedAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_fs_Alignment.pdf", sep = ""), output = "pdf")
      }
    }
  }

  # Print a table of protein categories for diagnostic purposes
  print(table(protC))

  # Add matching information to the 'proBed' dataframe
  proBed$match <- protC
  proBed$prop <- rep(pMatch, each = 2)
  proBed$matchType <- rep(alignType, each = 2)

  # Return a list containing the modified 'proBed' dataframe, match percentages, and alignment types
  return(list(proBed,
              pMatch,
              alignType))
}
