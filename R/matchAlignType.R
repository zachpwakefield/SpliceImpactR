matchAlignType <- function(proBed, protCode) {
  protAlign <- list()
  protC <- c()
  pMatch <- c()
  alignType <- c()
  cate <- c()
  for (i in seq(from=1,to=(length(protCode)-1), by=2)) {
    if (protCode[i] == "none" | protCode[i+1] == "none") {
      if (protCode[i] == "none" & protCode[i+1] != "none") {
        protC <- c(protC, "nonPC", "PC")
        protAlign[[i]] <- "none"
        protAlign[[i]] <- "onePC"
        alignType <- c(alignType, "onePC")
      } else if (protCode[i] != "none" & protCode[i+1] == "none") {
        protC <- c(protC, "PC", "nonPC")
        protAlign[[i]] <- "onePC"
        alignType <- c(alignType, "onePC")
      } else {
        protC <- c(protC, "nonPC", "nonPC")
        protAlign[[i]] <- "none"
        alignType <- c(alignType, "noPC")
      }
      pMatch <- c(pMatch, 0)
    } else if (protCode[i] == protCode[i+1]) {
      protC <- c(protC, "Same", "Same")
      protAlign[[i]] <- msa::msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])))
      pMatch <- c(pMatch, 1.04)
      alignType <- c(alignType, "Match")
    } else {
      protC <- c(protC, "Different", "Different")
      protAlign[[i]] <- msa::msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE)

      minPc <- min(nchar(protCode[i]), nchar(protCode[i+1]))
      pMatch <- c(pMatch, table(unlist(lapply(strsplit(msa::msaConsensusSequence(protAlign[[i]]), split = ""), function(x) x == "?")))[1]/min(nchar(protCode[i]), nchar(protCode[i+1])))
      if (nchar(paste(strsplit(msa::msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]][nchar(strsplit(msa::msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]]) > (.1*minPc)], collapse = "")) > .2*minPc)  {
        alignType <- c(alignType, "PartialMatch")
        # msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE
        # , file = paste(out_dir, "prettyAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_pm_prettyAlignment.pdf", sep = ""), output = "pdf")
      } else {
        alignType <- c(alignType, "FrameShift")
        # msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE
        # , file = paste(out_dir, "prettyAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_fs_prettyAlignment.pdf", sep = ""), output = "pdf")
      }
    }
  }
  print(table(protC))
  proBed$match <- protC
  proBed$prop <- rep(pMatch, each = 2)
  proBed$matchType <- rep(alignType, each = 2)
  return(list(proBed,
              pMatch,
              alignType))
}
