frameShiftDetectorInfer <- function(transDF, proBed) {
  pB_nuc <- left_join(proBed, transDF, by = c("transcriptID...2" = "transcriptID"), ) # transcript 1
  pB_nuc <- left_join(pB_nuc, transDF, by = c("transcriptID...6" = "transcriptID")) # transcript 2
  pB_nuc <- pB_nuc[pB_nuc$alignType != "onePC",]
  alignType <- parallel::mclapply(1:nrow(pB_nuc), mc.cores = 10, function(rowCount) {
    f1 <- as.numeric(Biostrings::pairwiseAlignment(pB_nuc$code.x[rowCount], pB_nuc$code.y[rowCount])@score)
    f2 <- as.numeric(Biostrings::pairwiseAlignment(substr(pB_nuc$code.x[rowCount], 2, nchar(pB_nuc$code.x[rowCount])), pB_nuc$code.y[rowCount])@score)
    f3 <- as.numeric(Biostrings::pairwiseAlignment(substr(pB_nuc$code.x[rowCount], 3, nchar(pB_nuc$code.x[rowCount])), pB_nuc$code.y[rowCount])@score)
    c("partialMatch", "frameShift", "frameShift")[which.max(c(f1, f2, f3))]
  })
}
frameShiftDetectorDirect <- function(seq1, seq2) {
  dnaSeq1 <- Biostrings::DNAString(seq1)
  dnaSeq2 <- Biostrings::DNAString(seq2)

  # Perform pairwise alignment
  alignment <- Biostrings::pairwiseAlignment(pattern = dnaSeq1, subject = dnaSeq2)

  patternAlignment <- as.character(alignment@pattern)
  subjectAlignment <- as.character(alignment@subject)

  # Count the number of '-' characters which represent gaps (indels)
  indelsPattern <- sum(strsplit(patternAlignment, "")[[1]] == "-")
  indelsSubject <- sum(strsplit(subjectAlignment, "")[[1]] == "-")
  if (indelsPattern > 0 || indelsSubject > 0) {
    totalIndels <- indelsPattern + indelsSubject
    if (totalIndels %% 3 == 0) {
      return("PartialMatch")
    } else {
      return("FrameShift")
    }
  } else {
    return("Match")
  }
}
