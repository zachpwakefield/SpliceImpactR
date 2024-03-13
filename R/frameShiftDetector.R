frameShiftDetector <- function(transDF, proBed) {
  pB_nuc <- left_join(proBed, transDF, by = c("transcriptID...2" = "transcriptID"), ) # transcript 1
  pB_nuc <- left_join(pB_nuc, transDF, by = c("transcriptID...6" = "transcriptID")) # transcript 2
  pB_nuc <- pB_nuc[pB_nuc$alignType != "onePC",]
  alignType <- parallel::mclapply(1:nrow(pB_nuc), mc.cores = 10, function(rowCount) {
    fsDirect(pB_nuc$code.x[rowCount], pB_nuc$code.y[rowCount])
  })
}

fsInfer <- function(seq1, seq2) {
  f1 <- as.numeric(Biostrings::pairwiseAlignment(seq1, seq2)@score)
  f2 <- as.numeric(Biostrings::pairwiseAlignment(substr(seq1, 2, nchar(seq1)), seq2)@score)
  f3 <- as.numeric(Biostrings::pairwiseAlignment(substr(seq1, 3, nchar(seq1)), seq2)@score)
  return(c("partialMatch", "frameShift", "frameShift")[which.max(c(f1, f2, f3))])
}

fsDirect <- function() {

}
fsDirect <- function(seq1, seq2) {
  # Perform pairwise alignment
  alignment <- Biostrings::pairwiseAlignment(pattern = Biostrings::DNAString(seq1),
                                             subject = Biostrings::DNAString(seq2))

  patternAlignment <- as.character(alignment@pattern)
  subjectAlignment <- as.character(alignment@subject)

  # Count the number of '-' characters which represent gaps (indels)
  totalIndels <- sum(strsplit(patternAlignment, "")[[1]] == "-")+sum(strsplit(subjectAlignment, "")[[1]] == "-")
  return(c("FrameShift", "Match", "PartialMatch")[c(totalIndels %% 3 == 0, totalIndels > 0, totalIndels == 0)])
}
