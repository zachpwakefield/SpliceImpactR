frameShiftDetector <- function(transDF, proBed) {
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
