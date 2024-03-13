frameShiftDetector <- function(transDF, proBed) {
  pB_nuc <- left_join(proBed, transDF, by = c("transcriptID...2" = "transcriptID"), ) # transcript 1
  pB_nuc <- left_join(pB_nuc, transDF, by = c("transcriptID...6" = "transcriptID")) # transcript 2
  pB_nuc <- pB_nuc[pB_nuc$alignType != "onePC",]
  alignType <- parallel::mclapply(1:nrow(pB_nuc), mc.cores = 10, function(rowCount) {
    fsDirect(pB_nuc$code.x[rowCount], pB_nuc$code.y[rowCount])
  })
  return(alignType)
}
findContinuousIndels <- function(indels) {
  diffs <- diff(indels)
  breaks <- which(diffs > 1)
  continuousRegions <- split(indels, cumsum(c(1, diff(indels) > 1)))
  continuousRegions
}
fsInfer <- function(seq1, seq2) {
  f1 <- as.numeric(Biostrings::pairwiseAlignment(seq1, seq2)@score)
  f2 <- as.numeric(Biostrings::pairwiseAlignment(substr(seq1, 2, nchar(seq1)), seq2)@score)
  f3 <- as.numeric(Biostrings::pairwiseAlignment(substr(seq1, 3, nchar(seq1)), seq2)@score)
  return(c("partialMatch", "frameShift", "frameShift")[which.max(c(f1, f2, f3))])
}

fsDirect <- function(seq1, seq2) {
  alignment <- pairwiseAlignment(pattern = DNAString(seq1), subject = DNAString(seq2))
  leading <- attr(regexpr("^-+", as.character(alignment)), "match.length")
  if (leading %% 3 == 0 | leading == -1) {
    checker <- consensusString(alignment)
    indels <- gregexpr("-", checker)[[1]]
    allIndels <- sort(indels)
    allIndels <- allIndels[allIndels != -1]
    continuousIndels <- findContinuousIndels(allIndels)
    frameShifts <- any(sapply(continuousIndels, function(region) length(region) %% 3 != 0))
    if (frameShifts) {
      return("FrameShift")
    } else {return("PartialMatch")}
  } else {return("FrameShift")}
}
# fsDirect <- function(seq1, seq2) {
#   alignment <- pairwiseAlignment(pattern = DNAString(seq1), subject = DNAString(seq2))
#   leading <- attr(regexpr("^-+", as.character(alignment)), "match.length")
#   if (leading %% 3 == 0 | leading == -1) {
#     patternAln <- as.character(alignment@pattern)
#     subjectAln <- as.character(alignment@subject)
#     indelsPattern <- gregexpr("-", patternAln)[[1]]
#     indelsSubject <- gregexpr("-", subjectAln)[[1]]
#     allIndels <- sort(c(indelsPattern, indelsSubject))
#     allIndels <- allIndels[allIndels != -1]
#     continuousIndels <- findContinuousIndels(allIndels)
#     frameShifts <- any(sapply(continuousIndels, function(region) length(region) %% 3 != 0))
#     if (frameShifts) {
#       return("FrameShift")
#     } else {return("PartialMatch")}
#   } else {return("FrameShift")}
# }
