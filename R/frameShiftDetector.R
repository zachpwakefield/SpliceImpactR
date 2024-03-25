frameShiftDetectorSum <- function(df, i) {
  unlist(c(fsDirectSpecific(df$code[i], df$code[i+1]),
             fsDirectShift(df$code[i], df$code[i+1])))
  return(pB_nuc)
}


findContinuousIndels <- function(indels) {
  diffs <- diff(indels)
  breaks <- which(diffs > 1)
  continuousRegions <- split(indels, cumsum(c(1, diff(indels) > 1)))
  continuousRegions
}

fsDirectSpecific <- function(seq1, seq2) {
  alignment <- msa::msaConsensusSequence(msa::msa(DNAStringSet(c(seq1, seq2)),
                                                  substitutionMatrix =
                                                    matrix(c(2, -1, -1, -1,
                                                             -1, 2, -1, -1,
                                                             -1, -1, 2, -1,
                                                             -1, -1, -1, 2),
                                                           byrow = TRUE, nrow = 4,
                                                           dimnames = list(c("A", "C", "G", "T"),
                                                                           c("A", "C", "G", "T"))),
                                                  gapOpening = 10, gapExtension = 1))
  alignScore <- sum(strsplit(alignment, "")[[1]] != "?")/max(nchar(seq1), nchar(seq2))
  leading <- attr(regexpr("^[?]+", alignment), "match.length")
  if (leading %% 3 == 0 | leading == -1) {
    indels <- gregexpr("[?]", alignment)[[1]]
    allIndels <- sort(indels)
    allIndels <- allIndels[allIndels != -1]
    continuousIndels <- findContinuousIndels(allIndels)
    frameShifts <- any(sapply(continuousIndels, function(region) length(region) %% 3 != 0))
    if (frameShifts) {
      return(c("FrameShift", alignScore))
    } else {return(c("PartialMatch", alignScore))}
  } else {return(c("FrameShift", alignScore))}
}

fsDirectShift <- function(seq1, seq2) {

  seqUse <- c(seq1, seq2)[which.min(c(nchar(seq1), nchar(seq2)))]
  seqMore <- c(seq1, seq2)[which.max(c(nchar(seq1), nchar(seq2)))]
  a1 <- Biostrings::pairwiseAlignment(pattern = Biostrings::DNAString(substr(seqUse, 1, nchar(seqUse))), subject = Biostrings::DNAString(seqMore), substitutionMatrix=matrix(c(2, -1, -1, -1,
                                                                                                                                                                               -1, 2, -1, -1,
                                                                                                                                                                               -1, -1, 2, -1,
                                                                                                                                                                               -1, -1, -1, 2),
                                                                                                                                                                             byrow = TRUE, nrow = 4,
                                                                                                                                                                             dimnames = list(c("A", "C", "G", "T"),
                                                                                                                                                                                             c("A", "C", "G", "T"))), gapOpening = 10, gapExtension = 1)
  a2 <- Biostrings::pairwiseAlignment(pattern = Biostrings::DNAString(substr(seqUse, 2, length(seqUse))), subject = Biostrings::DNAString(seqMore), substitutionMatrix=matrix(c(2, -1, -1, -1,
                                                                                                                                                                                -1, 2, -1, -1,
                                                                                                                                                                                -1, -1, 2, -1,
                                                                                                                                                                                -1, -1, -1, 2),
                                                                                                                                                                              byrow = TRUE, nrow = 4,
                                                                                                                                                                              dimnames = list(c("A", "C", "G", "T"),
                                                                                                                                                                                              c("A", "C", "G", "T"))), gapOpening = 10, gapExtension = 1)@score
  a3 <- Biostrings::pairwiseAlignment(pattern = Biostrings::DNAString(substr(seqUse, 3, length(seqUse))), subject = Biostrings::DNAString(seqMore), substitutionMatrix=matrix(c(2, -1, -1, -1,
                                                                                                                                                                                -1, 2, -1, -1,
                                                                                                                                                                                -1, -1, 2, -1,
                                                                                                                                                                                -1, -1, -1, 2),
                                                                                                                                                                              byrow = TRUE, nrow = 4,
                                                                                                                                                                              dimnames = list(c("A", "C", "G", "T"),
                                                                                                                                                                                              c("A", "C", "G", "T"))), gapOpening = 10, gapExtension = 1)@score

  a4 <- Biostrings::pairwiseAlignment(pattern = Biostrings::DNAString(substr(seqUse, ((nchar(seqUse)/2)-1.5), nchar(seqUse))), subject = Biostrings::DNAString(seqMore), substitutionMatrix=matrix(c(2, -1, -1, -1,
                                                                                                                                                                                                     -1, 2, -1, -1,
                                                                                                                                                                                                     -1, -1, 2, -1,
                                                                                                                                                                                                     -1, -1, -1, 2),
                                                                                                                                                                                                   byrow = TRUE, nrow = 4,
                                                                                                                                                                                                   dimnames = list(c("A", "C", "G", "T"),
                                                                                                                                                                                                                   c("A", "C", "G", "T"))), gapOpening = 10, gapExtension = 1)@score
  a5 <- Biostrings::pairwiseAlignment(pattern = Biostrings::DNAString(substr(seqUse, 1+((nchar(seqUse)/2)-1.5), nchar(seqUse))), subject = Biostrings::DNAString(seqMore), substitutionMatrix=matrix(c(2, -1, -1, -1,
                                                                                                                                                                                                       -1, 2, -1, -1,
                                                                                                                                                                                                       -1, -1, 2, -1,
                                                                                                                                                                                                       -1, -1, -1, 2),
                                                                                                                                                                                                     byrow = TRUE, nrow = 4,
                                                                                                                                                                                                     dimnames = list(c("A", "C", "G", "T"),
                                                                                                                                                                                                                     c("A", "C", "G", "T"))), gapOpening = 10, gapExtension = 1)@score
  a6 <- Biostrings::pairwiseAlignment(pattern = Biostrings::DNAString(substr(seqUse, 2+((nchar(seqUse)/2)-1.5), nchar(seqUse))), subject = Biostrings::DNAString(seqMore), substitutionMatrix=matrix(c(2, -1, -1, -1,
                                                                                                                                                                                                       -1, 2, -1, -1,
                                                                                                                                                                                                       -1, -1, 2, -1,
                                                                                                                                                                                                       -1, -1, -1, 2),
                                                                                                                                                                                                     byrow = TRUE, nrow = 4,
                                                                                                                                                                                                     dimnames = list(c("A", "C", "G", "T"),
                                                                                                                                                                                                                     c("A", "C", "G", "T"))), gapOpening = 10, gapExtension = 1)@score

  leading <- attr(regexpr("^-+", as.character(a1)), "match.length")
  if (leading %% 3 == 0 | leading == -1) {
    if (a1@score < a2 | a1@score < a3 | a4 < a5 | a4 < a6) {
      return("FrameShift")
    } else {return("PartialMatch")}
  } else {return("FrameShift")}
}
