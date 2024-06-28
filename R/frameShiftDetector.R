getFrameShiftInit <- function() {
  ensembl <- biomaRt::useEnsembl(biomart = "ensembl",
                                 dataset = "hsapiens_gene_ensembl")
  biomaRt::useEnsembl(biomart = "ensembl",
                               dataset = "hsapiens_gene_ensembl")
  attributes <- c("exon_chrom_start", "exon_chrom_end", "ensembl_exon_id",'cds_start', 'cds_end', 'phase', 'end_phase')
  exon_data <- biomaRt::getBM(attributes = attributes, mart = ensembl)
  exon_data <- exon_data[exon_data$ensembl_exon_id %in% newGTF$gtf$exonID,]
  exon_data$exon_biotype <- "coding"
  exon_data$exon_biotype[is.na(exon_data$cds_start) & is.na(exon_data$cds_end)] <- "noncoding"
  exon_data$cds_length<- exon_data$cds_end-exon_data$cds_start+1
  coding_exons <- exon_data[exon_data$exon_biotype =="coding",]
  exon_length_df <- exon_data[!(duplicated(exon_data[,c('ensembl_exon_id', 'cds_length', 'phase', 'end_phase')])),]
  return(list(exon_data = exon_data,
              coding_exons = coding_exons,
              exon_length_df = exon_length_df))
}


getFrameShift <- function(fC, et) {
  if (et %in% c("afe", "ale", "hfe", "hle")) {
    fs_out <- atheRead(addInf = fC, et)
  } else if (et %in% c('a3ss', 'a5ss')) {
    fs_out <- altsRead(addInf = fC)
  } else if (et == 'se') {
    fs_out <- seRead(addInf = fC)
  } else if (et == "ri") {
    fs_out <- irRead(addInf = fC)
  } else if (et == "mxe") {
    fs_out <- mxeRead(addInf = fC)
  }
  return(fs_out)
}

getFLOverlap <- function(transcript1, transcript2, ex, coding_exons = gfs_init$coding_exons) {

  df1 <- newGTF$gtf[newGTF$gtf$transcriptID %in% transcript1 & newGTF$gtf$exonID %in% coding_exons$ensembl_exon_id,]

  df2 <- newGTF$gtf[newGTF$gtf$transcriptID %in% transcript2 & newGTF$gtf$exonID %in% coding_exons$ensembl_exon_id,]

  # Find overlaps using vectorized operations
  overlap_matrix <- outer(df1$start, df2$stop, '<=') & outer(df1$stop, df2$start, '>=')
  if (dim(overlap_matrix)[1] == 0) {
    return('noOverlap')
  }
  # Get indices of the first overlap
  overlap_indices <- which(overlap_matrix, arr.ind = TRUE)
  if (dim(overlap_indices)[1] == 0) {
    return("noOverlap")
  } else {
    if (ex %in% c('afe', 'hfe')) {
      first_overlap_indices <- overlap_indices[1, ]

      # Extract the first overlapping pairs
      first_overlap_df1 <- df1[first_overlap_indices[1], ]
      first_overlap_df2 <- df2[first_overlap_indices[2], ]
      return(c(first_overlap_df1$exonID, first_overlap_df2$exonID, first_overlap_df1$classification, first_overlap_df2$classification))
    } else if (ex %in% c('ale', 'hle')) {
      last_overlap_indices <- overlap_indices[nrow(overlap_indices), ]

      # Extract the last overlapping pairs
      last_overlap_df1 <- df1[last_overlap_indices[1], ]
      last_overlap_df2 <- df2[last_overlap_indices[2], ]
      return(c(last_overlap_df1$exonID, last_overlap_df2$exonID, last_overlap_df1$classification, last_overlap_df2$classification))
    }


  }

}


atheRead <- function(addInf, et, coding_exons = gfs_init$coding_exons, exon_data = gfs_init$exon_data, exon_length_df = gfs_init$exon_length_df) {

  outReads <- unlist(lapply(seq(1, nrow(addInf), by=2), function(x) {
    if (addInf$prot[x] == "none" & addInf$prot[x+1] == "none") {
      return(paste0("noPC"))
    } else if (sum(c(addInf$prot[x] == "none", addInf$prot[x+1] == "none")) == 1) {
      return(paste0("onePC"))
    }
    overlappingExon <- getFLOverlap(addInf$transcript[x], addInf$transcript[x+1], ex=et)
    if (overlappingExon[1] == "noOverlap") {
      return(paste0("noOverlap"))
    } else {
      e1_over <- unique(exon_length_df$phase[exon_length_df$ensembl_exon_id == overlappingExon[1]])
      e2_over <- unique(exon_length_df$phase[exon_length_df$ensembl_exon_id == overlappingExon[2]])

      e1_under <- unique(exon_length_df$end_phase[exon_length_df$ensembl_exon_id == overlappingExon[1]])
      e2_under <- unique(exon_length_df$end_phase[exon_length_df$ensembl_exon_id == overlappingExon[2]])
      if ("first" %in% overlappingExon[3:4]) {
        return(paste0(ifelse(e1_under == e2_under, "PartialMatch", "FrameShift")))
      } else if ('last' %in% overlappingExon[3:4]) {
        return(paste0(ifelse(e1_over == e2_over, "PartialMatch", "FrameShift")))
      } else {
        return(paste0(ifelse(e1_over == e2_over & e1_under == e2_under, "PartialMatch", "FrameShift")))
      }
    }

  }))
  return(rep(outReads, each = 2))
}


irRead <- function(addInf) {
  outReads <- unlist(mclapply(addInf, mc.cores = 8, function(x) {
    codeVar <- codingCheck(x[3:(length(x)-1)])
    if (codeVar == "allNonCoding#") {
      paste0(codeVar, "PartialMatch")
    } else if (codeVar == "someNonCoding#") {
      paste0(codeVar, "PartialMatch")
    } else {
      ir <- exon_length_df$cds_length[exon_length_df$ensembl_exon_id %in% x[3]]
      sep_ex <- exon_length_df$cds_length[exon_length_df$ensembl_exon_id %in% x[4:(length(x)-1)]]
      sep_ex[is.na(sep_ex)] <- 0
      ir[is.na(ir)] <- 0
      if (abs(sum(ir)-sum(sep_ex)) %% 3 == 0) {
        paste0(codeVar, "PartialMatch")
      } else {
        paste0(codeVar, "FrameShift")
      }
    }
  }))
  return(outReads)
}

mxeRead <- function(addInf) {
  outReads <- unlist(mclapply(addInf, mc.cores = 8, function(x) {
    codeVar <- codingCheck(c(x[(which(x == "MXEI")+1):(which(x == "MXEJ")-1)], x[(which(x == "MXEJ")+1):(length(x)-1)]))
    if (codeVar == "allNonCoding#") {
      paste0(codeVar, "PartialMatch")
    } else if (codeVar == "someNonCoding#") {
      paste0(codeVar, "PartialMatch")
    } else {
      i_lengths <- exon_length_df$cds_length[exon_length_df$ensembl_exon_id %in% x[(which(x == "MXEI")+1):(which(x == "MXEJ")-1)]]
      j_lengths <- exon_length_df$cds_length[exon_length_df$ensembl_exon_id %in% x[(which(x == "MXEJ")+1):(length(x)-1)]]
      i_lengths[is.na(i_lengths)] <- 0
      j_lengths[is.na(j_lengths)] <- 0
      if (abs(sum(i_lengths)-sum(j_lengths)) %% 3 == 0) {
        paste0(codeVar,"PartialMatch")
      } else {
        paste0(codeVar,"FrameShift")
      }
    }
  }))
  return(outReads)
}

altsRead <- function(addInf) {
  outReads <- unlist(mclapply(addInf, mc.cores = 8, function(x) {
    codeVar <- codingCheck(x[c(3, 4)])
    if (codeVar == "allNonCoding#") {
      paste0(codeVar, "PartialMatch")
    } else if (codeVar == "someNonCoding#") {
      paste0(codeVar, "PartialMatch")
    } else {
      lengthsVers1 <- unique(exon_length_df$cds_length[exon_length_df$ensembl_exon_id == x[3]])
      lengthsVers2 <- unique(exon_length_df$cds_length[exon_length_df$ensembl_exon_id == x[4]])
      lengthsVers1[is.na(lengthsVers1)] <- 0
      lengthsVers2[is.na(lengthsVers2)] <- 0
      if (abs(max(lengthsVers1)-max(lengthsVers2)) %% 3 == 0) {
        paste0(codeVar, "PartialMatch")
      } else {
        paste0(codeVar, "FrameShift")
      }
    }
  }))
  return(outReads)
}

seRead <- function(addInf) {
  outReads <- unlist(lapply(seq(1, nrow(addInf), by=2), function(x) {
    if (addInf$prot[x] == "none" & addInf$prot[x+1] == "none") {
      return("noPC")
    } else if (sum(c(addInf$prot[x] == "none", addInf$prot[x+1] == "none")) == 1) {
      return("onePC")
    }
    se <- addInf$exonID[c(x, x+1)][addInf$exonID[c(x, x+1)] %in% newGTF$gtf$exonID[newGTF$gtf$classification == "internal"]]
    lengthsSE <- gfs_init$exon_length_df$cds_length[gfs_init$exon_length_df$ensembl_exon_id %in% se]
    print(lengthsSE)
    lengthsSE[is.na(lengthsSE)] <- 0
    sumLength <- sum(lengthsSE) %% 3
    if (sumLength == 0) {
      return("PartialMatch")
    } else {
      return("FrameShift")
    }
  }
  ))
  return(rep(outReads, each = 2))
}

#' Return both shift and specific results for a given alignment
#' @param df dataframe with protein code column
#' @param i row num
#' @return both specific and shift results
#' @export
frameShiftDetectorSum <- function(df, i) {
  pB_nuc <- unlist(c(fsDirectSpecific(df$code[i], df$code[i+1]),
             fsDirectShift(df$code[i], df$code[i+1])))
  return(pB_nuc)
}

#' Idenfity continuous indels
#' @param indels the indels from an alignment
#' @return continuous indels
#' @export
findContinuousIndels <- function(indels) {
  diffs <- diff(indels)
  breaks <- which(diffs > 1)
  continuousRegions <- split(indels, cumsum(c(1, diff(indels) > 1)))
  continuousRegions
}

#' Align and classify 2 sequences using specific method -- using nucleotides finding indels of non-3 multiples or leads or non-3
#'
#' @return both specific and shift results
#' @param seq1 first nuc sequence
#' @param seq2 second nuc sequence
#' @importFrom msa msa msaConsensusSequence
#' @importFrom Biostrings DNAStringSet
#' @return alignment type
#' @export
fsDirectSpecific <- function(seq1, seq2) {
  alignment <- msa::msaConsensusSequence(msa::msa(Biostrings::DNAStringSet(c(seq1, seq2)),
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

#' Align and classify 2 sequences using shift method,
#'
#' @return both specific and shift results
#' @importFrom msa msa msaConsensusSequence
#' @importFrom Biostrings DNAStringSet pairwiseAlignment
#' @return alignment type
#' @export
fsDirectShift <- function(seq1, seq2) {

  seqUse <- c(seq1, seq2)[which.min(c(nchar(seq1), nchar(seq2)))]
  seqMore <- c(seq1, seq2)[which.max(c(nchar(seq1), nchar(seq2)))]
  a1 <- Biostrings::pairwiseAlignment(pattern = Biostrings::DNAString(substr(seqUse, 1, nchar(seqUse))),
                                      subject = Biostrings::DNAString(seqMore),
                                      substitutionMatrix=matrix(c(2, -1, -1, -1,
                                                                  -1, 2, -1, -1,
                                                                  -1, -1, 2, -1,
                                                                  -1, -1, -1, 2),
                                                                byrow = TRUE, nrow = 4,
                                                                dimnames = list(c("A", "C", "G", "T"),
                                                                                c("A", "C", "G", "T"))),
                                      gapOpening = 10,
                                      gapExtension = 1)
  a2 <- Biostrings::pairwiseAlignment(pattern = Biostrings::DNAString(substr(seqUse, 2, length(seqUse))),
                                      subject = Biostrings::DNAString(seqMore),
                                      substitutionMatrix=matrix(c(2, -1, -1, -1,
                                                                  -1, 2, -1, -1,
                                                                  -1, -1, 2, -1,
                                                                  -1, -1, -1, 2),
                                                                byrow = TRUE, nrow = 4,
                                                                dimnames = list(c("A", "C", "G", "T"),
                                                                                c("A", "C", "G", "T"))),
                                      gapOpening = 10,
                                      gapExtension = 1)@score
  a3 <- Biostrings::pairwiseAlignment(pattern = Biostrings::DNAString(substr(seqUse, 3, length(seqUse))),
                                      subject = Biostrings::DNAString(seqMore),
                                      substitutionMatrix=matrix(c(2, -1, -1, -1,
                                                                  -1, 2, -1, -1,
                                                                  -1, -1, 2, -1,
                                                                  -1, -1, -1, 2),
                                                                byrow = TRUE, nrow = 4,
                                                                dimnames = list(c("A", "C", "G", "T"),
                                                                                c("A", "C", "G", "T"))),
                                      gapOpening = 10,
                                      gapExtension = 1)@score

  a4 <- Biostrings::pairwiseAlignment(pattern = Biostrings::DNAString(substr(seqUse, ((nchar(seqUse)/2)-1.5), nchar(seqUse))),
                                      subject = Biostrings::DNAString(seqMore),
                                      substitutionMatrix=matrix(c(2, -1, -1, -1,
                                                                  -1, 2, -1, -1,
                                                                  -1, -1, 2, -1,
                                                                  -1, -1, -1, 2),
                                                                byrow = TRUE, nrow = 4,
                                                                dimnames = list(c("A", "C", "G", "T"),
                                                                                c("A", "C", "G", "T"))),
                                      gapOpening = 10,
                                      gapExtension = 1)@score
  a5 <- Biostrings::pairwiseAlignment(pattern = Biostrings::DNAString(substr(seqUse, 1+((nchar(seqUse)/2)-1.5), nchar(seqUse))),
                                      subject = Biostrings::DNAString(seqMore),
                                      substitutionMatrix=matrix(c(2, -1, -1, -1,
                                                                  -1, 2, -1, -1,
                                                                  -1, -1, 2, -1,
                                                                  -1, -1, -1, 2),
                                                                byrow = TRUE, nrow = 4,
                                                                dimnames = list(c("A", "C", "G", "T"),
                                                                                c("A", "C", "G", "T"))),
                                      gapOpening = 10,
                                      gapExtension = 1)@score
  a6 <- Biostrings::pairwiseAlignment(pattern = Biostrings::DNAString(substr(seqUse, 2+((nchar(seqUse)/2)-1.5), nchar(seqUse))),
                                      subject = Biostrings::DNAString(seqMore), substitutionMatrix=matrix(c(2, -1, -1, -1,
                                                                                                            -1, 2, -1, -1,
                                                                                                            -1, -1, 2, -1,
                                                                                                            -1, -1, -1, 2),
                                                                                                          byrow = TRUE, nrow = 4,
                                                                                                          dimnames = list(c("A", "C", "G", "T"),
                                                                                                                          c("A", "C", "G", "T"))),
                                      gapOpening = 10,
                                      gapExtension = 1)@score

  leading <- attr(regexpr("^-+", as.character(a1)), "match.length")
  if (leading %% 3 == 0 | leading == -1) {
    if (a1@score < a2 | a1@score < a3 | a4 < a5 | a4 < a6) {
      return("FrameShift")
    } else {return("PartialMatch")}
  } else {return("FrameShift")}
}

alignmentScorer <- function(type, proBed) {
  values <- c(
    4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1, -1, -1, -4,
    -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1, -2,  0, -1, -4,
    -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  4, -3,  0, -1, -4,
    -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4, -3,  1, -1, -4,
    0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4,
    -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0, -2,  4, -1, -4,
    -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1, -3,  4, -1, -4,
    0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -4, -2, -1, -4,
    -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0, -3,  0, -1, -4,
    -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3,  3, -3, -1, -4,
    -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4,  3, -3, -1, -4,
    -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0, -3,  1, -1, -4,
    -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3,  2, -1, -1, -4,
    -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3,  0, -3, -1, -4,
    -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4,
    1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0, -2,  0, -1, -4,
    0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1, -1, -1, -4,
    -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -2, -2, -1, -4,
    -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -1, -2, -1, -4,
    0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3,  2, -2, -1, -4,
    -2, -1,  4,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4, -3,  0, -1, -4,
    -1, -2, -3, -3, -1, -2, -3, -4, -3,  3,  3, -3,  2,  0, -3, -2, -1, -2, -1,  2, -3,  3, -3, -1, -4,
    -1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -2, -2, -2,  0, -3,  4, -1, -4,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4,
    -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1
  )

  # Create the matrix
  matrix_labels <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "J", "Z", "X", "*")
  BLOSUM62 <- matrix(values, nrow = 25, ncol = 25, byrow = TRUE, dimnames = list(matrix_labels, matrix_labels))
  pB_nuc <- proBed
  NonMatch <- unlist(lapply(seq(1, length(pB_nuc$prot), by = 2), function(x) {
    if (sum(pB_nuc$prot[c(x, x+1)] == "none") == 0) {
    c(x, x+1)}
    }))
  pB_nuc2 <- pB_nuc[NonMatch,]
  alignType <- lapply(seq(1,nrow(pB_nuc2), by = 2), function(rowCount) {
    if (as.numeric(nchar(pB_nuc2$prot[rowCount]))*as.numeric(nchar(pB_nuc2$prot[rowCount+1])) >= 2000000000) {
      c(-1, -1)
    } else {
      msaA1 <- msa::msaConsensusSequence(msa::msa(Biostrings::AAStringSet(c(as.character(pB_nuc2$prot[rowCount]),
                                                                as.character(pB_nuc2$prot[rowCount+1]))),
                                                  substitutionMatrix = BLOSUM62))

      alignScore <- sum(strsplit(msaA1, "")[[1]] != "?")/nchar(msaA1)
      c(alignScore, alignScore, nchar(msaA1), nchar(msaA1))
    }
  })
  alignScore <- rep(0, length(pB_nuc$prot))
  alignLength <- rep(0, length(pB_nuc$prot))
  alignScore[NonMatch] <- unlist(lapply(alignType, "[[", 1))
  alignLength[NonMatch] <- unlist(lapply(alignType, "[[", 2))
  return(list(alignScore = alignScore,
         alignLength = alignLength))
}
