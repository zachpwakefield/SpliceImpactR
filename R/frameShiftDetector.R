getFrameShiftInit <- function() {
  ensembl <- biomaRt::useEnsembl(biomart = "ensembl",
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
  gfs_init <- getFrameShiftInit()
  if (et %in% c("AFE", "ALE", "HFE", "HLE")) {
    fs_out <- atheRead(addInf = fC, et,
                       coding_exons = gfs_init$coding_exons,
                       exon_data = gfs_init$exon_data,
                       exon_length_df = gfs_init$exon_length_df)
  } else if (et %in% c('A3SS', 'A5SS')) {
    fs_out <- altsRead(addInf = fC,
                       coding_exons = gfs_init$coding_exons,
                       exon_data = gfs_init$exon_data,
                       exon_length_df = gfs_init$exon_length_df)
  } else if (et == 'SE') {
    fs_out <- seRead(addInf = fC,
                     coding_exons = gfs_init$coding_exons,
                     exon_data = gfs_init$exon_data,
                     exon_length_df = gfs_init$exon_length_df)
  } else if (et == "RI") {
    fs_out <- irRead(addInf = fC,
                     coding_exons = gfs_init$coding_exons,
                     exon_data = gfs_init$exon_data,
                     exon_length_df = gfs_init$exon_length_df)
  } else if (et == "MXE") {
    fs_out <- mxeRead(addInf = fC,
                      coding_exons = gfs_init$coding_exons,
                      exon_data = gfs_init$exon_data,
                      exon_length_df = gfs_init$exon_length_df)
  }
  return(fs_out)
}

getFLOverlap <- function(transcript1, transcript2, ex, coding_exonsX) {

  df1 <- newGTF$gtf[newGTF$gtf$transcriptID %in% transcript1 & newGTF$gtf$exonID %in% coding_exonsX$ensembl_exon_id,]

  df2 <- newGTF$gtf[newGTF$gtf$transcriptID %in% transcript2 & newGTF$gtf$exonID %in% coding_exonsX$ensembl_exon_id,]

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
    if (ex %in% c('AFE', 'HFE')) {
      first_overlap_indices <- overlap_indices[1, ]

      # Extract the first overlapping pairs
      first_overlap_df1 <- df1[first_overlap_indices[1], ]
      first_overlap_df2 <- df2[first_overlap_indices[2], ]
      return(c(first_overlap_df1$exonID, first_overlap_df2$exonID, first_overlap_df1$classification, first_overlap_df2$classification))
    } else if (ex %in% c('ALE', 'HLE')) {
      last_overlap_indices <- overlap_indices[nrow(overlap_indices), ]

      # Extract the last overlapping pairs
      last_overlap_df1 <- df1[last_overlap_indices[1], ]
      last_overlap_df2 <- df2[last_overlap_indices[2], ]
      return(c(last_overlap_df1$exonID, last_overlap_df2$exonID, last_overlap_df1$classification, last_overlap_df2$classification))
    }


  }

}


atheRead <- function(addInf, et, coding_exons, exon_data, exon_length_df) {

  outReads <- unlist(lapply(seq(1, nrow(addInf), by=2), function(x) {
    if (addInf$prot[x] == "none" & addInf$prot[x+1] == "none") {
      return(paste0("noPC"))
    } else if (sum(c(addInf$prot[x] == "none", addInf$prot[x+1] == "none")) == 1) {
      return(paste0("onePC"))
    }
    overlappingExon <- getFLOverlap(addInf$transcript[x], addInf$transcript[x+1], ex=et, coding_exons)
    if (overlappingExon[1] == "noOverlap") {
      return(paste0("PartialMatch"))
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


irRead <- function(addInf, coding_exons, exon_data, exon_length_df) {
  outReads <- unlist(lapply(seq(1, nrow(addInf), by=2), function(x) {
    if (addInf$prot[x] == "none" & addInf$prot[x+1] == "none") {
      return("noPC")
    } else if (sum(c(addInf$prot[x] == "none", addInf$prot[x+1] == "none")) == 1) {
      return("onePC")
    }
      ir <- exon_length_df$cds_length[exon_length_df$ensembl_exon_id %in% addInf$exonID[x]]
      sep_ex <- exon_length_df$cds_length[exon_length_df$ensembl_exon_id %in% addInf$exonID[x+1]]
      sep_ex[is.na(sep_ex)] <- 0
      ir[is.na(ir)] <- 0
      if (abs(sum(ir)-sum(sep_ex)) %% 3 == 0) {
        return("PartialMatch")
      } else {
        return("FrameShift")
      }
  }))
  return(rep(outReads, each = 2))
}

mxeRead <- function(addInf, coding_exons, exon_data, exon_length_df) {
  outReads <- unlist(lapply(seq(1, nrow(addInf), by=2), function(x) {
    if (addInf$prot[x] == "none" & addInf$prot[x+1] == "none") {
      return("noPC")
    } else if (sum(c(addInf$prot[x] == "none", addInf$prot[x+1] == "none")) == 1) {
      return("onePC")
    }
      i_lengths <- exon_length_df$cds_length[exon_length_df$ensembl_exon_id %in% addInf$exonID[x]]
      j_lengths <- exon_length_df$cds_length[exon_length_df$ensembl_exon_id %in% addInf$exonID[x+1]]
      i_lengths[is.na(i_lengths)] <- 0
      j_lengths[is.na(j_lengths)] <- 0
      if (abs(sum(i_lengths)-sum(j_lengths)) %% 3 == 0) {
        return("PartialMatch")
      } else {
        return("FrameShift")
      }

  }))
  return(rep(outReads, each = 2))
}

altsRead <- function(addInf, coding_exons, exon_data, exon_length_df) {
  outReads <- unlist(lapply(seq(1, nrow(addInf), by=2), function(x) {
    if (addInf$prot[x] == "none" & addInf$prot[x+1] == "none") {
      return("noPC")
    } else if (sum(c(addInf$prot[x] == "none", addInf$prot[x+1] == "none")) == 1) {
      return("onePC")
    }
      lengthsVers1 <- unique(exon_length_df$cds_length[exon_length_df$ensembl_exon_id == addInf$exonID[x]])
      lengthsVers2 <- unique(exon_length_df$cds_length[exon_length_df$ensembl_exon_id == addInf$exonID[x+1]])
      lengthsVers1[is.na(lengthsVers1)] <- 0
      lengthsVers2[is.na(lengthsVers2)] <- 0
      if (abs(max(lengthsVers1)-max(lengthsVers2)) %% 3 == 0) {
        return("PartialMatch")
      } else {
        return("FrameShift")
      }
  }))
  return(rep(outReads, each = 2))
}

seRead <- function(addInf, coding_exons, exon_data, exon_length_df) {
  outReads <- unlist(lapply(seq(1, nrow(addInf), by=2), function(x) {
    if (addInf$prot[x] == "none" & addInf$prot[x+1] == "none") {
      return("noPC")
    } else if (sum(c(addInf$prot[x] == "none", addInf$prot[x+1] == "none")) == 1) {
      return("onePC")
    }
    se <- addInf$exonID[c(x, x+1)][addInf$exonID[c(x, x+1)] %in% newGTF$gtf$exonID[newGTF$gtf$classification == "internal"]]
    lengthsSE <- gfs_init$exon_length_df$cds_length[gfs_init$exon_length_df$ensembl_exon_id %in% se]
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
      list(c(-1, -1), c(-1, -1))
    } else {
      msaA1 <- msa::msaConsensusSequence(msa::msa(Biostrings::AAStringSet(c(as.character(pB_nuc2$prot[rowCount]),
                                                                as.character(pB_nuc2$prot[rowCount+1]))),
                                                  substitutionMatrix = BLOSUM62))

      alignScore <- sum(strsplit(msaA1, "")[[1]] != "?")/nchar(msaA1)
      list(c(alignScore, alignScore), c(nchar(msaA1), nchar(msaA1)))
    }
  })
  alignScore <- rep(0, length(pB_nuc$prot))
  alignLength <- rep(0, length(pB_nuc$prot))
  alignScore[NonMatch] <- unlist(lapply(alignType, "[[", 1))
  alignLength[NonMatch] <- unlist(lapply(alignType, "[[", 2))
  return(list(alignScore = alignScore,
         alignLength = alignLength))
}
