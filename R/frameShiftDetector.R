#' initialize get frame shift with data
#' @param newgtf from setupgtf
#' @param exon_data exon_data from earlier
#' @return multiple dfs with different info and organization related to coding exons
#'
#' @importFrom dplyr group_by mutate arrange select ungroup if_else
#' @keywords internal
getFrameShiftInit <- function(newgtf, exon_data) {
  exon_data <- exon_data[exon_data$ensembl_exon_id %in% newgtf$exonID,]
  exon_data$exon_biotype <- "coding"
  exon_data$exon_biotype[is.na(exon_data$cds_start) & is.na(exon_data$cds_end)] <- "noncoding"
  exon_data$cds_length<- exon_data$cds_end-exon_data$cds_start+1
  exon_data <- exon_data %>%
    dplyr::group_by(ensembl_transcript_id) %>%
    dplyr::mutate(
      adjusted_start = dplyr::if_else(strand == 1, exon_chrom_start, -exon_chrom_start)
    ) %>%
    dplyr::arrange(ensembl_transcript_id, adjusted_start) %>%
    dplyr::select(-adjusted_start) %>%
    dplyr::ungroup()


  coding_exons <- exon_data[exon_data$exon_biotype =="coding",]

  exon_length_df <- exon_data[!(duplicated(exon_data[,c('ensembl_transcript_id', 'ensembl_exon_id', 'cds_length', 'phase', 'end_phase')])),]
  exon_length_df <- left_join(exon_length_df, newgtf[,c('exonID', 'transcriptID', 'classification')], by = c("ensembl_exon_id" = "exonID", "ensembl_transcript_id" = "transcriptID"))

  exon_length_df <- exon_length_df[(!is.na(exon_length_df$cds_start) | !is.na(exon_length_df$cds_end)),]
  exon_length_df$prev_phase <- exon_length_df$phase
  exon_length_df$prev_end_phase <- exon_length_df$end_phase
  exon_length_df$phase <- (exon_length_df$cds_start-1) %% 3
  exon_length_df$end_phase <- (exon_length_df$cds_end-1) %% 3

  return(list(exon_data = exon_data,
              coding_exons = coding_exons,
              exon_length_df = exon_length_df))
}


#' identify frame shifts due to aRNAp events
#' @param fC addinf col
#' @param et aRNAp type
#' @param newgtf newgtf from setupGTF
#' @param exon_data exon data data frame
#' @return a character vector of types of alignmetns of each pair
#'
#' @keywords internal
getFrameShift <- function(fC, et, newgtf, exon_data) {
  gfs_init <- getFrameShiftInit(newgtf, exon_data)
  if (et %in% c("AFE", "HFE")) {
    fs_out <- afheRead(addInf = fC, et,
                       coding_exons = gfs_init$coding_exons,
                       exon_data = gfs_init$exon_data,
                       exon_length_df = gfs_init$exon_length_df,
                       newgtf)
  } else if (et %in% c("ALE", "HLE")) {
    fs_out <- alheRead(addInf = fC, et,
                       coding_exons = gfs_init$coding_exons,
                       exon_data = gfs_init$exon_data,
                       exon_length_df = gfs_init$exon_length_df,
                       newgtf)
  } else if (et %in% c('A3SS')) {
    fs_out <- alt3Read(addInf = fC,
                       coding_exons = gfs_init$coding_exons,
                       exon_data = gfs_init$exon_data,
                       exon_length_df = gfs_init$exon_length_df)
  } else if (et %in% c('A5SS')) {
    fs_out <- alt5Read(addInf = fC,
                       coding_exons = gfs_init$coding_exons,
                       exon_data = gfs_init$exon_data,
                       exon_length_df = gfs_init$exon_length_df)
  } else if (et == 'SE') {
    fs_out <- seRead(addInf = fC,
                     coding_exons = gfs_init$coding_exons,
                     exon_data = gfs_init$exon_data,
                     exon_length_df = gfs_init$exon_length_df,
                     newgtf)
  } else if (et == "RI") {
    fs_out <- irRead(addInf = fC,
                     coding_exons = gfs_init$coding_exons,
                     exon_data = gfs_init$exon_data,
                     exon_length_df = gfs_init$exon_length_df)
  } else if (et == "MXE") {
    fs_out <- mxeRead(addInf = fC,
                      coding_exons = gfs_init$coding_exons,
                      exon_data = gfs_init$exon_data,
                      exon_length_df = gfs_init$exon_length_df,
                      newgtf)
  }
  return(fs_out)
}


#' Get 1st, last, or general overlapping exons between 2 transcripts
#' @param transcript1 transcript 1
#' @param transcript2 transcript 2
#' @param ex exon_type
#' @param coding_exonsX coding exons from initFrameShift
#' @param eld exon_length_df from initFrameShift
#' @return a character vector of types of alignmetns of each pair
#'
#' @keywords internal
getFLOverlap <- function(transcript1, transcript2, ex, coding_exonsX, eld) {

  df1 <- eld[eld$ensembl_transcript_id %in% transcript1 & (!is.na(eld$genomic_coding_end) & !is.na(eld$genomic_coding_start)),]

  df2 <- eld[eld$ensembl_transcript_id %in% transcript2 & (!is.na(eld$genomic_coding_end) & !is.na(eld$genomic_coding_start)),]

  # Find overlaps using vectorized operations
  overlap_matrix <- outer(df1$genomic_coding_start, df2$genomic_coding_end, '<=') & outer(df1$genomic_coding_end, df2$genomic_coding_start, '>=')
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
      return(c(first_overlap_df1$ensembl_exon_id, first_overlap_df2$ensembl_exon_id))
    } else if (ex %in% c('ALE', 'HLE')) {
      last_overlap_indices <- overlap_indices[nrow(overlap_indices), ]

      # Extract the last overlapping pairs
      last_overlap_df1 <- df1[last_overlap_indices[1], ]
      last_overlap_df2 <- df2[last_overlap_indices[2], ]
      return(c(last_overlap_df1$ensembl_exon_id, last_overlap_df2$ensembl_exon_id))
    } else if (ex %in% c('SE', 'MXE')) {
      return(data.frame(t1 = df1$ensembl_exon_id[overlap_indices[,1]],
                        t2 = df2$ensembl_exon_id[overlap_indices[,2]]))
    }


  }

}

#' Get alignment status of last or hybrid last exons
#' @param addInf addInf column passed -- extra info
#' @param et exon_type
#' @param coding_exons coding exons from initFrameShift
#' @param exon_length exon_data from initFrameShift
#' @param exon_length_df exon_length_df from initFrameShift
#' @param newgtf from setup_gtf
#' @return a character vector of types of alignmetns of each pair
#'
#' @keywords internal
alheRead <- function(addInf, et, coding_exons, exon_data, exon_length_df, newgtf) {
  outReads <- unlist(lapply(seq(1, nrow(addInf), by=2), function(x) {
    if (addInf$prot[x] == "none" & addInf$prot[x+1] == "none") {
      return(c("noPC", "noRescue"))
    } else if (sum(c(addInf$prot[x] == "none", addInf$prot[x+1] == "none")) == 1) {
      return(c("onePC", "noRescue"))
    } else {
      return(c("PartialMatch", "noRescue"))
    }
  }))
  return(rep(outReads, each = 2))
}


#' Get alignment status of first or hybrid first exons
#' @param addInf addInf column passed -- extra info
#' @param et exon_type
#' @param coding_exons coding exons from initFrameShift
#' @param exon_length exon_data from initFrameShift
#' @param exon_length_df exon_length_df from initFrameShift
#' @param newgtf from setup_gtf
#' @return a character vector of types of alignmetns of each pair
#'
#' @keywords internal
afheRead <- function(addInf, et, coding_exons, exon_data, exon_length_df, newgtf) {

  outReads <- unlist(lapply(seq(1, nrow(addInf), by=2), function(x) {
    if (addInf$prot[x] == "none" & addInf$prot[x+1] == "none") {
      return(c("noPC", "noRescue"))
    } else if (sum(c(addInf$prot[x] == "none", addInf$prot[x+1] == "none")) == 1) {
      return(c("onePC", "noRescue"))
    }
    overlappingExon <- getFLOverlap(addInf$transcript[x], addInf$transcript[x+1], ex=et, coding_exons, exon_length_df)
    if (overlappingExon[1] == "noOverlap") {
      return(c("PartialMatch", "noRescue"))
    } else {
      if (exon_length_df$strand[exon_length_df$ensembl_exon_id == overlappingExon[1] &
                                exon_length_df$ensembl_transcript_id == addInf$transcript[x]][1] == 1) {
        ph <- c(exon_length_df$genomic_coding_start[exon_length_df$ensembl_exon_id == overlappingExon[1] &
                                                      exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                exon_length_df$genomic_coding_start[exon_length_df$ensembl_exon_id == overlappingExon[2] &
                                                      exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        min.ph <- which.min(ph)
        max.ph <- c(1, 2)[-min.ph]

        phaseHolder <- c(exon_length_df$phase[exon_length_df$ensembl_exon_id == overlappingExon[1] &
                                                exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                         exon_length_df$phase[exon_length_df$ensembl_exon_id == overlappingExon[2] &
                                                exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        phaseHolder[phaseHolder == -1] <- 0

        adj_max <- (((abs(ph[1]-ph[2])) %% 3) + phaseHolder[min.ph]) %% 3
        if (adj_max == phaseHolder[max.ph]) {
          c("PartialMatch", "noRescue")
        } else {
          c("FrameShift",
                 paste0((getRescue(addInf$transcript[x], addInf$transcript[x+1], overlappingExon[1], overlappingExon[2], exon_length_df)), collapse = "#"))
        }
      } else {
        ph <- c(exon_length_df$genomic_coding_end[exon_length_df$ensembl_exon_id == overlappingExon[1] &
                                                    exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                exon_length_df$genomic_coding_end[exon_length_df$ensembl_exon_id == overlappingExon[2] &
                                                    exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        min.ph <- which.min(ph)
        max.ph <- c(1, 2)[-min.ph]

        phaseHolder <- c(exon_length_df$phase[exon_length_df$ensembl_exon_id == overlappingExon[1] &
                                                exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                         exon_length_df$phase[exon_length_df$ensembl_exon_id == overlappingExon[2] &
                                                exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        phaseHolder[phaseHolder == -1] <- 0


        adj_max <- (((abs(ph[1]-ph[2])) %% 3) + phaseHolder[max.ph]) %% 3
        if (adj_max == phaseHolder[min.ph]) {
          c("PartialMatch", "noRescue")
        } else {
          c("FrameShift",
            paste0((getRescue(addInf$transcript[x], addInf$transcript[x+1], overlappingExon[1], overlappingExon[2], exon_length_df)), collapse = "#"))
        }

      }
    }

  }))
  return(rep(outReads, each = 2))
}

#' Get alignment status of introns
#' @param addInf addInf column passed -- extra info
#' @param coding_exons coding exons from initFrameShift
#' @param exon_length exon_data from initFrameShift
#' @param exon_length_df exon_length_df from initFrameShift
#' @return a character vector of types of alignmetns of each pair
#'
#' @keywords internal
irRead <- function(addInf, coding_exons, exon_data, exon_length_df) {
  outReads <- unlist(lapply(seq(1, nrow(addInf), by=2), function(x) {
    if (addInf$prot[x] == "none" & addInf$prot[x+1] == "none") {
      return(c("noPC", "noRescue"))
    } else if (sum(c(addInf$prot[x] == "none", addInf$prot[x+1] == "none")) == 1) {
      return(c("onePC", "noRescue"))
    }

    codeVal <- c(exon_length_df$cds_length[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]] > 0,
                 exon_length_df$cds_length[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]] > 0)
    codeVal[is.na(codeVal)] <- FALSE
    codeVar <- ifelse(sum(codeVal) > 0, ifelse(sum(codeVal) == 2, "allCoding", "someNonCoding"), "allNonCoding")
    if (codeVar == "allNonCoding") {
      return(c("PartialMatch", "noRescue"))
    } else if (codeVar == "someNonCoding") {
      return(c("PartialMatch", "noRescue"))
    } else {

    ir <- exon_length_df$cds_length[exon_length_df$ensembl_transcript_id %in% addInf$transcript[x] & exon_length_df$ensembl_exon_id %in% addInf$exonID[x]]
    sep_ex <- exon_length_df$cds_length[exon_length_df$ensembl_transcript_id %in% addInf$transcript[x+1] & exon_length_df$ensembl_exon_id %in% addInf$exonID[x+1]]
    sep_ex[is.na(sep_ex)] <- 0
    ir[is.na(ir)] <- 0
    if (abs(sum(ir)-sum(sep_ex)) %% 3 == 0) {
      return(c("PartialMatch", "noRescue"))
    } else {
      return(c("FrameShift",
             paste0((getRescue(addInf$transcript[x], addInf$transcript[x+1], addInf$exonID[x], addInf$exonID[x+1], exon_length_df, filterDownstream = TRUE)), collapse = "#")))
    }
      }
  }))
  return(rep(outReads, each = 2))
}

#' Get alignment status of introns
#' @importFrom methods is
#' @param addInf addInf column passed -- extra info
#' @param coding_exons coding exons from initFrameShift
#' @param exon_length exon_data from initFrameShift
#' @param exon_length_df exon_length_df from initFrameShift
#' @return a character vector of types of alignmetns of each pair
#'
#' @keywords internal
mxeRead <- function(addInf, coding_exons, exon_data, exon_length_df, newgtf) {
  outReads <- unlist(lapply(seq(1, nrow(addInf), by=2), function(x) {
    if (addInf$prot[x] == "none" & addInf$prot[x+1] == "none") {
      return(c("noPC", "noRescue"))
    } else if (sum(c(addInf$prot[x] == "none", addInf$prot[x+1] == "none")) == 1) {
      return(c("onePC", "noRescue"))
    }

    codeVal <- c(exon_length_df$cds_length[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]] > 0,
                 exon_length_df$cds_length[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]] > 0)
    codeVal[is.na(codeVal)] <- FALSE
    codeVar <- ifelse(sum(codeVal) > 0, ifelse(sum(codeVal) == 2, "allCoding", "someNonCoding"), "allNonCoding")
    if (codeVar == "allNonCoding") {
      return(c("PartialMatch", "noRescue"))
    } else if (codeVar == "someNonCoding") {
      return(c("PartialMatch", "noRescue"))
    } else {
      downstreamExon <- newgtf$exonID[(which(newgtf$exonID == addInf$exonID[x] &
                                               newgtf$transcriptID == addInf$transcript[x])+1)]

      if (length(downstreamExon) == 0) {
        return(c("PartialMatch", "noRescue"))
      }
      strand <- unique(newgtf$strand[newgtf$transcriptID == addInf$transcript[x]])
      seOv <- getFLOverlap(addInf$transcript[x], addInf$transcript[x+1], ex = "SE", coding_exonsX = coding_exons,
                            eld = exon_length_df)
      if (methods::is(seOv, "character")) {
        return(c("PartialMatch", "noRescue"))
      }
      fs_invest <- seOv[unique(which(seOv$t1 %in% downstreamExon+seOv$t2 %in% downstreamExon > 0)),]
      if (nrow(fs_invest) == 0) {
        return(c("PartialMatch", "noRescue"))
      }
      if (strand == '+') {
        ph <- c(exon_length_df$genomic_coding_start[exon_length_df$ensembl_exon_id == fs_invest$t1[1] &
                                                      exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                exon_length_df$genomic_coding_start[exon_length_df$ensembl_exon_id == fs_invest$t2[1] &
                                                      exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        min.ph <- which.min(ph)
        max.ph <- c(1, 2)[-min.ph]

        phaseHolder <- c(exon_length_df$phase[exon_length_df$ensembl_exon_id == fs_invest$t1[1] &
                                                exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                         exon_length_df$phase[exon_length_df$ensembl_exon_id == fs_invest$t2[1] &
                                                exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        phaseHolder[phaseHolder == -1] <- 0

        adj_max <- (((abs(ph[1]-ph[2])) %% 3) + phaseHolder[min.ph]) %% 3
        if (adj_max == phaseHolder[max.ph]) {
          return(c("PartialMatch", "noRescue"))
        } else {
          return(c("FrameShift",
                   paste0((getRescue(addInf$transcript[x], addInf$transcript[x+1], fs_invest$t1[1], fs_invest$t2[1], exon_length_df)), collapse = "#")))
        }
      } else {
        ph <- c(exon_length_df$genomic_coding_end[exon_length_df$ensembl_exon_id == fs_invest$t1[1] &
                                                    exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                exon_length_df$genomic_coding_end[exon_length_df$ensembl_exon_id == fs_invest$t2[1] &
                                                    exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        min.ph <- which.min(ph)
        max.ph <- c(1, 2)[-min.ph]

        phaseHolder <- c(exon_length_df$phase[exon_length_df$ensembl_exon_id == fs_invest$t1[1] &
                                                exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                         exon_length_df$phase[exon_length_df$ensembl_exon_id == fs_invest$t2[1] &
                                                exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        phaseHolder[phaseHolder == -1] <- 0


        adj_max <- (((abs(ph[1]-ph[2])) %% 3) + phaseHolder[max.ph]) %% 3
        if (adj_max == phaseHolder[min.ph]) {
          return(c("PartialMatch", "noRescue"))
        } else {
          return(c("FrameShift",
                   paste0((getRescue(addInf$transcript[x], addInf$transcript[x+1], fs_invest$t1[1], fs_invest$t2[1], exon_length_df)), collapse = "#")))
        }

      }
    }

  }))
  return(rep(outReads, each = 2))
}

#' Find overlap between two sequences
#' @param start1 sequence 1 start
#' @param end1 sequence 1 end
#' @param start2 sequence 2 start
#' @param end2 sequence 2 end
#' @return which sequences are overlapping
#'
#' @keywords internal
overlap <- function(start1, end1, start2, end2) {
  start_min <- min(start1, end1)
  end_max <- max(start1, end1)
  temp_start2 <- pmin(start2, end2)
  temp_end2 <- pmax(start2, end2)
  return(pmax(start_min, start2) <= pmin(end_max, end2))
}

#' Get alignment status of a3ss
#' @param addInf addInf column passed -- extra info
#' @param coding_exons coding exons from initFrameShift
#' @param exon_length exon_data from initFrameShift
#' @param exon_length_df exon_length_df from initFrameShift
#' @return a character vector of types of alignmetns of each pair
#'
#' @keywords internal
alt3Read <- function(addInf, coding_exons, exon_data, exon_length_df) {
  outReads <- unlist(lapply(seq(1, nrow(addInf), by=2), function(x) {
    if (addInf$prot[x] == "none" & addInf$prot[x+1] == "none") {
      return(c("noPC", "noRescue"))
    } else if (sum(c(addInf$prot[x] == "none", addInf$prot[x+1] == "none")) == 1) {
      return(c("onePC", "noRescue"))
    }


    codeVal <- c(exon_length_df$cds_length[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]] > 0,
                 exon_length_df$cds_length[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]] > 0)
    codeVal[is.na(codeVal)] <- FALSE
    codeVar <- ifelse(sum(codeVal) > 0, ifelse(sum(codeVal) == 2, "allCoding", "someNonCoding"), "allNonCoding")
    if (codeVar == "allNonCoding") {
      return(c("PartialMatch", "noRescue"))
    } else if (codeVar == "someNonCoding") {
      return(c("PartialMatch", "noRescue"))
    } else {


    if (overlap(exon_length_df$genomic_coding_start[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                exon_length_df$genomic_coding_end[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                exon_length_df$genomic_coding_start[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]],
                exon_length_df$genomic_coding_end[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]]))
      {

      if (exon_length_df$strand[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]][1] == 1) {
        ph <- c(exon_length_df$genomic_coding_start[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                exon_length_df$genomic_coding_start[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])

        min.ph <- which.min(ph)
        max.ph <- c(1, 2)[-min.ph]

        phaseHolder <- c(exon_length_df$phase[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                         exon_length_df$phase[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        phaseHolder[phaseHolder == -1] <- 0

        adj_max <- (((abs(ph[1]-ph[2])) %% 3) + phaseHolder[min.ph]) %% 3
        if (adj_max == phaseHolder[max.ph]) {
          return(c("PartialMatch", "noRescue"))
        } else {
          return(c("FrameShift",
                 paste0((getRescue(addInf$transcript[x], addInf$transcript[x+1], addInf$exonID[x], addInf$exonID[x+1], exon_length_df)), collapse = "#")))
        }
      } else {
        ph <- c(exon_length_df$genomic_coding_end[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                exon_length_df$genomic_coding_end[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        min.ph <- which.min(ph)
        max.ph <- c(1, 2)[-min.ph]

        phaseHolder <- c(exon_length_df$phase[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                         exon_length_df$phase[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        phaseHolder[phaseHolder == -1] <- 0


        adj_max <- (((abs(ph[1]-ph[2])) %% 3) + phaseHolder[max.ph]) %% 3
        if (adj_max == phaseHolder[min.ph]) {
          return(c("PartialMatch", "noRescue"))
        } else {
          return(c("FrameShift",
                   paste0((getRescue(addInf$transcript[x], addInf$transcript[x+1], addInf$exonID[x], addInf$exonID[x+1], exon_length_df)), collapse = "#")))
        }

      }

    } else {c("PartialMatch", "noRescue")}
      }
  }))
  return(rep(outReads, each = 2))
}


#' Get alignment status of a5ss
#' @param addInf addInf column passed -- extra info
#' @param coding_exons coding exons from initFrameShift
#' @param exon_length exon_data from initFrameShift
#' @param exon_length_df exon_length_df from initFrameShift
#' @return a character vector of types of alignmetns of each pair
#'
#' @keywords internal
alt5Read <- function(addInf, coding_exons, exon_data, exon_length_df) {
  outReads <- unlist(lapply(seq(1, nrow(addInf), by=2), function(x) {
    if (addInf$prot[x] == "none" & addInf$prot[x+1] == "none") {
      return(c("noPC", "noRescue"))
    } else if (sum(c(addInf$prot[x] == "none", addInf$prot[x+1] == "none")) == 1) {
      return(c("onePC", "noRescue"))
    }

    codeVal <- c(exon_length_df$cds_length[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]] > 0,
                 exon_length_df$cds_length[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]] > 0)
    codeVal[is.na(codeVal)] <- FALSE
    codeVar <- ifelse(sum(codeVal) > 0, ifelse(sum(codeVal) == 2, "allCoding", "someNonCoding"), "allNonCoding")
    if (codeVar == "allNonCoding") {
      return(c("PartialMatch", "noRescue"))
    } else if (codeVar == "someNonCoding") {
      return(c("PartialMatch", "noRescue"))
    } else {

    if (overlap(exon_length_df$genomic_coding_start[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                exon_length_df$genomic_coding_end[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                exon_length_df$genomic_coding_start[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]],
                exon_length_df$genomic_coding_end[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1
                ]])) {

      if (exon_length_df$strand[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]][1] == 1) {
        ph <- c(exon_length_df$genomic_coding_end[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                exon_length_df$genomic_coding_end[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])

        min.ph <- which.min(ph)
        max.ph <- c(1, 2)[-min.ph]

        phaseHolder <- c(exon_length_df$end_phase[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                         exon_length_df$end_phase[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        phaseHolder[phaseHolder == -1] <- 0

        adj_max <- (((abs(ph[1]-ph[2])) %% 3) + phaseHolder[min.ph]) %% 3
        if (adj_max == phaseHolder[max.ph]) {
          return(c("PartialMatch", "noRescue"))
        } else {
          return(c("FrameShift",
                 paste0((getRescue(addInf$transcript[x], addInf$transcript[x+1], addInf$exonID[x], addInf$exonID[x+1], exon_length_df)), collapse = "#")))
        }
      } else {
        ph <- c(exon_length_df$genomic_coding_start[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                exon_length_df$genomic_coding_start[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        min.ph <- which.min(ph)
        max.ph <- c(1, 2)[-min.ph]

        phaseHolder <- c(exon_length_df$end_phase[exon_length_df$ensembl_exon_id == addInf$exonID[x] & exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                         exon_length_df$end_phase[exon_length_df$ensembl_exon_id == addInf$exonID[x+1] & exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        phaseHolder[phaseHolder == -1] <- 0


        adj_max <- (((abs(ph[1]-ph[2])) %% 3) + phaseHolder[max.ph]) %% 3
        if (adj_max == phaseHolder[min.ph]) {
          return(c("PartialMatch", "noRescue"))
        } else {
          return(c("FrameShift",
                   paste0((getRescue(addInf$transcript[x], addInf$transcript[x+1], addInf$exonID[x], addInf$exonID[x+1], exon_length_df)), collapse = "#")))
        }

      }

    } else {return(c("PartialMatch", "noRescue"))}
      }
  }))
  return(rep(outReads, each = 2))
}

#' Get alignment status of se
#' @importFrom methods is
#' @param addInf addInf column passed -- extra info
#' @param coding_exons coding exons from initFrameShift
#' @param exon_length exon_data from initFrameShift
#' @param exon_length_df exon_length_df from initFrameShift
#' @return a character vector of types of alignmetns of each pair
#'
#' @keywords internal
seRead <- function(addInf, coding_exons, exon_data, exon_length_df, newgtf) {
  outReads <- unlist(lapply(seq(1, nrow(addInf), by=2), function(x) {
    if (addInf$prot[x] == "none" & addInf$prot[x+1] == "none") {
      return(c("noPC", "noRescue"))
    } else if (sum(c(addInf$prot[x] == "none", addInf$prot[x+1] == "none")) == 1) {
      return(c("onePC", "noRescue"))
    }
    seExon <- x+grep("Inclusion", addInf$add_inf[c(x, x+1)])-1
    codeVal <- c(exon_length_df$cds_length[exon_length_df$ensembl_exon_id == addInf$exonID[seExon] & exon_length_df$ensembl_transcript_id == addInf$transcript[seExon]] > 0)
    codeVal[is.na(codeVal)] <- FALSE
    codeVal <- ifelse(length(codeVal) == 0, FALSE, codeVal)
    codeVar <- ifelse(sum(codeVal) == 0, "allNonCoding", "allCoding")
    if (codeVar == "allNonCoding") {
      return(c("PartialMatch", "noRescue"))
    } else if (codeVar == "someNonCoding") {
      return(c("PartialMatch", "noRescue"))
    } else {
      downstreamExon <- newgtf$exonID[(which(newgtf$exonID == addInf$exonID[seExon] &
                                               newgtf$transcriptID == addInf$transcript[seExon])+1)]
      if (length(downstreamExon) == 0) {
        return(c("PartialMatch", "noRescue"))
      }
      strand <- unique(newgtf$strand[newgtf$transcriptID == addInf$transcript[seExon]])
      seOv <- getFLOverlap(addInf$transcript[x], addInf$transcript[x+1], ex = "SE", coding_exonsX = coding_exons,
                            eld = exon_length_df)
      if (methods::is(seOv, "character")) {
        return(c("PartialMatch", "noRescue"))
      }
      fs_invest <- seOv[unique(which(seOv$t1 %in% downstreamExon+seOv$t2 %in% downstreamExon > 0)),]
      if (nrow(fs_invest) == 0) {
        return(c("PartialMatch", "noRescue"))
      }
      if (strand == '+') {
        ph <- c(exon_length_df$genomic_coding_start[exon_length_df$ensembl_exon_id == fs_invest$t1[1] &
                                                      exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                exon_length_df$genomic_coding_start[exon_length_df$ensembl_exon_id == fs_invest$t2[1] &
                                                      exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        min.ph <- which.min(ph)
        max.ph <- c(1, 2)[-min.ph]

        phaseHolder <- c(exon_length_df$phase[exon_length_df$ensembl_exon_id == fs_invest$t1[1] &
                                                exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                         exon_length_df$phase[exon_length_df$ensembl_exon_id == fs_invest$t2[1] &
                                                exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        phaseHolder[phaseHolder == -1] <- 0

        adj_max <- (((abs(ph[1]-ph[2])) %% 3) + phaseHolder[min.ph]) %% 3
        if (adj_max == phaseHolder[max.ph]) {
          return(c("PartialMatch", "noRescue"))
        } else {
          return(c("FrameShift",
                   paste0((getRescue(addInf$transcript[x], addInf$transcript[x+1], fs_invest$t1[1], fs_invest$t2[1], exon_length_df)), collapse = "#")))
        }
      } else {
        ph <- c(exon_length_df$genomic_coding_end[exon_length_df$ensembl_exon_id == fs_invest$t1[1] &
                                                    exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                exon_length_df$genomic_coding_end[exon_length_df$ensembl_exon_id == fs_invest$t2[1] &
                                                    exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        min.ph <- which.min(ph)
        max.ph <- c(1, 2)[-min.ph]

        phaseHolder <- c(exon_length_df$phase[exon_length_df$ensembl_exon_id == fs_invest$t1[1] &
                                                exon_length_df$ensembl_transcript_id == addInf$transcript[x]],
                         exon_length_df$phase[exon_length_df$ensembl_exon_id == fs_invest$t2[1] &
                                                exon_length_df$ensembl_transcript_id == addInf$transcript[x+1]])
        phaseHolder[phaseHolder == -1] <- 0


        adj_max <- (((abs(ph[1]-ph[2])) %% 3) + phaseHolder[max.ph]) %% 3
        if (adj_max == phaseHolder[min.ph]) {
          return(c("PartialMatch", "noRescue"))
        } else {
          return(c("FrameShift",
                   paste0((getRescue(addInf$transcript[x], addInf$transcript[x+1], fs_invest$t1[1], fs_invest$t2[1], exon_length_df)), collapse = "#")))
        }

      }


    }
  }
  ))
  return(rep(outReads, each = 2))
}

#' score alignments of proteins
#'
#' @param type event type
#' @param proBed proBed paired output
#' @return align scores and lengths
#' @importFrom msa msaPrettyPrint msa msaClustalW
#' @keywords internal
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
                                                  substitutionMatrix = BLOSUM62, method = "ClustalW"))

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

#' get next downstream exon for MXE, SE
#'
#' @param df1 info related to first transcript
#' @param df2 info related to second transcript
#' @param e1 first exon
#' @param e2 second exon
#' @param eld exon_length_df from previous function
#' @return either new exon set for limiting at least one df
#' @keywords internal
getNextOverlap <- function(df1, df2, e1, e2, eld) {
  if (length(which(df1$ensembl_exon_id == e1)) == 0) {
    e1_out <- e1
  } else {
    e1_out <- df1$ensembl_exon_id[which(df1$ensembl_exon_id == e1)[c(length(which(df1$ensembl_exon_id == e1)))]+1]
  }
  if (length(which(df2$ensembl_exon_id == e2)) == 0) {
    e2_out <- e2
  } else {
    e2_out <- df2$ensembl_exon_id[which(df2$ensembl_exon_id == e2)[c(length(which(df2$ensembl_exon_id == e2)))]+1]
  }

  return(c(e1_out, e2_out))
}

#' check for rescue of the frame shifts
#'
#' @param transcript1 first transcript
#' @param transcript2 second transcript
#' @param e1 first exon
#' @param e2 second exon
#' @param eld exon_length_df from previous function
#' @param filterDownstream param to choose whether or not to filter to downstream of exons
#' @return either noRescues or the transcript and exons that contain rescue
#' @importFrom dplyr mutate relocate
#' @keywords internal
getRescue <- function(transcript1, transcript2, e1, e2, eld, filterDownstream = FALSE) {

  df1 <- eld[eld$ensembl_transcript_id %in% transcript1 &
               (!is.na(eld$genomic_coding_end) & !is.na(eld$genomic_coding_start)),]

  df2 <- eld[eld$ensembl_transcript_id %in% transcript2 &
               (!is.na(eld$genomic_coding_end) & !is.na(eld$genomic_coding_start)),]


  # Find overlaps using vectorized operations
  overlap_matrix <- outer(df1$genomic_coding_start, df2$genomic_coding_end, '<=') &
    outer(df1$genomic_coding_end, df2$genomic_coding_start, '>=')

  if (dim(overlap_matrix)[1] == 0) {
    return('noOverlap')
  }
  # Get indices of the first overlap
  overlap_indices <- which(overlap_matrix, arr.ind = TRUE)

  overlapping_df1 <- df1[overlap_indices[,1], ]
  overlapping_df2 <- df2[overlap_indices[,2], ]

  if (filterDownstream) {
    exons <- getNextOverlap(df1[df1$ensembl_exon_id %in% c(e1, overlapping_df1$ensembl_exon_id),],
                            df2[df2$ensembl_exon_id %in% c(e2, overlapping_df2$ensembl_exon_id),],
                            e1, e2, eld)
    e1 <- exons[1]
    e2 <- exons[2]
    if (any(is.na(c(e1, e2)))) {
      return("noRescue")
    }
  }
  colnames(overlapping_df1) <- paste0(colnames(overlapping_df1), ".x")
  colnames(overlapping_df2) <- paste0(colnames(overlapping_df2), ".y")

  # Combine the dataframes
  overlapping_combined <- cbind(overlapping_df1, overlapping_df2)


  expanded <- overlapping_combined[which(overlapping_combined$ensembl_exon_id.x == e1 &
                                           overlapping_combined$ensembl_exon_id.y == e2)[1]:nrow(overlapping_combined), ]

  expanded$overlap_start <- pmax(expanded$genomic_coding_start.x, expanded$genomic_coding_start.y)
  expanded$overlap_stop <- pmin(expanded$genomic_coding_end.x, expanded$genomic_coding_end.y)


  if (expanded$strand.y[1] == 1) {
    expanded$phase_overlap.x <- (expanded$phase.x+(expanded$overlap_start-expanded$genomic_coding_start.x)) %% 3
    expanded$end_phase_overlap.x <- (expanded$phase_overlap.x + (expanded$overlap_stop-expanded$overlap_start)) %% 3
    expanded$phase_overlap.y <- (expanded$phase.y+(expanded$overlap_start-expanded$genomic_coding_start.y)) %% 3
    expanded$end_phase_overlap.y <- (expanded$phase_overlap.y + (expanded$overlap_stop-expanded$overlap_start)) %% 3
  } else {
    expanded$phase_overlap.x <- (expanded$phase.x + abs(expanded$overlap_stop - expanded$genomic_coding_end.x)) %% 3
    expanded$end_phase_overlap.x <- (expanded$phase_overlap.x + abs(expanded$overlap_start - expanded$overlap_stop)) %% 3
    expanded$phase_overlap.y <- (expanded$phase.y + abs(expanded$overlap_stop - expanded$genomic_coding_end.y)) %% 3
    expanded$end_phase_overlap.y <- (expanded$phase_overlap.y + abs(expanded$overlap_start - expanded$overlap_stop)) %% 3
  }

  expanded$overlap_status_start <- ifelse(expanded$phase_overlap.x == expanded$phase_overlap.y, "same", "diff")
  expanded$overlap_status_stop <- ifelse(expanded$end_phase_overlap.x == expanded$end_phase_overlap.y, "same", "diff")
  expanded$status_start <- ifelse(expanded$phase.x == expanded$phase.y, "same", "diff")
  expanded$status_stop <- ifelse(expanded$end_phase.x == expanded$end_phase.y, "same", "diff")
  gwo <- expanded %>% dplyr::relocate(status_start, .before = overlap_status_start) %>%
    mutate(summarized_fs = case_when(
      status_start == "same" & overlap_status_start == "same" & overlap_status_stop == "same" & status_stop == "diff" ~ "samediff",
      status_start == "same" & overlap_status_start == "same" & overlap_status_stop == "diff" & status_stop == "diff" ~ "samediff",
      status_start == "same" & overlap_status_start == "diff" & overlap_status_stop == "diff" & status_stop == "diff" ~ "samediff",

      status_start == "diff" & overlap_status_start == "same" & overlap_status_stop == "same" & status_stop == "same" ~ "diffsame",
      status_start == "diff" & overlap_status_start == "diff" & overlap_status_stop == "same" & status_stop == "same" ~ "diffsame",
      status_start == "diff" & overlap_status_start == "diff" & overlap_status_stop == "diff" & status_stop == "same" ~ "diffsame",

      status_start == "same" & overlap_status_start == "same" & overlap_status_stop == "same" & status_stop == "same" ~ "same",
      status_start == "diff" & overlap_status_start == "diff" & overlap_status_stop == "diff" & status_stop == "diff" ~ "diff",
      TRUE ~ "other"))

  gwo$summarized_fs[(gwo$prev_end_phase.x == -1 | gwo$prev_end_phase.y == -1) &
                      gwo$status_start == 'same' &
                      gwo$overlap_status_start == 'same' &
                      gwo$overlap_status_stop == 'same' &
                      gwo$status_stop == 'diff' ] <- 'same'

  gwo$summarized_fs[(gwo$prev_end_phase.x == -1 | gwo$prev_end_phase.y == -1) &
                      gwo$status_start == 'diff' &
                      gwo$overlap_status_start == 'diff' &
                      gwo$overlap_status_stop == 'diff' &
                      gwo$status_stop == 'same' ] <- 'diff'

  gwo$summarized_fs[(gwo$prev_end_phase.x == -1 | gwo$prev_end_phase.y == -1) &
                      gwo$status_start == 'same' &
                      gwo$overlap_status_start == 'diff' &
                      gwo$overlap_status_stop == 'diff' &
                      gwo$status_stop == 'same' ] <- 'samediff'


  boolFS <- grepl('same', gwo$summarized_fs)
  boolFS[1] <- ifelse(gwo$summarized_fs[1] == 'samediff', FALSE, boolFS[1])
  if (sum(boolFS) == 0) {
    return("noRescue")
  } else {
    return(c(gwo$ensembl_transcript_id.x[which(boolFS)[1]],
             gwo$ensembl_exon_id.x[which(boolFS)[1]],
             gwo$ensembl_transcript_id.y[which(boolFS)[1]],
             gwo$ensembl_exon_id.y[which(boolFS)[1]]))
  }
}
