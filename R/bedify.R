#' Get bed file format from matched transcripts / exons / events
#' for foreground and background separately
#'
#' @param matched from matcher function
#' @param output_location location where everything is being saved
#' @param cores # of cores requested
#' @param gtf gtf dataframe from setup_gtf
#' @return a bed file
#' @importFrom dplyr filter select rename left_join mutate
#' @keywords internal
bedifyForeground <- function(matched, outname, cores, gtf) {

    gtf_exons_limited <- gtf %>%
        dplyr::filter(.data$type == "exon") %>%
        dplyr::select(.data$start, .data$stop, .data$transcriptID)

    matched_carbon <- matched %>%
        dplyr::rename(matchedStart = .data$start) %>%
        dplyr::rename(matchedStop = .data$stop)

    ## Left join matched data with gtf_exons based on transcriptID
    bed <- dplyr::left_join(matched_carbon, gtf_exons_limited,
                            by = "transcriptID")

    ## Manipulate into BED format
    bed <- bed %>%
      dplyr::mutate(
        score = 0,
        name = paste(.data$transcriptID, "#", .data$input_id, sep = ""),
        exonID = .data$exonID,
        chromStart = ifelse(.data$strand == "+", bed$start - 1, bed$start),
        chrom = bed$chr,
        chromEnd =ifelse(.data$strand == "+", bed$stop, bed$stop + 1)
      ) %>%
      dplyr::select(
        .data$chrom, .data$chromStart, .data$chromEnd, .data$name, .data$score,
        .data$strand, .data$delta.psi, .data$p.adj, add_inf = .data$add_inf,
        exonID = .data$exonID
      )


    return(bed)
}

#' Get bed file format from matched transcripts / exons / events
#' for foreground and background separately
#'
#' @param matched from matcher function
#' @param output_location location where everything is being saved
#' @param cores # of cores requested
#' @param gtf gtf dataframe from setup_gtf
#' @return a bed file
#' @importFrom dplyr filter select rename left_join mutate
#' @keywords internal
bedifyBackground <- function(matched, outname, cores, gtf) {

    gtf_exons_limited <- gtf %>%
        dplyr::filter(.data$type == "exon") %>%
        dplyr::select(.data$start, .data$stop, .data$transcriptID)

    matched_carbon <- matched %>%
        dplyr::rename(matchedStart = .data$start) %>%
        dplyr::rename(matchedStop = .data$stop)

    ## Left join matched data with gtf_exons based on transcriptID
    bed <- dplyr::left_join(matched_carbon, gtf_exons_limited,
                            by = "transcriptID")

    ## Manipulate into BED format
    bed <- bed %>%
        dplyr::mutate(
            score = 0,
            name = paste(.data$transcriptID, "#", .data$input_id, sep = ""),
            chromStart = ifelse(.data$strand == "+", bed$start - 1, bed$start),
            chrom = bed$chr,
            chromEnd =ifelse(.data$strand == "+", bed$stop, bed$stop + 1)
    ) %>%
    dplyr::select(
        .data$chrom, .data$chromStart, .data$chromEnd, .data$name,
        .data$score, .data$strand
    )

    return(bed)
}
