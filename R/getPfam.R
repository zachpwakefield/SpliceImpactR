#' get the pfam domains for each protein
#'
#' @param foreground output from getForeground
#' @param background output from getBackground
#' @param pdir directory of package
#' @param output_location location to make background directory
#' @param cores number of requested cores
#' @return figures and dataframes with paired data
#' @importFrom readr read_tsv write_tsv
#' @importFrom parallel mclapply
#' @importFrom biomaRt useEnsembl getBM
#' @export
getPfam <- function(background, foreground, pdir, output_location, cores = 1, biomart_data) {

  iPFAM <- initPFAM(biomart_data)

  pfam_hg38_exon <- data.frame(transcriptID = iPFAM$exon_level$ensembl_transcript_id,
                               exonID = iPFAM$exon_level$ensembl_exon_id,
                               geneID = iPFAM$exon_level$ensembl_gene_id,
                               domains = iPFAM$exon_level$interpro_description,
                               domains_ipID = iPFAM$exon_level$domainIDs,
                               domains_pfamID = iPFAM$exon_level$pfam)

  pfam_hg38_transcript <- data.frame(transcriptID = iPFAM$transcript_level$ensembl_transcript_id,
                                     geneID = iPFAM$transcript_level$ensembl_gene_id,
                                     domains = iPFAM$transcript_level$interpro_description,
                                     domains_ipID = iPFAM$transcript_level$domainIDs,
                                     domains_pfamID = iPFAM$transcript_level$pfam)

  # Process foreground
  foreground$proBed$X1 <- paste0(foreground$proBed$transcript, "#",
                                 foreground$proBed$gene, ";",
                                 foreground$proBed$chr, ":",
                                 foreground$proBed$start, "-",
                                 foreground$proBed$stop, ";",
                                 foreground$proBed$strand)
  fg_out <- dplyr::left_join(foreground$proBed, pfam_hg38_exon, by = c('transcript' = 'transcriptID', 'exonID'))
  fg_out <- fg_out[!is.na(fg_out$domains),]


  # Process background data similarly to the foreground
  background$proBed$X1 <- paste0(background$proBed$transcript, "#",
                                 background$proBed$gene, ";",
                                 background$proBed$chr, ":",
                                 background$proBed$start, "-",
                                 background$proBed$stop, ";",
                                 background$proBed$strand)
  bg_out <- dplyr::left_join(background$proBed, pfam_hg38_transcript, by = c('transcript' = 'transcriptID'))
  bg_out <- bg_out[!is.na(bg_out$domains),]

  bg_out <- rbind(bg_out, fg_out[!(fg_out$transcript %in% bg_out$transcript),-c(which(colnames(fg_out) %in% c('delta.psi', 'p.adj', 'add_inf', 'exonID')))])


  # Write the processed foreground and background data to TSV files
  write_tsv(fg_out, paste0(output_location, "Foreground/", "fgoutFast.fa.tsv"))
  write_tsv(bg_out, paste0(output_location, "Foreground/", "bgoutFast.fa.tsv"))
  return(list(fg_out = fg_out,
              bg_out = bg_out))
}


#' get the pfam domains and ip conversions
#' @return exon and transcript level domains
#' @importFrom AnnotationDbi mappedkeys
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom PFAM.db PFAMINTERPRO
#' @importFrom dplyr inner_join left_join
#' @export
initPFAM <- function(biomart_data) {

  x <- PFAM.db::PFAMINTERPRO
  mapped_keys <- AnnotationDbi::mappedkeys(x)
  xx <- as.list(x[mapped_keys])
  pfam2ipscan <- data.frame(pfam_id = names(xx),
                            domainIDs = unlist(xx))

  ip_convert <- biomart_data$ip[!duplicated(biomart_data$ip[,c(4, 5)]) & biomart_data$ip$interpro_description != "" & biomart_data$ip$interpro != "",c(4, 5)]
  pfam_to_ip <- dplyr::inner_join(pfam2ipscan, ip_convert, by = c('domainIDs' = 'interpro'))

  e <- biomart_data$code_regions
  e$cds_start_aa <- ceiling(e$cds_start/3)
  e$cds_end_aa <- ceiling(e$cds_end/3)


  be8 <- left_join(biomart_data$pfam_exon_level, e, by = c("ensembl_transcript_id", "ensembl_exon_id"))

  be8 <- be8[!is.na(be8$cds_start) & !is.na(be8$cds_end),]
  be8 <- be8[!is.na(be8$pfam_start) & !is.na(be8$pfam_end),]

  pfam_hg38 <- dplyr::inner_join(be8, pfam_to_ip, by = c('pfam' = 'pfam_id'))

  transcript_level <- pfam_hg38[!duplicated(pfam_hg38[,c(8, 2, 5, 6, 7, 13, 14)]), c(8, 2, 5, 6, 7, 13, 14)]
  exon_level <- pfam_hg38[with(pfam_hg38, (pfam_start <= cds_end_aa) & (pfam_end >= cds_start_aa)),c(8, 1, 2, 5, 13, 14)]
  return(list(transcript_level = transcript_level,
              exon_level = exon_level))
}
