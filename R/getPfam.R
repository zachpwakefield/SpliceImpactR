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
getPfam <- function(background, foreground, pdir, output_location, cores = 1) {

  iPFAM <- initPFAM()

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
  fg_out <- left_join(foreground$proBed, pfam_hg38_transcript, by = c('transcript' = 'transcriptID'))
  fg_out <- fg_out$domains[!is.na(fg_out$domains)]


  # Process background data similarly to the foreground
  background$proBed$X1 <- paste0(background$proBed$transcript, "#",
                                 background$proBed$gene, ";",
                                 background$proBed$chr, ":",
                                 background$proBed$start, "-",
                                 background$proBed$stop, ";",
                                 background$proBed$strand)
  bg_out <- left_join(background$proBed, pfam_hg38_transcript, by = c('transcript' = 'transcriptID'))
  bg_out <- bg_out$domains[!is.na(bg_out$domains)]

  bg_out <- rbind(bg_out, fg_out[!(fg_out$transcriptID %in% bg_out$transcriptID),])


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
initPFAM <- function() {

  x <- PFAM.db::PFAMINTERPRO
  mapped_keys <- AnnotationDbi::mappedkeys(x)
  xx <- as.list(x[mapped_keys])
  pfam2ipscan <- data.frame(pfam_id = names(xx),
                            domainIDs = unlist(xx))

  attributes <- c("ensembl_transcript_id", "chromosome_name",
                  "transcript_biotype","interpro_description", "interpro")
  ip <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(1:23, "X", "Y"), "protein_coding"), filters = c('chromosome_name', "transcript_biotype"))
  ip_convert <- ip[!duplicated(ip[,c(4, 5)]) & ip$interpro_description != "" & ip$interpro != "",c(4, 5)]
  pfam_to_ip <- dplyr::inner_join(pfam2ipscan, ip_convert, by = c('domainIDs' = 'interpro'))

  attributes <- c("ensembl_transcript_id", "ensembl_exon_id", "pfam", 'chromosome_name', 'transcript_biotype', 'pfam_start', 'pfam_end')
  b1 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(1:2), "protein_coding"), filters = c('chromosome_name', "transcript_biotype"))
  b2 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(3:4), "protein_coding"), filters = c('chromosome_name', "transcript_biotype"))
  b3 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(5:7), "protein_coding"), filters = c('chromosome_name', "transcript_biotype"))
  b4 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(8:10), "protein_coding"), filters = c('chromosome_name', "transcript_biotype"))
  b5 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(11:13), "protein_coding"), filters = c('chromosome_name', "transcript_biotype"))
  b6 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(14:16), "protein_coding"), filters = c('chromosome_name', "transcript_biotype"))
  b7 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(17:20), "protein_coding"), filters = c('chromosome_name', "transcript_biotype"))
  b8 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(21:23), "protein_coding"), filters = c('chromosome_name', "transcript_biotype"))
  b9 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c("X", "Y"), "protein_coding"), filters = c('chromosome_name', "transcript_biotype"))

  attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id", "cds_start", "cds_end")
  e <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list("protein_coding"), filters = c("transcript_biotype"))
  e$cds_start_aa <- ceiling(e$cds_start/3)
  e$cds_end_aa <- ceiling(e$cds_end/3)

  pfam_exon_level <- do.call(rbind, list(b1, b2, b3, b4, b5, b6, b7, b8, b9))

  be8 <- left_join(pfam_exon_level, e, by = c("ensembl_transcript_id", "ensembl_exon_id"))

  be8 <- be8[!is.na(be8$cds_start) & !is.na(be8$cds_end),]
  be8 <- be8[!is.na(be8$pfam_start) & !is.na(be8$pfam_end),]

  pfam_hg38 <- dplyr::inner_join(be8, pfam_to_ip, by = c('pfam' = 'pfam_id'))

  transcript_level <- pfam_hg38[!duplicated(pfam_hg38[,c(8, 2, 5, 6, 7, 13, 14)]), c(8, 2, 5, 6, 7, 13, 14)]
  exon_level <- pfam_hg38[with(pfam_hg38, (pfam_start <= cds_end_aa) & (pfam_end >= cds_start_aa)),c(8, 1, 2, 5, 13, 14)]
  return(list(transcript_level = transcript_level,
              exon_level = exon_level))
}
