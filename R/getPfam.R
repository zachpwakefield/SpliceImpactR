#' get the pfam domains for each protein
#'
#' @param foreground output from getForeground
#' @param background output from getBackground
#' @param pdir directory of package
#' @param output_location location to make background directory
#' @param cores number of requested cores
#' @param biomart_data biomart_data from setup stage
#' @return matched pfam data for the foreground and background previously identified
#' @importFrom readr read_tsv write_tsv
#' @importFrom parallel mclapply
#' @importFrom biomaRt useEnsembl getBM
#'
#' @examples
#'
#' pdir <- system.file("extdata", package="SpliceImpactR")
#' dataDirectory <- paste0(pdir, "/")
#' test_group <- paste0(dataDirectory, "rawData/", c("test1","test2", "test3"))
#' control_group <- paste0(dataDirectory, "rawData/", c("control1", "control2", "control3"))
#' data_df <- data.frame(
#'     sample_names = c(control_group, test_group),
#'     phenotype_names = c(
#'       rep("control", length(control_group)),
#'       rep("test", length(test_group))
#'      ),
#'    stringsAsFactors = FALSE
#'   )
#' data_df$utc <- "control"
#' data_df$utc[data_df$phenotype_names == unique(data_df$phenotype_names)[2]] <- "test"
#'
#' transDF <- readr::read_csv(paste0(dataDirectory, "transcripts_limited_transDF.csv"))
#' c_trans <- readr::read_lines(paste0(dataDirectory, "transcripts_limited_c_trans.csv"))
#'
#' transcripts_sample <- list(transDF = transDF,
#'                            c_trans = c_trans)
#'
#' gtf_sample <- list(gtf = readr::read_csv(paste0(dataDirectory, "gtf_limited.csv")),
#'             transcript_gtf = readr::read_csv(paste0(dataDirectory, "transcript_gtf_limited.csv")))
#' translations_sample <- readr::read_lines(paste0(dataDirectory, "translations_limited.csv"))
#'
#' ip <- readr::read_csv(paste0(dataDirectory, "biomart_ip.csv"))
#' code_regions <- readr::read_csv(paste0(dataDirectory, "biomart_code_regions.csv"))
#' pfam_exon_level <- readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon_level.csv"))
#' fsd_exon_data <- readr::read_csv(paste0(dataDirectory, "biomart_data_sample.csv"))
#' pfam_data = readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon.csv"))
#' biomart_data <- list(ip = ip,
#'                      code_regions = code_regions,
#'                      fsd_exon_data = fsd_exon_data,
#'                      pfam_exon_level = pfam_exon_level,
#'                      pfam_data = pfam_data)
#'
#' result <- differential_inclusion_HITindex(test_names = test_group,
#'                                           control_names = control_group,
#'                                           et = "AFE",
#'                                           outlier_threshold = "Inf",
#'                                           minReads = 10,
#'                                           min_prop_samples = 0,
#'                                           chosen_method = "qbGLM"
#'                                           )
#'
#' fg <- getForeground(input = result,
#'                             test_names = test_group,
#'                             control_names = control_group,
#'                             thresh = .1,
#'                             fdr = .05,
#'                             mOverlap = .1,
#'                             exon_type = "AFE",
#'                             output_location = NULL,
#'                             cores = 1,
#'                             gtf = gtf_sample,
#'                             max_zero_prop = 1,
#'                             min_prop_samples = 0,
#'                             translations = translations_sample)
#'
#' bg <- getBackground(input=c(test_group, control_group),
#'                     mOverlap = 0.1,
#'                     cores = 1,
#'                     exon_type = "AFE",
#'                     output_location = NULL, gtf_sample, translations_sample)
#' library(msa)
#' pfg <- getPaired(foreground = fg$proBed,
#'           et = "AFE",
#'           nucleotides = transcripts_sample,
#'           newGTF = gtf_sample,
#'           cores = 1,
#'           output_location = NULL,
#'           saveAlignments = FALSE,
#'           exon_data = biomart_data$fsd_exon_data)
#'
#' pfamData <- getPfam(background = bg,
#'                     foreground = fg,
#'                     pdir,
#'                     output_location = NULL,
#'                     cores = 1,
#'                     biomart_data)
#'
#' @export
getPfam <- function(background, foreground, pdir, output_location = NULL, cores = 1, biomart_data) {

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
  bg_out <- dplyr::left_join(background$proBed, pfam_hg38_transcript, by = c('transcript' = 'transcriptID'), relationship = "many-to-many")
  bg_out <- bg_out[!is.na(bg_out$domains),]

  bg_out <- rbind(bg_out, fg_out[!(fg_out$transcript %in% bg_out$transcript),-c(which(colnames(fg_out) %in% c('delta.psi', 'p.adj', 'add_inf', 'exonID')))])


  if (!is.null(output_location)) {
    write_tsv(fg_out, paste0(output_location, "Foreground/", "fgoutFast.fa.tsv"))
    write_tsv(bg_out, paste0(output_location, "Foreground/", "bgoutFast.fa.tsv"))
  }

  return(list(fg_out = fg_out,
              bg_out = bg_out))
}


#' get the pfam domains and ip conversions
#' @return exon and transcript level domains
#' @importFrom AnnotationDbi mappedkeys
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom PFAM.db PFAMINTERPRO
#' @importFrom dplyr inner_join left_join
#' @keywords internal
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
