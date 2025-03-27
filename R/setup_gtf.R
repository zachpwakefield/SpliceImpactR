#' setting up gtf with exon labels, first, internal, last, hybrid, etc
#' @importFrom dplyr filter group_by mutate ungroup inner_join select case_when n join_by
#' @importFrom biomaRt useEnsembl getBM
#' @return annotated gtf and various hybrid info. also gene to transcript to protein naming dataframe
#' @keywords internal

setupAnnotation <- function(biomart_data) {

  exon_data <- biomart_data$setup_gtf_exon_data
  exon_data_exon_label <- exon_data %>% dplyr::filter(chromosome_name %in% c(seq_len(23), "X", "Y")) %>%
    dplyr::group_by(ensembl_gene_id, ensembl_transcript_id) %>%
    dplyr::mutate(exon_label = dplyr::case_when(
      rank == min(rank) ~ "first",
      rank == max(rank) ~ "last",
      TRUE ~ "internal"
    )) %>%
    dplyr::ungroup()


  first_exons <- exon_data_exon_label %>% dplyr::filter(exon_label == "first")
  internal_exons <- exon_data_exon_label %>% dplyr::filter(exon_label == "internal")
  last_exons <- exon_data_exon_label %>% dplyr::filter(exon_label == "last")

  # Join the first and internal exons to find overlaps within the same gene
  hybrid_first_exons <- dplyr::inner_join(first_exons, exon_data_exon_label, by = "ensembl_gene_id", suffix = c("_first", "_internal")) %>%
    dplyr::filter(chromosome_name_first == chromosome_name_internal &
                    ensembl_transcript_id_first != ensembl_transcript_id_internal &
                    ((exon_chrom_start_first <= exon_chrom_end_internal & exon_chrom_end_first >= exon_chrom_start_internal) |
                       (exon_chrom_start_internal <= exon_chrom_end_first & exon_chrom_end_internal >= exon_chrom_start_first))) %>%
    dplyr::group_by(ensembl_transcript_id_first, ensembl_transcript_id_internal) %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::ungroup() %>% filter(exon_label_internal == "internal")

  # Select relevant columns to display the hybrid exons
  hybrid_first_exons <- hybrid_first_exons %>% dplyr::select(ensembl_gene_id, ensembl_transcript_id_first, exon_chrom_start_first, exon_chrom_end_first, ensembl_exon_id_first,
                                                             ensembl_transcript_id_internal, exon_chrom_start_internal, exon_chrom_end_internal, ensembl_exon_id_internal)

  # Join the first and internal exons to find overlaps within the same gene
  hybrid_last_exons <- dplyr::inner_join(last_exons, exon_data_exon_label, by = "ensembl_gene_id", suffix = c("_last", "_internal")) %>%
    dplyr::filter(chromosome_name_last == chromosome_name_internal &
                    ensembl_transcript_id_last != ensembl_transcript_id_internal &
                    ((exon_chrom_start_last <= exon_chrom_end_internal & exon_chrom_end_last >= exon_chrom_start_internal) |
                       (exon_chrom_start_internal <= exon_chrom_end_last & exon_chrom_end_internal >= exon_chrom_start_last))) %>%
    dplyr::group_by(ensembl_transcript_id_last, ensembl_transcript_id_internal) %>%
    dplyr::filter(n() == 1) %>%
    dplyr::ungroup() %>% filter(exon_label_internal == "internal")

  # Select relevant columns to display the hybrid exons
  hybrid_last_exons <- hybrid_last_exons %>% dplyr::select(ensembl_gene_id, ensembl_transcript_id_last, exon_chrom_start_last, exon_chrom_end_last, ensembl_exon_id_last,
                                                           ensembl_transcript_id_internal, exon_chrom_start_internal, exon_chrom_end_internal, ensembl_exon_id_internal)

  ptg_conversion <- biomart_data$ptg_init[biomart_data$ptg_init$chromosome_name %in% c(seq_len(23), "X", "Y"),]

  ptg_conversion$chromosome_name <- NULL

  ptg_conversion$external_gene_name[is.na(ptg_conversion$external_gene_name)] <- ptg_conversion$ensembl_gene_id[is.na(ptg_conversion$external_gene_name)]
  ptg_conversion$external_transcript_name[is.na(ptg_conversion$external_transcript_name)] <- ptg_conversion$ensembl_gene_id[is.na(ptg_conversion$external_transcript_name)]
  ptg_conversion$ensembl_peptide_id[is.na(ptg_conversion$ensembl_peptide_id)] <- "none"

  tgp_biomart <- data.frame(
    gene_id = ptg_conversion$ensembl_gene_id,
    transcript_id = ptg_conversion$ensembl_transcript_id,
    gene_name = ptg_conversion$external_gene_name,
    transcript_name = ptg_conversion$external_transcript_name,
    protein_id = ptg_conversion$ensembl_peptide_id)

  gtf <- data.frame(
    geneID = exon_data_exon_label$ensembl_gene_id,
    transcriptID = exon_data_exon_label$ensembl_transcript_id,
    geneName = exon_data_exon_label$external_gene_name,
    exonID = exon_data_exon_label$ensembl_exon_id,
    classification = exon_data_exon_label$exon_label,
    rownum = seq_len(nrow(exon_data_exon_label)),
    chr = paste0("chr", exon_data_exon_label$chromosome_name),
    start = exon_data_exon_label$exon_chrom_start,
    stop = exon_data_exon_label$exon_chrom_end,
    strand = exon_data_exon_label$strand,
    type = rep("exon", nrow(exon_data_exon_label)),
    gpc = exon_data_exon_label$gene_biotype,
    tpc = exon_data_exon_label$transcript_biotype
  )

  gtf <- gtf[gtf$gpc == "protein_coding",]

  merge_transcript_exon <- gtf[gtf$classification == "first",c("transcriptID", "rownum")]
  merge_transcript_exon <- merge_transcript_exon[!duplicated(merge_transcript_exon),]
  transcript_gtf <- dplyr::left_join(biomart_data$transcript_data, merge_transcript_exon, by = dplyr::join_by('ensembl_transcript_id' == 'transcriptID'))
  transcript_gtf <- transcript_gtf[transcript_gtf$ensembl_transcript_id %in% gtf$transcriptID,]
  transcript_gtf$chromosome_name <- paste0('chr', transcript_gtf$chromosome_name)
  colnames(transcript_gtf) <- c('transcriptID', 'chr', 'strand', 'start', 'stop', 'rownum')
  transcript_gtf <- transcript_gtf %>% dplyr::mutate(strand = ifelse(strand == -1, "-", "+"))
  gtf <- gtf %>% dplyr::mutate(strand = ifelse(strand == -1, "-", "+"))

  hleList <- unique(c(paste0(hybrid_last_exons$ensembl_transcript_id_last, ';', hybrid_last_exons$ensembl_transcript_id_internal),
                      paste0(hybrid_last_exons$ensembl_transcript_id_internal, ';', hybrid_last_exons$ensembl_transcript_id_last)))
  hfeList <- unique(c(paste0(hybrid_first_exons$ensembl_transcript_id_first, ';', hybrid_first_exons$ensembl_transcript_id_internal),
                      paste0(hybrid_first_exons$ensembl_transcript_id_internal, ';', hybrid_first_exons$ensembl_transcript_id_first)))

  gtf <- gtf %>% dplyr::group_by(transcriptID) %>%
    dplyr::mutate(
      adjusted_start = ifelse(strand == "+", start, -start)
    ) %>%
    dplyr::arrange(transcriptID, adjusted_start) %>%
    dplyr::select(-adjusted_start) %>%
    ungroup()


  # gtf <- gtf[gtf$transcriptID %in% use_transcripts,]
  # tgp_biomart <- tgp_biomart[tgp_biomart$transcript_id %in% use_transcripts,]
  # transcript_gtf <- transcript_gtf[transcript_gtf$ensembl_transcript_id %in% use_transcripts,]

  return(list(gtf = gtf,
              transcript_gtf = transcript_gtf,
              hybrid_last_extract = hybrid_last_exons,
              hybrid_first_extract = hybrid_first_exons,
              hybrid_first_extract_transcripts = hfeList,
              hybrid_last_extract_transcripts = hleList,
              tgp_biomart = tgp_biomart))
}


#' import / load + save biomaRt data
#' @param save_location path to save output
#' @importFrom biomaRt useEnsembl getBM
#' @return various biomaRt results for use downstream in pipeline
#'
#' @examples
#' \donttest{
#' library(biomaRt)
#' new_config <- httr::config(ssl_verifypeer = FALSE)
#' httr::set_config(new_config, override = FALSE)
#' ensembl <- biomaRt::useEnsembl(biomart = "ensembl",
#'                                dataset = "hsapiens_gene_ensembl",
#'                                mirror = 'useast')
#' pdir <- system.file("extdata", package = "SpliceImpactR")
#' dataDirectory <- paste0(pdir, "/")
#' biomart_data <- setupBiomart(dataDirectory)
#' }
#' @export

setupBiomart <- function(save_location) {
  ##setup_gtf
  if (!(file.exists(paste0(save_location, 'setup_gtf_exon_data.csv')))) {
    attributes <- c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand",
                    "ensembl_gene_id", "ensembl_transcript_id", "gene_biotype",
                    "transcript_biotype", "external_gene_name", "rank", "ensembl_exon_id")
    setup_gtf_exon_data <- biomaRt::getBM(attributes = attributes, mart = ensembl)
    write_csv(setup_gtf_exon_data, paste0(save_location, 'setup_gtf_exon_data.csv'))
  } else {
    setup_gtf_exon_data <- read_csv(paste0(save_location, 'setup_gtf_exon_data.csv'))
  }

  if (!(file.exists(paste0(save_location, 'ptg_init.csv')))) {
    attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name",
                    "external_transcript_name", "ensembl_peptide_id", "chromosome_name")
    ptg_init <- biomaRt::getBM(attributes = attributes, mart = ensembl)
    write_csv(ptg_init, paste0(save_location, 'ptg_init.csv'))
  } else {
    ptg_init <- read_csv(paste0(save_location, 'ptg_init.csv'))
  }

  if (!(file.exists(paste0(save_location, 'transcript_data.csv')))) {
    attributes <- c("ensembl_transcript_id", "chromosome_name",
                    "strand", "transcript_start", "transcript_end")
    transcript_data <- biomaRt::getBM(attributes = attributes, mart = ensembl,
                                      values = c(seq_len(23), "X", "Y"), filters = 'chromosome_name')
    write_csv(transcript_data, paste0(save_location, 'transcript_data.csv'))
  } else {
    transcript_data <- read_csv(paste0(save_location, 'transcript_data.csv'))
  }

  ##getTTI
  if (!(file.exists(paste0(save_location, 'pfam_data.csv')))) {
    attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "pfam",
                    "transcript_biotype")
    pfam_data <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(seq_len(23), "X", "Y"), "protein_coding"),
                                filters = c('chromosome_name', "transcript_biotype"))
    write_csv(pfam_data, paste0(save_location, 'pfam_data.csv'))
  } else {
    pfam_data <- read_csv(paste0(save_location, 'pfam_data.csv'))
  }

  ##getPFAM
  if (!(file.exists(paste0(save_location, 'ip.csv')))) {
    attributes <- c("ensembl_transcript_id", "chromosome_name",
                    "transcript_biotype","interpro_description", "interpro")
    ip <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(seq_len(23), "X", "Y"), "protein_coding"),
                         filters = c('chromosome_name', "transcript_biotype"))
    write_csv(ip, paste0(save_location, 'ip.csv'))
  } else {
    ip <- read_csv(paste0(save_location, 'ip.csv'))
  }

  if (!(file.exists(paste0(save_location, 'pfam_exon_level.csv')))) {
    attributes <- c("ensembl_transcript_id", "ensembl_exon_id", "pfam", 'chromosome_name', 'transcript_biotype', 'pfam_start', 'pfam_end')
    b1 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(1, 2), "protein_coding"),
                         filters = c('chromosome_name', "transcript_biotype"))
    b2 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(3, 4), "protein_coding"),
                         filters = c('chromosome_name', "transcript_biotype"))
    b3 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(5, 6, 7), "protein_coding"),
                         filters = c('chromosome_name', "transcript_biotype"))
    b4 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(8, 9, 10), "protein_coding"),
                         filters = c('chromosome_name', "transcript_biotype"))
    b5 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(11, 12, 13), "protein_coding"),
                         filters = c('chromosome_name', "transcript_biotype"))
    b6 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(14, 15, 16), "protein_coding"),
                         filters = c('chromosome_name', "transcript_biotype"))
    b7 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(17, 18, 19, 20), "protein_coding"),
                         filters = c('chromosome_name', "transcript_biotype"))
    b8 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c(21, 22, 23), "protein_coding"),
                         filters = c('chromosome_name', "transcript_biotype"))
    b9 <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list(c("X", "Y"), "protein_coding"),
                         filters = c('chromosome_name', "transcript_biotype"))
    pfam_exon_level <- do.call(rbind, list(b1, b2, b3, b4, b5, b6, b7, b8, b9))
    write_csv(pfam_exon_level, paste0(save_location, 'pfam_exon_level.csv'))
  } else {
    pfam_exon_level <- read_csv(paste0(save_location, 'pfam_exon_level.csv'))
  }

  if (!(file.exists(paste0(save_location, 'code_regions.csv')))) {
    attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id", "cds_start", "cds_end")
    code_regions <- biomaRt::getBM(attributes = attributes, mart = ensembl, values = list("protein_coding"), filters = c("transcript_biotype"))
    write_csv(code_regions, paste0(save_location, 'code_regions.csv'))
  } else {
    code_regions <- read_csv(paste0(save_location, 'code_regions.csv'))
  }

  ##frameShiftDetector
  if (!(file.exists(paste0(save_location, 'fsd_exon_data.csv')))) {
    attributes <- c('ensembl_transcript_id', "exon_chrom_start", "exon_chrom_end", "ensembl_exon_id",'cds_start', 'cds_end', 'phase', 'end_phase',
                    'genomic_coding_start', 'genomic_coding_end', 'strand')
    fsd_exon_data <- biomaRt::getBM(attributes = attributes, mart = ensembl)
    write_csv(fsd_exon_data, paste0(save_location, 'fsd_exon_data.csv'))
  } else {
    fsd_exon_data <- read_csv(paste0(save_location, 'fsd_exon_data.csv'))
  }

  return(list(
    setup_gtf_exon_data = setup_gtf_exon_data,
    ptg_init = ptg_init,
    transcript_data = transcript_data,
    pfam_data = pfam_data,
    ip = ip,
    pfam_exon_level = pfam_exon_level,
    code_regions = code_regions,
    fsd_exon_data = fsd_exon_data
  ))
}
