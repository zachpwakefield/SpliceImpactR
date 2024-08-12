#' setting up gtf with exon labels, first, internal, last, hybrid, etc
#' @importFrom dplyr filter group_by mutate ungroup inner_join select case_when
#' @importFrom biomaRt useEnsembl getBM
#' @return annotated gtf and various hybrid info. also gene to transcript to protein naming dataframe
#' @keywords internal

setupAnnotation <- function() {
  ensembl <- biomaRt::useEnsembl(biomart = "ensembl",
                                 dataset = "hsapiens_gene_ensembl")
  attributes <- c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand", "ensembl_gene_id", "ensembl_transcript_id", "gene_biotype",
                  "transcript_biotype", "external_gene_name", "rank", "ensembl_exon_id")
  exon_data <- biomaRt::getBM(attributes = attributes, mart = ensembl)

  exon_data_exon_label <- exon_data %>% dplyr::filter(chromosome_name %in% c(1:23, "X", "Y")) %>%
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
    dplyr::filter(n() == 1) %>%
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



  attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name", "external_transcript_name", "ensembl_peptide_id", "chromosome_name")
  ptg_init <- biomaRt::getBM(attributes = attributes, mart = ensembl, uniqueRows = TRUE)
  ptg_conversion <- ptg_init[ptg_init$chromosome_name %in% c(1:23, "X", "Y"),]
  ptg_conversion$chromosome_name <- NULL
  ptg_conversion$external_gene_name[ptg_conversion$external_gene_name == ""] <- ptg_conversion$ensembl_gene_id[ptg_conversion$external_gene_name == ""]
  ptg_conversion$external_transcript_name[ptg_conversion$external_transcript_name == ""] <- ptg_conversion$ensembl_gene_id[ptg_conversion$external_transcript_name == ""]
  ptg_conversion$ensembl_peptide_id[ptg_conversion$ensembl_peptide_id == ""] <- "none"

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
    rownum = 1:nrow(exon_data_exon_label),
    chr = paste0("chr", exon_data_exon_label$chromosome_name),
    start = exon_data_exon_label$exon_chrom_start,
    stop = exon_data_exon_label$exon_chrom_end,
    strand = exon_data_exon_label$strand,
    type = rep("exon", nrow(exon_data_exon_label)),
    gpc = exon_data_exon_label$gene_biotype,
    tpc = exon_data_exon_label$transcript_biotype
  )
  gtf <- gtf[gtf$gpc == "protein_coding",]


  transcript_attributes <- c("ensembl_transcript_id", "chromosome_name", "strand", "transcript_start", "transcript_end")
  transcript_data <- biomaRt::getBM(attributes = transcript_attributes, mart = ensembl, values = c(1:23, "X", "Y"), filters = 'chromosome_name')
  merge_transcript_exon <- gtf[gtf$classification == "first",c("transcriptID", "rownum")]
  merge_transcript_exon <- merge_transcript_exon[!duplicated(merge_transcript_exon),]
  transcript_gtf <- dplyr::left_join(transcript_data, merge_transcript_exon, by = join_by('ensembl_transcript_id' == 'transcriptID'))
  transcript_gtf <- transcript_gtf[transcript_gtf$ensembl_transcript_id %in% gtf$transcriptID,]
  transcript_gtf$chromosome_name <- paste0('chr', transcript_gtf$chromosome_name)
  colnames(transcript_gtf) <- c('transcriptID', 'chr', 'strand', 'start', 'stop', 'rownum')
  transcript_gtf <- transcript_gtf %>% dplyr::mutate(strand = ifelse(strand == -1, "-", "+"))
  gtf <- gtf %>% dplyr::mutate(strand = ifelse(strand == -1, "-", "+"))



  hleList <- unique(c(paste0(hybrid_last_exons$ensembl_transcript_id_last, ';', hybrid_last_exons$ensembl_transcript_id_internal),
                      paste0(hybrid_last_exons$ensembl_transcript_id_internal, ';', hybrid_last_exons$ensembl_transcript_id_last)))
  hfeList <- unique(c(paste0(hybrid_first_exons$ensembl_transcript_id_first, ';', hybrid_first_exons$ensembl_transcript_id_internal),
                      paste0(hybrid_first_exons$ensembl_transcript_id_internal, ';', hybrid_first_exons$ensembl_transcript_id_first)))


  return(list(gtf = gtf,
              transcript_gtf = transcript_gtf,
              hybrid_last_extract = hybrid_last_exons,
              hybrid_first_extract = hybrid_first_exons,
              hybrid_first_extract_transcripts = hfeList,
              hybrid_last_extract_transcripts = hleList,
              tgp_biomart = tgp_biomart))
}
