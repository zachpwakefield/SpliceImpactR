#' setting up translations
#' @param translations_location location of translations file
#' @importFrom rtracklayer import
#' @importFrom readr read_lines
#' @importFrom utils download.file
#' @importFrom R.utils gunzip
#' @return modified protein code fasta
#'
#' @examples
#'
#' pdir <- system.file("extdata", package = "SpliceImpactR")
#' dataDirectory <- paste0(pdir, "/")
#' translations <- getTranslations(translations_location = dataDirectory)
#'
#' @export
getTranslations <- function(translations_location) {
  if (!(file.exists(paste0(translations_location, "gencode.v45.pc_translations.fa")))) {
    url <- "ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.pc_translations.fa.gz"
    utils::download.file(url,paste0(translations_location, "gencode.v45.pc_translations.fa.gz"))
    R.utils::gunzip(paste0(translations_location, "gencode.v45.pc_translations.fa.gz"), remove=FALSE, overwrite=TRUE)
  }
  c_trans_pre <- readr::read_lines(paste0(translations_location, "gencode.v45.pc_translations.fa"))
  transcript_title <- c(grep(">", c_trans_pre), (length(c_trans_pre)+1))
  tts <- unlist(lapply(strsplit(unlist(lapply(strsplit(c_trans_pre[transcript_title[seq_len(length(transcript_title)-1)]], split = "[|]"), "[[", 2)), split = "[.]"), "[[", 1))
  c_trans <- unlist(lapply(seq_len(length(transcript_title)-1), function(x) {
    c(tts[x], paste(c_trans_pre[(transcript_title[x]+1):(transcript_title[x+1]-1)], collapse = ""))
  }))
  return(c_trans)
}

#' setting up transcriptions
#' @param transcripts_location location of transcripts file / where to save the transcript file if not already imported
#' @importFrom rtracklayer import
#' @importFrom readr read_lines
#' @importFrom utils download.file
#' @importFrom R.utils gunzip
#' @return modified nucleotide code fasta
#'
#' @examples
#'
#' pdir <- system.file("extdata", package = "SpliceImpactR")
#' dataDirectory <- paste0(pdir, "/")
#' transcripts <- getTranscripts(transcripts_location = dataDirectory)
#'
#' @export
getTranscripts <- function(transcripts_location) {
  if (!(file.exists(paste0(transcripts_location, "gencode.v45.pc_transcripts.fa")))) {
    url <- "ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.pc_transcripts.fa.gz"
    utils::download.file(url,paste0(transcripts_location, "gencode.v45.pc_transcripts.fa.gz"))
    R.utils::gunzip(paste0(transcripts_location, "gencode.v45.pc_transcripts.fa.gz"), remove=FALSE, overwrite=TRUE)
  }
  nc <- readr::read_lines(paste0(transcripts_location, "gencode.v45.pc_transcripts.fa"))
  transcript_title <- c(grep(">", nc), (length(nc)+1))
  splitF <- strsplit(nc[(transcript_title[seq_len(length(transcript_title)-1)])], split = "[|]")
  cdsL <- unlist(lapply(splitF, function(x) grep("CDS:", x)))
  tts <- gsub(">", "", unlist(lapply(strsplit(unlist(lapply(strsplit(nc[transcript_title[seq_len(length(transcript_title)-1)]],
                                                                     split = "[|]"), "[[", 1)), split = "[.]"), "[[", 1)))
  cds_locs <- unlist(lapply(seq_len(length(transcript_title)-1), function(x) {
    strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(nc[transcript_title[x]], split = "[|]"), "[[", cdsL[x])),
                                    split = ":"), "[[", 2)), split = "[-]")
  }), recursive = FALSE)
  c_trans <- unlist(lapply(seq_len(length(transcript_title)-1), function(x) {
    c(tts[x], substr(paste(nc[(transcript_title[x]+1):(transcript_title[x+1]-1)], collapse = ""),
                     as.numeric(cds_locs[[x]][1]), as.numeric(cds_locs[[x]][2]) ))
  }))
  transDF <- data.frame(transcriptID = c_trans[seq(1, length(c_trans), by = 2)],
                        code = c_trans[seq(2, length(c_trans), by = 2)])
  # transDF$length <- nchar(transDF$code)
  return(list(c_trans = c_trans, transDF = transDF))
}

#' wrapper for setup_gtf, extract annotation information
#' @param biomart_data output from setup_biomart
#' @return annotated gtf and others from setup_gtf
#'
#' @examples
#'
#' pdir <- system.file("extdata", package = "SpliceImpactR")
#' dataDirectory <- paste0(pdir, "/")
#'
#' setup_gtf_exon_data <- readr::read_csv(paste0(dataDirectory, "biomart_setup_gtf_exon_data.csv"))
#' ptg_init <- readr::read_csv(paste0(dataDirectory, "biomart_ptg_init.csv"))
#' transcript_data <- readr::read_csv(paste0(dataDirectory, "biomart_transcript_data.csv"))
#'
#' biomart_data_sample <- list(setup_gtf_exon_data = setup_gtf_exon_data,
#'                             ptg_init = ptg_init,
#'                             transcript_data = transcript_data)
#' gtf <- getAnnotation(biomart_data = biomart_data_sample)
#'
#' @export
getAnnotation <- function(biomart_data) {
  cg <- setupAnnotation(biomart_data)
  return(cg)
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    # data/table column-like variables frequently used in your code
    "classification", "geneID", "chr", "jaccard", ".data", "exon", "strand",
    "transcriptID", "gtf_transcripts", "exclusion_rownum", "transcript_starts",
    "nUP", "nDOWN", "nLE", "nFE", "nDiff", "has_nonzero_psi", "valid_group",
    "psi", "valid_reads", "sample_name", "gene", "exon_id", "exon_chrom_start",
    "exon_chrom_end", "delta.psi", "delta.psi.x", "delta.psi.y", "delta_HIT",
    "control_HIT", "test_HIT", "control_average_psi", "test_average_psi",
    "median_sum_IJC_SJC", "zero_count", "p.val", "cooks_d", "inclusion",
    "exclusion", "total", "type", "domain", "domain2", "domains", "status",
    "condition", "condition_sample", "HITindex", "avg_HITindex", "mean_control",
    "mean_HITindex", "gene_exon", "control", "reads", "dataDirectory",
    "status_start", "overlap_status_start", "adjusted_start", "pos_exon_id",
    "neg_exon_id", "pos_start", "neg_start", "group", "size_adjusted",
    "x_start", "x_end", "y", "size", "id", "IncLevel1", "IncLevel2",
    "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IJC",
    "SJC", "IncFormLen", "SkipFormLen", "psi_adjusted", "exon_label",
    "exon_label_internal", "input_id", "chromosome_name",
    "chromosome_name_first", "chromosome_name_internal", "chromosome_name_last",
    "ensembl_gene_id", "ensembl_transcript_id", "ensembl_transcript_id_first",
    "ensembl_transcript_id_internal", "ensembl_transcript_id_last",
    "ensembl_exon_id_first", "ensembl_exon_id_internal", "ensembl_exon_id_last",
    "X1stExonStart_0base", "X1stExonEnd", "X2ndExonStart_0base", "X2ndExonEnd",
    "add_inf", "skip_exon", "median", "p_value", "vals", "val",
    "normAScount", "normASpg", "onePC", "Match", "MedianScore", "PartialMatch",
    "PropRescue", "DomainCount", "FrameShift", "AlignmentScores", "Category",
    "Proportion", "RelativeDomainProp", "event", "present", "exonID", "strands",
    "exon_label", "ensembl", "reg",

    # data.table “special” symbols
    ".N", ".SD",

    # Tidyverse “special” helpers or columns
    ".",  # e.g. from dplyr pipelines or ggplot2
    "starts_with", "contains", "write_csv", "from",

    # Additional or newly discovered missing variables
    "transcript", "X1", "types2", "types", "p.adj", "coords", "dens", "count",
    "fdr", "test", "name", "Median",
    "exon_chrom_start_first", "exon_chrom_end_first",
    "exon_chrom_start_internal", "exon_chrom_end_internal",
    "exon_chrom_start_last", "exon_chrom_end_last"
  ))
}
