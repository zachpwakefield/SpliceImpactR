#' setting up translations
#' @param translations_location location of translations file
#' @importFrom rtracklayer import
#' @importFrom readr read_lines
#' @importFrom utils download.file
#' @importFrom R.utils gunzip
#' @return modified protein code fasta
#' @export
getTranslations <- function(translations_location) {
  if (!(file.exists(paste0(translations_location, "gencode.v45.pc_translations.fa")))) {
    url = "ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.pc_translations.fa.gz"
    utils::download.file(url,paste0(translations_location, "gencode.v45.pc_translations.fa.gz"))
    R.utils::gunzip(paste0(translations_location, "gencode.v45.pc_translations.fa.gz"), remove=FALSE, overwrite=TRUE)
  }
  c_trans_pre <- readr::read_lines(paste0(translations_location, "gencode.v45.pc_translations.fa"))
  transcript_title <- c(grep(">", c_trans_pre), (length(c_trans_pre)+1))
  tts <- unlist(lapply(strsplit(unlist(lapply(strsplit(c_trans_pre[transcript_title[1:(length(transcript_title)-1)]], split = "[|]"), "[[", 2)), split = "[.]"), "[[", 1))
  c_trans <- unlist(lapply(1:(length(transcript_title)-1), function(x) {
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
#' @export
getTranscripts <- function(transcripts_location) {
  if (!(file.exists(paste0(transcripts_location, "gencode.v45.pc_transcripts.fa")))) {
    url = "ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.pc_transcripts.fa.gz"
    utils::download.file(url,paste0(transcripts_location, "gencode.v45.pc_transcripts.fa.gz"))
    R.utils::gunzip(paste0(transcripts_location, "gencode.v45.pc_transcripts.fa.gz"), remove=FALSE, overwrite=TRUE)
  }
  nc <- readr::read_lines(paste0(transcripts_location, "gencode.v45.pc_transcripts.fa"))
  transcript_title <- c(grep(">", nc), (length(nc)+1))
  splitF <- strsplit(nc[(transcript_title[1:(length(transcript_title)-1)])], split = "[|]")
  cdsL <- unlist(lapply(splitF, function(x) grep("CDS:", x)))
  tts <- gsub(">", "", unlist(lapply(strsplit(unlist(lapply(strsplit(nc[transcript_title[1:(length(transcript_title)-1)]],
                                                                     split = "[|]"), "[[", 1)), split = "[.]"), "[[", 1)))
  cds_locs <- unlist(lapply(1:(length(transcript_title)-1), function(x) {
    strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(nc[transcript_title[x]], split = "[|]"), "[[", cdsL[x])),
                                    split = ":"), "[[", 2)), split = "[-]")
  }), recursive = FALSE)
  c_trans <- unlist(lapply(1:(length(transcript_title)-1), function(x) {
    c(tts[x], substr(paste(nc[(transcript_title[x]+1):(transcript_title[x+1]-1)], collapse = ""),
                     as.numeric(cds_locs[[x]][1]), as.numeric(cds_locs[[x]][2]) ))
  }))
  transDF <- data.frame(transcriptID = c_trans[seq(1, length(c_trans), by = 2)],
                        code = c_trans[seq(2, length(c_trans), by = 2)])
  # transDF$length <- nchar(transDF$code)
  return(list(c_trans = c_trans, transDF = transDF))
}

#' wrapper for setup_gtf
#' @return annotated gtf and others from setup_gtf
#' @export
getAnnotation <- function() {
  cg <- setupAnnotation()
  return(cg)
}
