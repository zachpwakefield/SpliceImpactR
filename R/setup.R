# get_c_trans
# input: directory containing f_protein_code_from_gencodev43.txt, default should be package directory or: system.file(package="domainEnrichment")
# output: f_protein_code_from_gencodev43.txt lines

get_c_trans <- function(translations_location) {
  c_trans_pre <- readr::read_lines(translations_location)
  transcript_title <- c(grep(">", c_trans_pre), (length(c_trans_pre)+1))
  tts <- unlist(lapply(strsplit(unlist(lapply(strsplit(c_trans_pre[transcript_title[1:(length(transcript_title)-1)]], split = "[|]"), "[[", 2)), split = "[.]"), "[[", 1))
  c_trans <- unlist(lapply(1:(length(transcript_title)-1), function(x) {
    c(tts[x], paste(c_trans_pre[(transcript_title[x]+1):(transcript_title[x+1]-1)], collapse = ""))
  }))
  return(c_trans)
}

get_c_nucs <- function(transcripts_location) {
  nc <- readr::read_lines("/projectnb2/evolution/zwakefield/Annotations/hg38_gencode/gencode.v45.pc_transcripts.fa")
  transcript_title <- c(grep(">", nc), (length(nc)+1))
  splitF <- strsplit(nc[(transcript_title[1:(length(transcript_title)-1)])], split = "[|]")
  cdsL <- unlist(lapply(splitF, function(x) grep("CDS:", x)))
  tts <- gsub(">", "", unlist(lapply(strsplit(unlist(lapply(strsplit(nc[transcript_title[1:(length(transcript_title)-1)]], split = "[|]"), "[[", 1)), split = "[.]"), "[[", 1)))
  cds_locs <- unlist(lapply(1:(length(transcript_title)-1), function(x) {
    strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(nc[transcript_title[x]], split = "[|]"), "[[", cdsL[x])), split = ":"), "[[", 2)), split = "[-]")
  }), recursive = F)
  c_trans <- unlist(lapply(1:(length(transcript_title)-1), function(x) {
    c(tts[x], substr(paste(nc[(transcript_title[x]+1):(transcript_title[x+1]-1)], collapse = ""),
                     as.numeric(cds_locs[[x]][1]), as.numeric(cds_locs[[x]][2]) ))
  }))
  transDF <- data.frame(transcriptID = c_trans[seq(1, length(c_trans), by = 2)],
                        code = c_trans[seq(2, length(c_trans), by = 2)])
  # transDF$length <- nchar(transDF$code)
  return(list(c_trans = c_trans, transDF = transDF))
}

# get_gtf
# input: directory containing all the gtf files, default should be package directory or: system.file(package="domainEnrichment")
# output: cg csv

get_gtf <- function(location, cores) {
  cg <- setup_gtf(location, cores = cores)
  return(cg)
}
