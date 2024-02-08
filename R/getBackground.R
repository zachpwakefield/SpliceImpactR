getBackground <- function(input, mOverlap, cores, nC, nE, exon_type, pdir, output_location) {
  ## extract all first exons and create combined data.frame with gene, location
  files <- paste(input, unlist(lapply(input, function(x) list.files(x)[grep('[.]exon', list.files(x))])), sep = "")

  if (exon_type == "AFE") {
    lim <- c("first")
  } else if (exon_type == "ALE") {
    lim <- c("last")
  } else if (exon_type == "SE") {
    lim <- c("internal")
  }
  first_exons <- unique(unlist(lapply(files, function(x) {
    in_file <- read.delim(x)
    in_file <- in_file[in_file$ID == lim,]
    paste(in_file$gene, ';', in_file$exon, ';',  in_file$strand, sep = "")})))
  redExon <- data.frame(geneR = unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 1)), split = '[.]'), "[[", 1)),
                        chr = unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 1)), split = ":"), "[[", 1)),
                        start = as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 1)), split = ":"), "[[", 2))),
                        stop = as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 2))
                        )
  )

  ## Remove duplicate rows
  redExon <- redExon[!duplicated(redExon),]

  print("exon loaded...")

  ## Use getTranscriptBackground() to extract total and (if background = F) paired transcripts matched to input exons
  matched <- getTranscriptBackground(gtf = gtf, redExon = redExon, ex_type = exon_type, minOverlap = mOverlap, cores = cores)
  print("exons matched, bed-ifying...")

  ## Use bedifyBackground() to extract the total or matched (if background = F) bed file

  bed <- bedifyBackground(matched, outname = output_location, cores = cores)
  print("done bed-ifying...")

  ## extract unqiue transcript names as trans and all trancript names as possT
  trans <- unlist(lapply(strsplit(unique(bed$name), "#"), "[[", 1))
  possT <- unlist(lapply(strsplit(bed$name, "#"), "[[", 1))


  print("Finding annotated proteins...")
  ## Find annotated proteins for transcripts if possible
  protCode <- unlist(parallel::mclapply(trans, mc.cores = cores, function(x) {

    c_trans[which(c_trans == x)[1]+1]

  }))
  protCode[is.na(protCode)]<- "none"

  print("Making output data...")
  proBed <- data.frame(id = unique(bed$name), strand = unlist(lapply(unique(bed$name), function(x) unique(bed$strand[bed$name == x][1]))),
                       prot = protCode) %>%
    tidyr::separate(id, c("transcript", "id"), "#") %>%
    tidyr::separate("id", c("gene", "chr"), ";") %>%
    tidyr::separate('chr', c('chr', 'coords'), ':') %>%
    tidyr::separate('coords', c('start', 'stop'), '-')

  ## Make fasta file with id & strand in first line and protein code
  proFast <- c()
  for (i in 1:length(proBed[,1])) {
    proFast <- c(proFast, paste(">", proBed$transcript[i], "#", proBed$gene[i], ";", proBed$chr[i], ":", proBed$start[i], "-", proBed$stop[i], ";", proBed$strand[i], sep = ""),
                 proBed$prot[i])
  }


  write_csv(proBed, paste0(output_location, "bgoutBed.csv"))
  write_lines(proFast, paste0(output_location, "bgoutFast.fa"))
  write_csv(matched,  paste0(output_location, "bgmatched.csv"))
  write_csv(bed,  paste0(output_location, "bgbed.csv"))
  return(list(matched = matched,
              bed = bed,
              proBed = proBed,
              proFast = proFast))
}
