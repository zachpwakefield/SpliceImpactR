getBackground <- function(input, mOverlap, cores, nC, nE, exon_type, pdir, output_location) {
  ## extract all first exons and create combined data.frame with gene, location
  files <- paste(input, unlist(lapply(input, function(x) list.files(x)[grep('[.]exon', list.files(x))])), sep = "")

  if (exon_type == "AFE" | exon_type == "HFE") {
    lim <- c("first")
  } else if (exon_type == "ALE" | exon_type == "HLE") {
    lim <- c("last")
  } else if (!(exon_type %in% c("AFE", "ALE"))) {
    lim <- c("internal")
  }
  first_exons <- unique(unlist(lapply(files, function(x) {
    in_file <- read.delim(x)
    in_file <- in_file[in_file$ID == lim,]
    paste(in_file$gene, ';', in_file$exon, ';',  in_file$strand, sep = "")})))
  first_exons <- first_exons[grepl('[-]', first_exons) & grepl(';', first_exons)]
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

  ## extract unique transcript names as trans and all transcript names as possT
  trans <- str_extract(unique(bed$name), "^[^#]*")

  # For all names
  possT <- str_extract(bed$name, "^[^#]*")


  print("Finding annotated proteins...")
  ## Find annotated proteins for transcripts if possible
  protCode <- c_trans[match(trans, c_trans) + 1]
  protCode[is.na(protCode)]<- "none"

  merger <- data.frame(name = unique(bed$name),
                       prot = protCode)

  print("Making output data...")
  ## Make dataframe proBed for output of matched transcripts with protein code
  bed_summary <- bed %>%
    dplyr::group_by(.data$name) %>%
    dplyr::summarise(strand = dplyr::first(.data$strand),
                     .groups = 'drop')



  # Assuming 'protCode' is already defined and in the correct order for the unique 'bed$name'
  # If 'protCode' needs to be matched by 'bed$name', ensure that step is handled accordingly

  # Create the initial 'proBed' data frame without the need for 'lapply' or 'unique'
  proBed <- bed_summary %>% dplyr::left_join(merger, by='name') %>%
    tidyr::separate(.data$name, c("transcript", "id"), "#") %>%
    tidyr::separate(.data$id, c("gene", "chr"), ";") %>%
    tidyr::separate(.data$chr, c("chr", "coords"), ':') %>%
    tidyr::separate(.data$coords, c("start", "stop"), '-')

  # Create the FASTA headers using vectorized paste function
  fasta_headers <- paste0(">", proBed$transcript, "#", proBed$gene, ";", proBed$chr, ":", proBed$start, "-", proBed$stop, ";", proBed$strand)

  # Assuming 'proBed$prot' contains the protein sequences
  fasta_sequences <- proBed$prot

  # Interleave headers and sequences
  proFast <- paste(rbind(fasta_headers, fasta_sequences))

  system(paste0("mkdir ", output_location, "Background/"))
  write_csv(proBed, paste0(output_location, "Background/", "bgoutBed.csv"))
  write_lines(proFast, paste0(output_location, "Background/", "bgoutFast.fa"))
  write_csv(matched,  paste0(output_location, "Background/", "bgmatched.csv"))
  write_csv(bed,  paste0(output_location, "Background/", "bgbed.csv"))

  return(list(matched = matched,
              bed = bed,
              proBed = proBed,
              proFast = proFast))
}
