getForeground <- function(input, test_names, control_names, thresh, fdr, mOverlap, nC, nE, exon_type, pdir, output_location, cores) {

  ## If using foreground set, read in diExon file and extract differentially included exons using lfc()
  df <- input
  df.l <- lfc(df, numCont = nC, numExp = nE, exon_type = exon_type,
              cores = cores, test_names = test_names, control_names = control_names)

  ## Filter df.l for paired analysis
  bdf.l <- df.l[abs(df.l$delta.psi) >= thresh & df.l$p.adj <= fdr,]

  ## Make volcano plot with make_lfcPlot()
  lfcPlot <- make_lfcPlot(df.l, pdir = pdir)

  ## Make data.frame with gene, location of each exon on total foreground
  redExon <- data.frame(geneR = unlist(lapply(strsplit(bdf.l$gene, split = "[.]"), "[[", 1)),
                        chr = sapply(strsplit(bdf.l$exon, split = ":"), "[[", 1),
                        start = as.numeric(sapply(strsplit(sapply(strsplit(bdf.l$exon, split = ":"), "[[", 2), split = "[-]"), "[[", 1)),
                        stop = as.numeric(sapply(strsplit(sapply(strsplit(bdf.l$exon, split = ":"), "[[", 2), split = "[-]"), "[[", 2)),
                        delta.psi = bdf.l$delta.psi,
                        p.adj = bdf.l$p.adj,
                        add_inf = bdf.l$add_inf
  )


  print("exon loaded...")

  ## Use getTranscriptForeground() to extract total and (if background = F) paired transcripts matched to input exons
  matched <- getTranscriptForeground(gtf = gtf, redExon = redExon, ex_type = exon_type, minOverlap = mOverlap, cores = cores)
  print("exons matched, bed-ifying...")

  ## Use bedifyForeground() to extract the  bed file

  bed <- bedifyForeground(matched, outname = output_location, cores = cores)
  print("done bed-ifying...")

  ## extract uniqiue transcript names as trans and all trancript names as possT
  trans2 <- str_extract(unique(bed$name), "^[^#]*")

  # For all names
  possT2 <- str_extract(bed$name, "^[^#]*")

  print("Finding annotated proteins...")
  ## Find annotated proteins for transcripts if possible


  protCode <- c_trans[match(trans, c_trans) + 1]
  protCode[is.na(protCode)]<- "none"

  print("Making output data...")
  ## Make dataframe proBed for output of matched transcripts with protein code
  bed_summary <- bed %>%
    group_by(name) %>%
    summarise(strand = first(strand),
              delta.psi = first(delta.psi),
              p.adj = first(p.adj),
              add_inf = first(add_inf),
              .groups = 'drop')

  # Assuming 'protCode' is already defined and in the correct order for the unique 'bed$name'
  # If 'protCode' needs to be matched by 'bed$name', ensure that step is handled accordingly

  # Create the initial 'proBed' data frame without the need for 'lapply' or 'unique'
  proBed <- bed_summary %>%
    mutate(prot = protCode) %>%
    separate(name, c("transcript", "id"), "#") %>%
    separate(id, c("gene", "chr"), ";") %>%
    separate(chr, c("chr", "coords"), ':') %>%
    separate(coords, c("start", "stop"), '-')


  # Create the FASTA headers using vectorized paste function
  fasta_headers <- paste0(">", proBed$transcript, "#", proBed$gene, ";", proBed$chr, ":", proBed$start, "-", proBed$stop, ";", proBed$strand)

  # Assuming 'proBed$prot' contains the protein sequences
  fasta_sequences <- proBed$prot

  # Interleave headers and sequences
  proFast2 <- paste(rbind(fasta_headers, fasta_sequences))

  write_csv(proBed, paste0(output_location, "fgoutBed.csv"))
  write_lines(proFast, paste0(output_location, "fgoutFast.fa"))
  write_csv(matched,  paste0(output_location, "fgmatched.csv"))
  write_csv(df.l,  paste0(output_location, "fglfc.csv"))
  write_csv(bed,  paste0(output_location, "fgexonBed.csv"))


  return(list(proBed = proBed,
              proFast = proFast,
              protCode = protCode,
              matched = matched))
}
