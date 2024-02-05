getForeground <- function(input, thresh, fdr, mOverlap, nC, nE, exon_type, pdir, output_location, cores) {

  ## If using foreground set, read in diExon file and extract differentially included exons using lfc()
  df <- input
  df.l <- lfc(df, numCont = nC, numExp = nE, exon_type = exon_type, cores = cores)

  ## Filter df.l for paired analysis
  bdf.l <- df.l[abs(df.l$delta.psi) >= thresh & df.l$p.adj <= fdr,]

  ## Make volcano plot with make_lfcPlot()
  lfcPlot <- make_lfcPlot(df.l)

  ## Make data.frame with gene, location of each exon on total foreground
  redExon <- data.frame(geneR = unlist(lapply(strsplit(bdf.l$gene, split = "[.]"), "[[", 1)),
                        chr = sapply(strsplit(bdf.l$exon, split = ":"), "[[", 1),
                        start = as.numeric(sapply(strsplit(sapply(strsplit(bdf.l$exon, split = ":"), "[[", 2), split = "[-]"), "[[", 1)),
                        stop = as.numeric(sapply(strsplit(sapply(strsplit(bdf.l$exon, split = ":"), "[[", 2), split = "[-]"), "[[", 2)),
                        delta.psi = bdf.l$delta.psi,
                        p.adj = bdf.l$p.adj
  )


  print("exon loaded...")

  ## Use getTranscript() to extract total and (if background = F) paired transcripts matched to input exons
  matched <- getTranscript2(gtf = gtf, redExon = redExon, ex_type = exon_type, minOverlap = mOverlap, cores = cores)
  print("exons matched, bed-ifying...")

  ## Use bedify() to extract the  bed file

  bed <- bedify2(matched, outname = output_location, cores = cores)
  print("done bed-ifying...")

  ## extract unqiue transcript names as trans and all trancript names as possT
  trans <- unlist(lapply(strsplit(unique(bed$name), "#"), "[[", 1))
  possT <- unlist(lapply(strsplit(bed$name, "#"), "[[", 1))


  ## Find annotated proteins for transcripts if possible
  protCode <- unlist(parallel::mclapply(trans, mc.cores = 8, function(x) {
    rc <- c_trans[which(c_trans == x)+1]
    if (length(rc) > 0) {
      rc[1]
    } else {"none"}
  }))

  ## Make dataframe proBed for output of matched transcripts with protein code
  proBed <- data.frame(id = unique(bed$name),
                       strand = unlist(lapply(unique(bed$name), function(x) unique(bed$strand[bed$name == x][1]))),
                       prot = protCode,
                       delta.psi = unlist(lapply(unique(bed$name), function(x) unique(bed$delta.psi[bed$name == x][1]))),
                       p.adj = unlist(lapply(unique(bed$name), function(x) unique(bed$p.adj[bed$name == x][1])))) %>%
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
