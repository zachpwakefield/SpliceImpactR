#' get Foreground differentially included transcripts in samples
#'
#' @param input diExon info from the differentially included functions
#' @param mOverlap overlap to identify a match to annotation
#' @param exon_type placeholder for other functions
#' @param thresh delta psi threshold to filter
#' @param fdr adj p to filter for
#' @param gtf gtf dataframe from setup_gtf
#' @param pdir location of the package
#' @param output_location location to make background directory
#' @return matched : matched transcripts dataframe, bed : bed file of the matched transcripts
#' proBed : output for further functions with protein code and protein info,
#' proFast : fasta file of proteins identified in proBed
#' @importFrom dplyr arrange first left_join group_by summarise
#' @importFrom tidyr separate
#' @export
getForeground <- function(input, test_names, control_names, thresh, fdr,
                          mOverlap,exon_type, pdir,
                          output_location, cores = 1, gtf) {

  ## If using foreground set, read in diExon file and extract differentially included exons using diff_info()
  df <- input
  df.l <- diColor(df, color_thresh = .2)

  ## Filter df.l for paired analysis
  qdf.l <- qualityFilter(df.l, nT = length(test_names), nC = length(control_names))
  bdf.l <- significanceFilter(qdf.l)

  ## Make volcano plot with make_lfcPlot()
  lfcPlot <- make_dPsiPlot(df.l, thresh = thresh, pdir = pdir)

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
    dplyr::group_by(name) %>%
    dplyr::summarise(strand = dplyr::first(strand),
                     delta.psi = dplyr::first(delta.psi),
                     p.adj = dplyr::first(p.adj),
                     add_inf = dplyr::first(add_inf),
                     .groups = 'drop')



  # Assuming 'protCode' is already defined and in the correct order for the unique 'bed$name'
  # If 'protCode' needs to be matched by 'bed$name', ensure that step is handled accordingly

  # Create the initial 'proBed' data frame without the need for 'lapply' or 'unique'
  proBed <- bed_summary %>% dplyr::left_join(merger, by='name') %>%
    tidyr::separate(name, c("transcript", "id"), "#") %>%
    tidyr::separate(id, c("gene", "chr"), ";") %>%
    tidyr::separate(chr, c("chr", "coords"), ':') %>%
    tidyr::separate(coords, c("start", "stop"), '-')


  # Create the FASTA headers using vectorized paste function
  fasta_headers <- paste0(">", proBed$transcript, "#", proBed$gene, ";", proBed$chr, ":", proBed$start, "-", proBed$stop, ";", proBed$strand)

  # Assuming 'proBed$prot' contains the protein sequences
  fasta_sequences <- proBed$prot

  # Interleave headers and sequences
  proFast <- paste(rbind(fasta_headers, fasta_sequences))

  system(paste0("mkdir ", output_location, "Foreground/"))
  write_csv(proBed, paste0(output_location, "Foreground/", "fgoutBed.csv"))
  write_lines(proFast, paste0(output_location, "Foreground/", "fgoutFast.fa"))
  write_csv(matched,  paste0(output_location, "Foreground/", "fgmatched.csv"))
  write_csv(df.l,  paste0(output_location, "Foreground/", "fglfc.csv"))
  write_csv(bed,  paste0(output_location, "Foreground/", "fgexonBed.csv"))

  pdf(paste0(output_location, "Foreground/", "delta_psi_plot.pdf"))
  print(lfcPlot)
  dev.off()

  return(list(proBed = proBed,
              proFast = proFast,
              protCode = protCode,
              matched = matched))
}
