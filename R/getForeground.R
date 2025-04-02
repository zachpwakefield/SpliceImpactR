#' get Foreground from differentially included transcripts in samples
#'
#' @param input diExon info from the differentially included functions
#' @param mOverlap overlap to identify a match to annotation
#' @param exon_type placeholder for other functions
#' @param thresh delta psi threshold to filter
#' @param fdr adj p to filter for
#' @param gtf full output from setup_gtf
#' @param max_zero_prop max number of zeros in to count an exon
#' @param min_prop_samples min proportion of samples needed to not be counted as outliers
#' @param output_location location to make background directory
#' @param translations from getTranslations
#' @param test_names path to test files without .{event}PSI
#' @param control_names path to control files without .{event}PSI
#' @param cores number of cores to assign
#' @return matched : matched transcripts dataframe, bed : bed file of the matched transcripts
#' proBed : output for further functions with protein code and protein info,
#' proFast : fasta file of proteins identified in proBed
#' @importFrom dplyr arrange first left_join group_by summarise
#' @importFrom tidyr separate
#' @importFrom readr read_csv read_lines
#'
#'
#' @examples
#' pdir <- system.file("extdata", package="SpliceImpactR")
#' dataDirectory <- paste0(pdir, "/")
#' test_group <- paste0(dataDirectory, "rawData/", c("test1","test2", "test3"))
#' control_group <- paste0(dataDirectory, "rawData/", c("control1", "control2", "control3"))
#'
#' result <- differential_inclusion_HITindex(test_names = test_group,
#'                                           control_names = control_group,
#'                                           et = "AFE",
#'                                           outlier_threshold = "Inf",
#'                                           minReads = 10,
#'                                           min_prop_samples = 0,
#'                                           chosen_method = "qbGLM"
#'                                           )
#' gtf_sample <- list(gtf = readr::read_csv(paste0(dataDirectory, "gtf_limited.csv")),
#'             transcript_gtf = readr::read_csv(paste0(dataDirectory, "transcript_gtf_limited.csv")))
#' translations_sample <- readr::read_lines(paste0(dataDirectory, "translations_limited.csv"))
#'
#' fg <- getForeground(input = result,
#'                             test_names = test_group,
#'                             control_names = control_group,
#'                             thresh = .1,
#'                             fdr = .05,
#'                             mOverlap = .1,
#'                             exon_type = "AFE",
#'                             output_location = NULL,
#'                             cores = 1,
#'                             gtf = gtf_sample,
#'                             max_zero_prop = 1,
#'                             min_prop_samples = 0,
#'                             translations = translations_sample)
#'
#'
#' @export
getForeground <- function(input, test_names, control_names, thresh = .1, fdr = 0.05,
                          mOverlap,exon_type,
                          output_location = NULL,
                          cores = 1,
                          gtf,
                          max_zero_prop = .5,
                          min_prop_samples = .5,
                          translations) {

  ## If using foreground set, read in diExon file and extract differentially included exons using diff_info()
  df <- input
  df.l <- diColor(df, color_thresh = .2)

  ## Filter df.l for paired analysis
  qdf.l <- qualityFilter(df.l, nT = length(test_names),
                         nC = length(control_names),
                         min_prop_samples = min_prop_samples,
                         max_zero_prop = max_zero_prop)
  bdf.l <- significanceFilter(qdf.l, fdr = fdr, d.psi = thresh)

  ## Make volcano plot with make_lfcPlot()
  lfcPlot <- make_dPsiPlot(qdf.l, gtf$gtf, thresh = thresh)

  ## Make data.frame with gene, location of each exon on total foreground
  redExon <- data.frame(geneR = unlist(lapply(strsplit(bdf.l$gene, split = "[.]"), "[[", 1)),
                        chr = unlist(lapply(strsplit(bdf.l$exon, split = ":"), "[[", 1)),
                        start = as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(bdf.l$exon, split = ":"),
                                                                  "[[", 2)), split = "[-]"), "[[", 1))),
                        stop = as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(bdf.l$exon, split = ":"),
                                                                 "[[", 2)), split = "[-]"), "[[", 2))),
                        delta.psi = bdf.l$delta.psi,
                        p.adj = bdf.l$p.adj,
                        add_inf = bdf.l$add_inf
  )


  ## Use getTranscriptForeground() to extract total and (if background = F) paired transcripts matched to input exons
  matched <- getTranscriptForeground(gtf$gtf, redExon, ex_type = exon_type, minOverlap = mOverlap, cores, gtf$transcript_gtf)

  ## Use bedifyForeground() to extract the  bed file
  bed <- bedifyForeground(matched, outname = output_location, cores = cores, gtf=gtf$gtf)

  ## extract uniqiue transcript names as trans and all trancript names as possT
  trans <- str_extract(unique(bed$name), "^[^#]*")

  # For all names
  possT <- str_extract(bed$name, "^[^#]*")

  ## Find annotated proteins for transcripts if possible
  protCode <- translations[match(trans, translations) + 1]
  protCode[is.na(protCode)]<- "none"

  merger <- data.frame(name = unique(bed$name),
                       prot = protCode)

  ## Make dataframe proBed for output of matched transcripts with protein code
  bed_summary <- bed %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(strand = dplyr::first(strand),
                     delta.psi = dplyr::first(delta.psi),
                     p.adj = dplyr::first(p.adj),
                     add_inf = dplyr::first(add_inf),
                     exonID = dplyr::first(exonID),
                     .groups = 'drop')


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

  if (!is.null(output_location)) {
    system(paste0("mkdir ", output_location, "Foreground/"))
    readr::write_csv(proBed, paste0(output_location, "Foreground/", "fgoutBed.csv"))
    readr::write_lines(proFast, paste0(output_location, "Foreground/", "fgoutFast.fa"))
    readr::write_csv(matched,  paste0(output_location, "Foreground/", "fgmatched.csv"))
    readr::write_csv(df.l,  paste0(output_location, "Foreground/", "fglfc.csv"))
    readr::write_csv(bed,  paste0(output_location, "Foreground/", "fgexonBed.csv"))

    pdf(paste0(output_location, "Foreground/", "delta_psi_plot.pdf"))
    print(lfcPlot)
    dev.off()
  }

  return(list(proBed = proBed,
              proFast = proFast,
              protCode = protCode,
              matched = matched))
}
