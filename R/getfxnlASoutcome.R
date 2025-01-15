#' Full run of the pipeline for a single alternative splicing event type
#'
#' @param test_group the paths of one phenotype
#' @param control_group the paths of one phenotype
#' @param data_df dataframe with sample paths, type, and phenotype
#' @param mOverlap overlap to identify a match to annotation
#' @param exon_type type of exon being investigated
#' @param output_location location to make  directory
#' @param cutoff cutoff for significance of differential inclusion
#' @param bg whether bg preran or needs making
#' @param plotAlignments whether to output paired alignments
#' @param gtf output from setup_gtf
#' @param tti_location location of previuously made tti or ""
#' @param full_pipe if TRUE, doesn't save output to R, only writes to output_location
#' @param translations from getTranslations
#' @param transcripts from getTranscripts
#' @param max_zero_prop max prop of samples that can be 0
#' @param min_prop_samples min prop of samples from each phenotype required to show a specific event
#' @return nothing or all output from pipeline
#' @export
getfxnlASoutcome <- function(output_location,
                             test_group,control_group,data_df,
                             exon_type, cutoff = .1, outlier_handle = "Inf",
                             cores = 1,
                             tti_location = "", full_pipe = TRUE,
                             bg = NA, mOverlap = .05, gtf,
                             plotAlignments = FALSE,
                             transcripts, translations,
                             biomart_data,
                             max_zero_prop = .5,
                             min_prop_samples = .5) {
  system(paste0("mkdir ",  output_location))
  pdir <- system.file(package="SpliceImpactR")

  if (exon_type %in% c("AFE", "HFE")) {
      diHIT <- differential_inclusion_HITindex(test_names = test_group, control_names = control_group, et = "AFE",
                                               outlier_threshold = outlier_handle, minReads = 10,
                                               min_prop_samples)

      diAS <- diHIT[diHIT$type == "AFE",]
  } else if (exon_type %in% c("ALE", "HLE")) {
      diHIT <- differential_inclusion_HITindex(test_names = test_group, control_names = control_group, et = "ALE",
                                               outlier_threshold = outlier_handle, minReads = 10,
                                               min_prop_samples)

      diAS <- diHIT[diHIT$type == "ALE",]
    } else {
    diAS <- differential_inclusion_rMATS(test_names = test_group, control_names = control_group,
                                         et = exon_type, outlier_threshold = outlier_handle,
                                         minReads = 10,
                                         min_prop_samples)
  }


  fg <- getForeground(input=diAS,
                     test_names = test_group,
                     control_names = control_group,
                     thresh = cutoff,
                     fdr=.05,
                     mOverlap=mOverlap,
                     cores=cores,
                     exon_type=exon_type,
                     pdir=pdir,
                     output_location=output_location, gtf=gtf,
                     max_zero_prop, min_prop_samples, translations)

  if (!(exon_type %in% c("HFE", "HLE"))) {
    initial_comparison <- getOverviewComparison(data_df, exon_type, output_location)
  }

  pfg <- getPaired(foreground = fg$proBed, et = exon_type, nucleotides = transcripts,
                   output_location = output_location, newGTF = gtf, saveAlignments = plotAlignments,
                   exon_data = biomart_data$fsd_exon_data)

  if (nrow(pfg$paired_proBed) > 1) {
    length_comparison <- getLengthComparison(data_df, pfg$paired_proBed, output_location)
  }


  bg_param <- suppressWarnings(is.na(bg))
  if (bg_param) {
    bg_input <- gsub("[^/]*$", "", c(control_group, test_group))
    bg <- getBackground(input=bg_input,
                        mOverlap = mOverlap,
                        cores = cores,
                        exon_type = exon_type,
                        pdir = pdir,
                        output_location = output_location, gtf=gtf, translations)
  }


  #####
  pfam <- getPfam(foreground = fg, background = bg, pdir = pdir,
                  cores = cores, output_location = output_location, biomart_data)


  gD <- getDomainData(fg = fg, bg = bg, pfg=pfg, cores = cores, pfam = pfam, output_location = output_location,
                fdr_use = .05, min_sample_success = 2, engine = "Pfam", topViz = 15)

  if (nrow(pfg$paired_proBed) > 1) {
    tti <- getTTI(paired_foreground = pfg$paired_proBed, background = bg$proBed,
                  pdir = pdir,
                  steps=1,
                  max_vertices_for_viz = 2000,
                  fdr = .05,
                  plot_bool = TRUE,
                  ppidm_class = "Gold_Standard",
                  write_igraphs = TRUE,
                  output_location = output_location,
                  tti_location = tti_location, tgp = gtf$tgp_biomart)
  } else {
    tti <- NA
  }

  if (full_pipe) {
    return(NA)
  } else {
    return(list(diAS = diAS,
                fg = fg,
                pfg = pfg,
                bg = bg,
                pfam = pfam,
                gD = gD,
                tti = tti))
  }

}
