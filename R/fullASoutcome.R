#' Extract differentially included AFE, ALE from HIT Index data
#'
#' @param as_types a vector of the type of alternative splicing events that are being investigated
#' @param output_directory the path (ending with '/') that the output is desired
#' @param data_directory the path to the data, in the format specified by organizeSamples
#' @param outlier_handle value to threshold the outlier detection
#' @param cutoff diInclusion cutoff to use to identify significance
#' @param cutoff diInclusion cutoff to use to identify significance
#' @param plotAlignments whether to plot alignments from getPaired/matchAlignType
#' @param bg_pre if bg was made earlier, param to give premade bg
#' @param gtf output from setup_gtf
#' @param tti_location location of previously made transcript-transcript interactions network
#' @param mOverlap minimum overlap to call an exon as matched to annotation
#' @param translations from getTranslations
#' @param transcripts from getTranscripts
#' @return nothing in R, output to the output_directory
#' @export
fullASoutcome <- function(as_types = c("AFE", "ALE", "HFE", "HLE", "SE", "MXE", "RI", "A5SS", "A3SS"),
                          output_directory, data_directory,
                          data_df, outlier_handle,
                          cutoff = .1, cores = 1, bg_pre = NA,
                          tti_location = "/projectnb/evolution/zwakefield/allison_mettl/analysis/sir/",
                          mOverlap = .05, s_gtf, plotAlignments = FALSE, transcripts, translations) {
  system(paste0("mkdir ",  output_directory))
  pdir <- system.file(package="SpliceImpactR")
  ##get bg for all classes

  # Annotate control / test groups
  data_df$utc <- "control"
  data_df$utc[data_df$phenotype_names == unique(data_df$phenotype_names)[2]] <- "test"
  messageControl <- paste0(unique(data_df$phenotype_names)[1], ": control group")
  messageTest <- paste0(unique(data_df$phenotype_names)[2], ": test group")
  message(messageControl)
  message(messageTest)


  control_group <- data_df$sample_names[data_df$utc == "control"]
  test_group <- data_df$sample_names[data_df$utc == "test"]

  bg_param <- suppressWarnings(is.na(bg_pre))
  if (bg_param) {
    bg_input <- gsub("[^/]*$", "", c(control_group, test_group))
    bg <- getBackground(input=bg_input,
                        mOverlap,
                        cores,
                        exon_type = as_types[1],
                        pdir,
                        output_location = output_directory, s_gtf, translations)
  } else {
    bg <- bg_pre
  }
  lapply(as_types, function(x) {
    messageOut <- paste0(x, " analysis...")
    message(messageOut)
    system(paste0("mkdir ",  paste0(output_directory, x, "/")))

    fAS <- getfxnlASoutcome(output_location = paste0(output_directory, x, "/"),
                            test_group = test_group,
                            control_group = control_group,
                            data_df = data_df,
                            exon_type = x,
                            cutoff = cutoff,
                            outlier_handle = outlier_handle,
                            cores = cores,
                            tti_location = tti_location,
                            full_pipe = TRUE,
                            mOverlap = mOverlap,
                            bg = bg,
                            gtf=s_gtf,
                            plotAlignments,
                            transcripts, translations)
  })

}
