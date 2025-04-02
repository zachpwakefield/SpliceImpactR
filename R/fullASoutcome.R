#' Perform full analysis on various alternative RNA processing events
#'
#' @param as_types a vector of the type of alternative splicing events that are being investigated
#' @param output_directory the path (ending with '/') that the output is desired
#' @param data_directory the path to the data, in the format specified by organizeSamples
#' @param outlier_handle value to threshold the outlier detection
#' @param cutoff diInclusion cutoff to use to identify significance
#' @param plotAlignments whether to plot alignments from getPaired/matchAlignType
#' @param bg_pre if bg was made earlier, param to give premade bg
#' @param s_gtf output from getAnnotation
#' @param tti_location location of previously made transcript-transcript interactions network
#' @param mOverlap minimum overlap to call an exon as matched to annotation
#' @param translations from getTranslations
#' @param transcripts from getTranscripts
#' @param max_zero_prop max prop of samples that can be 0
#' @param min_prop_samples min prop of samples from each phenotype required to show a specific event
#' @param initTTI the edge list from initTTI
#' @param data_df data frame containing one column of the control_group and test_group (sample_names) and one column of phenotype names (phenotype_names)
#' @param cores number of cores to allocate
#' @param biomart_data from setupBiomart
#' @param chosen_method choice of staistical model (qbGLM, nbGLM, zinbGLM, wilcox)
#' @return nothing in R, output to the output_directory
#'
#' @examples
#' pdir <- system.file("extdata", package="SpliceImpactR")
#' dataDirectory <- paste0(pdir, "/")
#' test_group <- paste0(dataDirectory, "rawData/", c("test1","test2", "test3"))
#' control_group <- paste0(dataDirectory, "rawData/", c("control1", "control2", "control3"))
#' data_df <- data.frame(
#'     sample_names = c(control_group, test_group),
#'     phenotype_names = c(
#'       rep("control", length(control_group)),
#'       rep("test", length(test_group))
#'      ),
#'    stringsAsFactors = FALSE
#'   )
#'
#' transDF <- readr::read_csv(paste0(dataDirectory, "transcripts_limited_transDF.csv"))
#' c_trans <- readr::read_lines(paste0(dataDirectory, "transcripts_limited_c_trans.csv"))
#'
#' transcripts_sample <- list(transDF = transDF,
#'                            c_trans = c_trans)
#'
#' gtf_sample <- list(gtf = readr::read_csv(paste0(dataDirectory, "gtf_limited.csv")),
#'             transcript_gtf = readr::read_csv(paste0(dataDirectory, "transcript_gtf_limited.csv")),
#'             tgp_biomart = readr::read_csv(paste0(dataDirectory, "tgp_biomart_limited"))
#'             )
#' translations_sample <- readr::read_lines(paste0(dataDirectory, "translations_limited.csv"))
#'
#' ip <- readr::read_csv(paste0(dataDirectory, "biomart_ip.csv"))
#' code_regions <- readr::read_csv(paste0(dataDirectory, "biomart_code_regions.csv"))
#' pfam_exon_level <- readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon_level.csv"))
#' fsd_exon_data <- readr::read_csv(paste0(dataDirectory, "biomart_data_sample.csv"))
#' pfam_data = readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon.csv"))
#'
#' biomart_data_sample <- list(ip = ip,
#'                      code_regions = code_regions,
#'                      fsd_exon_data = fsd_exon_data,
#'                      pfam_exon_level = pfam_exon_level,
#'                      pfam_data = pfam_data)
#'
#' twoASfullRun <- fullASoutcome(as_types = c("AFE", "SE", "HIT"),
#'                               output_directory = NULL,
#'                               data_directory = dataDirectory,
#'                               data_df,
#'                               outlier_handle = "Inf",
#'                               cutoff = .1,
#'                               cores = 1,
#'                               bg_pre = NA,
#'                               tti_location = NULL,
#'                               initTTI = NULL,
#'                               mOverlap = .05,
#'                               s_gtf = gtf_sample,
#'                               plotAlignments = FALSE,
#'                               transcripts = transcripts_sample,
#'                               translations = translations_sample,
#'                               biomart_data = biomart_data_sample,
#'                               max_zero_prop = 1,
#'                               min_prop_samples = 0,
#'                               chosen_method = 'qbGLM')
#'
#'
#'
#'
#' @export
#'

fullASoutcome <- function(as_types = c("AFE", "ALE", "HFE", "HLE", "SE", "MXE", "RI", "A5SS", "A3SS", "HIT"),
                          output_directory,
                          data_directory,
                          data_df,
                          outlier_handle = "Inf",
                          cutoff = .1,
                          cores = 1,
                          bg_pre = NA,
                          tti_location = NULL,
                          initTTI = NULL,
                          mOverlap = .05,
                          s_gtf,
                          plotAlignments = FALSE,
                          transcripts,
                          translations,
                          biomart_data,
                          max_zero_prop = .5,
                          min_prop_samples = .5,
                          chosen_method = 'nbGLM') {
  if (!is.null(output_directory)) {
    system(paste0("mkdir ",  output_directory))
  }
  pdir <- system.file(package="SpliceImpactR")
  ##get bg for all classes

  # Annotate control / test groups
  data_df$utc <- "control"
  data_df$utc[data_df$phenotype_names == unique(data_df$phenotype_names)[2]] <- "test"

  control_group <- data_df$sample_names[data_df$utc == "control"]
  test_group <- data_df$sample_names[data_df$utc == "test"]

  bg_param <- sum(is.na(bg_pre)) == length(bg_pre)
  if (bg_param) {
    bg <- getBackground(input=c(control_group, test_group),
                        mOverlap,
                        cores,
                        exon_type = as_types[1],
                        output_location = output_directory, s_gtf, translations)
  } else {
    bg <- bg_pre
  }

  if (is.null(tti_location) & is.null(initTTI)) {
    initTTI <- init_ddi(pdir = data_directory,
                        output_location = output_directory,
                        ppidm_class = "Gold_Standard",
                        cores = 1,
                        removeDups = TRUE,
                        biomart_data$pfam_data)
  }

  fullASresult <- lapply(as_types, function(x) {
    messageOut <- paste0(x, " analysis...")
    message(messageOut)
    if (!is.null(output_directory)) {
      system(paste0("mkdir ",  paste0(output_directory, x, "/")))
    }
    if (x == 'HIT') {
      if (!is.null(output_directory)) {
        system(paste0("mkdir ",  paste0(output_directory, "/HITindex")))
      }
      hitCompare <- getHitCompare(data_df,
                                  if (is.null(output_directory)) {NULL} else {paste0(output_directory, "/HITindex/")},
                                  .25)
      hitCompare
    } else {

    fAS <- getfxnlASoutcome(output_location = if (is.null(output_directory)) {NULL} else {paste0(output_directory, x, "/")},
                            test_group = test_group,
                            control_group = control_group,
                            data_df = data_df,
                            exon_type = x,
                            cutoff = cutoff,
                            outlier_handle = outlier_handle,
                            cores = cores,
                            tti_location = tti_location,
                            full_pipe = FALSE,
                            mOverlap = mOverlap,
                            bg = bg,
                            gtf=s_gtf,
                            plotAlignments,
                            transcripts, translations,
                            biomart_data = biomart_data,
                            max_zero_prop,
                            min_prop_samples,
                            chosen_method,
                            initTTI = initTTI
                            )
    fAS
    }
  })
  fullASresult[[length(fullASresult) + 1]] <- bg
  names(fullASresult) <- c(as_types, "Background")
  return(fullASresult)
}
