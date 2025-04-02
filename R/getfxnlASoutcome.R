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
#' @param initTTI takes in output from the respective function for passing to getTTI
#' @param outlier_handle method for handling outliers
#' @param cores number of cores to assign
#' @param chosen_method statistical method to choose (qbGLM, nbGLM, zingGLM, wilcox)
#' @param biomart_data data straight from setup_biomart
#' @return nothing or all output from pipeline
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
#' data_df$utc <- "control"
#' data_df$utc[data_df$phenotype_names == unique(data_df$phenotype_names)[2]] <- "test"
#'
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
#'
#' ip <- readr::read_csv(paste0(dataDirectory, "biomart_ip.csv"))
#' code_regions <- readr::read_csv(paste0(dataDirectory, "biomart_code_regions.csv"))
#' pfam_exon_level <- readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon_level.csv"))
#' fsd_exon_data <- readr::read_csv(paste0(dataDirectory, "biomart_data_sample.csv"))
#' pfam_data = readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon.csv"))
#' biomart_data_sample <- list(ip = ip,
#'                      code_regions = code_regions,
#'                      fsd_exon_data = fsd_exon_data,
#'                      pfam_exon_level = pfam_exon_level,
#'                      pfam_data = pfam_data)
#'
#' initDDI <- init_ddi(pdir = dataDirectory,
#'                     output_location = NULL,
#'                     ppidm_class = c("Gold_Standard", "Gold", "Silver", "Bronze")[1],
#'                     removeDups = TRUE,
#'                     cores = 1,
#'                     pfam_data = biomart_data_sample$pfam_data)
#' oneASrun <- getfxnlASoutcome(output_location = NULL,
#'                              test_group,
#'                              control_group,
#'                              data_df,
#'                              exon_type = "AFE",
#'                              cutoff = .1,
#'                              outlier_handle = "Inf",
#'                              cores = 1,
#'                              tti_location = NULL,
#'                              full_pipe = FALSE,
#'                              bg = NA,
#'                              mOverlap = .05,
#'                              gtf_sample,
#'                              plotAlignments = FALSE,
#'                              transcripts = transcripts_sample,
#'                              translations = translations_sample,
#'                              biomart_data = biomart_data_sample,
#'                              max_zero_prop = 1,
#'                              min_prop_samples = 0,
#'                              chosen_method = 'qbGLM',
#'                              initTTI = initDDI)
#' @export
getfxnlASoutcome <- function(output_location = NULL,
                             test_group,control_group,data_df,
                             exon_type, cutoff = .1, outlier_handle = "Inf",
                             cores = 1,
                             tti_location = NULL, full_pipe = TRUE,
                             bg = NA, mOverlap = .05, gtf,
                             plotAlignments = FALSE,
                             transcripts, translations,
                             biomart_data,
                             max_zero_prop = .5,
                             min_prop_samples = .5,
                             chosen_method,
                             initTTI = NULL) {
  if (!is.null(output_location)) {
    system(paste0("mkdir ",  output_location))
  }
  pdir <- system.file(package="SpliceImpactR")

  if (exon_type %in% c("AFE", "HFE")) {
      diHIT <- differential_inclusion_HITindex(test_names = test_group, control_names = control_group, et = "AFE",
                                               outlier_threshold = outlier_handle, minReads = 10,
                                               min_prop_samples, chosen_method)

      diAS <- diHIT[diHIT$type == "AFE",]
  } else if (exon_type %in% c("ALE", "HLE")) {
      diHIT <- differential_inclusion_HITindex(test_names = test_group, control_names = control_group, et = "ALE",
                                               outlier_threshold = outlier_handle, minReads = 10,
                                               min_prop_samples, chosen_method)

      diAS <- diHIT[diHIT$type == "ALE",]
    } else {
      if (chosen_method == 'zinbGLM' | chosen_method == 'nbGLM') {
        chosen_method_rmats <- 'qbGLM'
      } else {
        chosen_method_rmats <- chosen_method
      }
    diAS <- differential_inclusion_rMATS(test_names = test_group, control_names = control_group,
                                         et = exon_type, outlier_threshold = outlier_handle,
                                         minReads = 10,
                                         min_prop_samples, chosen_method_rmats)
  }


  fg <- getForeground(input=diAS,
                     test_names = test_group,
                     control_names = control_group,
                     thresh = cutoff,
                     fdr=.05,
                     mOverlap=mOverlap,
                     cores=cores,
                     exon_type=exon_type,
                     output_location=output_location, gtf=gtf,
                     max_zero_prop, min_prop_samples, translations)

  if (!(exon_type %in% c("HFE", "HLE"))) {
    initial_comparison <- getOverviewComparison(data_df, exon_type, output_location)
  }

  pfg <- getPaired(foreground = fg$proBed, et = exon_type, nucleotides = transcripts,
                   output_location = output_location, newGTF = gtf, saveAlignments = plotAlignments,
                   exon_data = biomart_data$fsd_exon_data)

  if (length(pfg) > 1) {
    if (exon_type %in% c("AFE", "ALE")) {
      proxPlot <- getProximalShift(exon_type, pfg$exon_pairs, pfg$paired_proBed, output_location)
    }
    length_comparison <- getLengthComparison(data_df, pfg$paired_proBed, output_location)



  bg_param <- sum(is.na(bg)) == length(bg)
  if (bg_param) {
    bg <- getBackground(input=c(control_group, test_group),
                        mOverlap = mOverlap,
                        cores = cores,
                        exon_type = exon_type,
                        output_location = output_location,
                        gtf=gtf,
                        translations)
  }


  #####
  pfam <- getPfam(foreground = fg,
                  background = bg,
                  pdir = pdir,
                  cores = cores,
                  output_location = output_location,
                  biomart_data)


  gD <- getDomainData(fg = fg,
                      bg = bg,
                      pfg=pfg,
                      cores = cores,
                      pfam = pfam,
                      output_location = output_location,
                      fdr_use = .05,
                      min_sample_success = 2,
                      engine = "Pfam",
                      topViz = 15)

  if (is.null(tti_location) & is.null(initTTI)) {
    initTTI <- init_ddi(pdir = dataDirectory,
                        output_location = output_location,
                        ppidm_class = "Gold_Standard",
                        cores = 1,
                        removeDups = TRUE,
                        biomart_data$pfam_data)
  }

  if (nrow(pfg$paired_proBed) > 1) {
    tti <- getTTI(paired_foreground = pfg$paired_proBed,
                  background = bg$proBed,
                  steps=1,
                  max_vertices_for_viz = 2000,
                  fdr = .05,
                  ppidm_class = "Gold_Standard",
                  write_igraphs = FALSE,
                  output_location = output_location,
                  tti_location = tti_location,
                  tgp = gtf$tgp_biomart,
                  init_edgelist = initTTI$edgelist)
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
                tti = tti,
                proxPlot = proxPlot,
                length_comparison = length_comparison,
                initial_comparison = initial_comparison))
  }
  } else {
    return(list(diAS = diAS,
                fg = fg))
                 }

}
