test_that("getIntegratedResults produces summaryTable and plots as expected", {
  pdir <- system.file("extdata", package="SpliceImpactR")
  dataDirectory <- paste0(pdir, "/")
  test_group <- paste0(dataDirectory, "rawData/", c("test1","test2", "test3"))
  control_group <- paste0(dataDirectory, "rawData/", c("control1", "control2", "control3"))
  data_df <- data.frame(
    sample_names = c(control_group, test_group),
    phenotype_names = c(
      rep("control", length(control_group)),
      rep("test", length(test_group))
    ),
    stringsAsFactors = FALSE
  )

  transcripts_sample <- list(transDF = readr::read_csv(paste0(dataDirectory, "transcripts_limited_transDF.csv")),
                             c_trans = readr::read_lines(paste0(dataDirectory, "transcripts_limited_c_trans.csv")))

  gtf_sample <- list(gtf = readr::read_csv(paste0(dataDirectory, "gtf_limited.csv")),
                     transcript_gtf = readr::read_csv(paste0(dataDirectory, "transcript_gtf_limited.csv")),
                     tgp_biomart = readr::read_csv(paste0(dataDirectory, "tgp_biomart_limited"))
  )
  translations_sample <- readr::read_lines(paste0(dataDirectory, "translations_limited.csv"))
  biomart_data_sample <- list(ip = readr::read_csv(paste0(dataDirectory, "biomart_ip.csv")),
                              code_regions = readr::read_csv(paste0(dataDirectory, "biomart_code_regions.csv")),
                              pfam_exon_level = readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon_level.csv")),
                              fsd_exon_data = readr::read_csv(paste0(dataDirectory, "biomart_data_sample.csv")),
                              pfam_data = readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon.csv")))


  initDDI <- init_ddi(pdir = dataDirectory,
                      output_location = NULL,
                      ppidm_class = c("Gold_Standard", "Gold", "Silver", "Bronze")[1],
                      removeDups = TRUE,
                      cores = 1,
                      pfam_data = biomart_data_sample$pfam_data)
  twoASfullRun <- fullASoutcome(as_types = c("AFE", "SE", "HIT"),
                                output_directory = NULL,
                                data_directory = dataDirectory,
                                data_df,
                                outlier_handle = "Inf",
                                cutoff = .1,
                                cores = 1,
                                bg_pre = NA,
                                tti_location = NULL,
                                initTTI = initDDI,
                                mOverlap = .05,
                                s_gtf = gtf_sample,
                                plotAlignments = FALSE,
                                transcripts = transcripts_sample,
                                translations = translations_sample,
                                biomart_data = biomart_data_sample,
                                max_zero_prop = 1,
                                min_prop_samples = 0,
                                chosen_method = 'qbGLM')

  fg_list <- list("AFE" = twoASfullRun$AFE$fg$proBed)
  pfg_list <- list("AFE" = twoASfullRun$AFE$pfg$paired_proBed)
  domain_list <- list("AFE" = twoASfullRun$AFE$gD$data)

  integrated <- getIntegratedResults(fg_list,
                                     pfg_list,
                                     domain_list)
  expect_type(integrated, "list")
  expect_length(integrated, 2)

  expected_top_names <- c("summaryTable", "plots")
  expect_true(all(expected_top_names %in% names(integrated)))

  expect_s3_class(integrated$summaryTable, "data.frame")

  plot_list <- integrated$plots
  expect_type(plot_list, "list")
  expect_length(plot_list, 9)

  expected_plot_names <- c(
    "RelativeUsePlot", "MedianAlignScorePlot", "PropRescuePlot", "DomainCountPlot",
    "RelativeDomainPlot", "pairedEventTypePlot", "geneUpSet", "domainUpSet", "alignmentScorePlot"
  )
  expect_true(all(expected_plot_names %in% names(plot_list)))


  expect_s3_class(plot_list$RelativeUsePlot,      "ggplot")
  expect_s3_class(plot_list$MedianAlignScorePlot, "ggplot")
  expect_s3_class(plot_list$PropRescuePlot,       "ggplot")
  expect_s3_class(plot_list$DomainCountPlot,      "ggplot")
  expect_s3_class(plot_list$RelativeDomainPlot,   "ggplot")
  expect_s3_class(plot_list$pairedEventTypePlot,  "ggplot")
  expect_s3_class(plot_list$alignmentScorePlot,   "ggplot")

  expect_type(plot_list$geneUpSet,   "character")
  expect_type(plot_list$domainUpSet, "character")
})
