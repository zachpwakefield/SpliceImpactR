test_that("fullASoutcome returns a list of length 4 with correct sub-structures", {
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


  initTTI <- init_ddi(pdir = dataDirectory,
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
                                initTTI = initTTI,
                                mOverlap = .05,
                                s_gtf = gtf_sample,
                                plotAlignments = FALSE,
                                transcripts = transcripts_sample,
                                translations = translations_sample,
                                biomart_data = biomart_data_sample,
                                max_zero_prop = 1,
                                min_prop_samples = 0,
                                chosen_method = 'qbGLM')
  expect_type(twoASfullRun, "list")
  expect_length(twoASfullRun, 4)

  expected_top_names <- c("AFE", "SE", "HIT", "Background")
  expect_true(all(expected_top_names %in% names(twoASfullRun)))

  AFE <- twoASfullRun$AFE

  expect_type(AFE, "list")
  expect_length(AFE, 10)
  expected_AFE_names <- c(
    "diAS", "fg", "pfg", "bg", "pfam",
    "gD", "tti", "proxPlot", "length_comparison", "initial_comparison"
  )
  expect_true(all(expected_AFE_names %in% names(AFE)))


  expect_s3_class(AFE$diAS, "data.frame")
  expect_type(AFE$fg, "list")
  expect_type(AFE$pfg, "list")

  SE <- twoASfullRun$SE
  expect_type(SE, "list")
  expect_length(SE, 2)
  expect_true(all(c("diAS", "fg") %in% names(SE)))
  expect_s3_class(SE$diAS, "data.frame")
  fg_SE <- SE$fg
  expect_type(fg_SE, "list")
  expect_length(fg_SE, 4)
  expect_s3_class(fg_SE$proBed, "data.frame")
  expect_type(fg_SE$proFast,    "character")
  expect_type(fg_SE$protCode,   "character")
  expect_s3_class(fg_SE$matched,"data.frame")

  HIT <- twoASfullRun$HIT
  expect_type(HIT, "list")
  expect_length(HIT, 3)
  expect_s3_class(HIT$diHIT, "data.frame")

  diPlot <- HIT$diPlot
  expect_type(diPlot, "list")
  expect_length(diPlot, 2)
  expect_s3_class(diPlot[[1]], "ggplot")
  expect_s3_class(diPlot[[2]], "ggplot")
  expect_s3_class(HIT$meanHeatmap, "pheatmap")

  Background <- twoASfullRun$Background
  expect_type(Background, "list")
  expect_length(Background, 4)
  expect_s3_class(Background$matched, "data.frame")
  expect_s3_class(Background$bed,     "data.frame")
  expect_s3_class(Background$proBed,  "data.frame")
  expect_type(Background$proFast,     "character")
})
