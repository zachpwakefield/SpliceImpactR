test_that("getProximalShift works with sample data", {

  pdir <- system.file("extdata", package="SpliceImpactR")
  dataDirectory <- paste0(pdir, "/")
  test_group <- paste0(dataDirectory, "rawData/", c("test1","test2", "test3"))
  control_group <- paste0(dataDirectory, "rawData/", c("control1", "control2", "control3"))

  transcripts_sample <- list(
    transDF  = readr::read_csv(paste0(dataDirectory, "transcripts_limited_transDF.csv")),
    c_trans = readr::read_lines(paste0(dataDirectory, "transcripts_limited_c_trans.csv"))
  )

  gtf_sample <- list(
    gtf           = readr::read_csv(paste0(dataDirectory, "gtf_limited.csv")),
    transcript_gtf = readr::read_csv(paste0(dataDirectory, "transcript_gtf_limited.csv"))
  )

  translations_sample <- readr::read_lines(paste0(dataDirectory, "translations_limited.csv"))
  biomart_data_sample <- readr::read_csv(paste0(dataDirectory, "biomart_data_sample.csv"))

  result <- differential_inclusion_HITindex(
    test_names        = test_group,
    control_names     = control_group,
    et                = "AFE",
    outlier_threshold = "Inf",
    minReads          = 10,
    min_prop_samples  = 0,
    chosen_method     = "qbGLM"
  )

  fg <- getForeground(
    input           = result,
    test_names      = test_group,
    control_names   = control_group,
    thresh          = 0.1,
    fdr             = 0.05,
    mOverlap        = 0.1,
    exon_type       = "AFE",
    output_location = NULL,
    cores           = 1,
    gtf             = gtf_sample,
    max_zero_prop   = 1,
    min_prop_samples = 0,
    translations    = translations_sample
  )

  pfg <- getPaired(
    foreground      = fg$proBed,
    et              = "AFE",
    nucleotides     = transcripts_sample,
    newGTF          = gtf_sample,
    cores           = 1,
    output_location = NULL,
    saveAlignments  = FALSE,
    exon_data       = biomart_data_sample
  )

  proxShift <- getProximalShift(
    type       = "AFE",
    exon_pairs      = pfg$exon_pairs,
    ep_supp   = pfg$paired_proBed,
    output_location = NULL
  )

  expect_type(proxShift, "list")
  expect_length(proxShift, 2)
  expect_s3_class(proxShift[[1]], "ggplot")
  expect_s3_class(proxShift[[2]], "table")
})
