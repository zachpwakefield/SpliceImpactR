test_that("getPaired returns expected output structure", {

  pdir <- system.file("extdata", package="SpliceImpactR")
  dataDirectory <- paste0(pdir, "/")
  skip_if_not(dir.exists(dataDirectory), "Data directory not found.")


  test_group <- paste0(dataDirectory, "rawData/", c("test1", "test2", "test3"))
  control_group <- paste0(dataDirectory, "rawData/", c("control1", "control2", "control3"))

  transcripts_transDF <- file.path(dataDirectory, "transcripts_limited_transDF.csv")
  transcripts_c_trans <- file.path(dataDirectory, "transcripts_limited_c_trans.csv")
  gtf_limited <- file.path(dataDirectory, "gtf_limited.csv")
  transcript_gtf_limited <- file.path(dataDirectory, "transcript_gtf_limited.csv")
  translations_limited <- file.path(dataDirectory, "translations_limited.csv")
  biomart_data <- file.path(dataDirectory, "biomart_data_sample.csv")

  skip_if_not(file.exists(transcripts_transDF), "Missing transcripts_limited_transDF.csv")
  skip_if_not(file.exists(transcripts_c_trans), "Missing transcripts_limited_c_trans.csv")
  skip_if_not(file.exists(gtf_limited), "Missing gtf_limited.csv")
  skip_if_not(file.exists(transcript_gtf_limited), "Missing transcript_gtf_limited.csv")
  skip_if_not(file.exists(translations_limited), "Missing translations_limited.csv")
  skip_if_not(file.exists(biomart_data), "Missing biomart_data_sample.csv")

  transcripts_sample <- list(
    transDF = readr::read_csv(transcripts_transDF),
    c_trans = readr::read_lines(transcripts_c_trans)
  )
  gtf_sample <- list(
    gtf = readr::read_csv(gtf_limited),
    transcript_gtf = readr::read_csv(transcript_gtf_limited)
  )
  translations_sample <- readr::read_lines(translations_limited)
  biomart_data_sample <- readr::read_csv(biomart_data)

  result <- differential_inclusion_HITindex(
    test_names = test_group,
    control_names = control_group,
    et = "AFE",
    outlier_threshold = "Inf",
    minReads = 10,
    min_prop_samples = 0,
    chosen_method = "qbGLM"
  )

  fg <- getForeground(
    input = result,
    test_names = test_group,
    control_names = control_group,
    thresh = 0.1,
    fdr = 0.05,
    mOverlap = 0.1,
    exon_type = "AFE",
    output_location = NULL,
    cores = 1,
    gtf = gtf_sample,
    max_zero_prop = 1,
    min_prop_samples = 0,
    translations = translations_sample
  )

  pfg <- getPaired(
    foreground = fg$proBed,
    et = "AFE",
    nucleotides = transcripts_sample,
    newGTF = gtf_sample,
    cores = 1,
    output_location = NULL,
    saveAlignments = FALSE,
    exon_data = biomart_data_sample
  )

  expect_true(is.list(pfg), "Expected a list from getPaired()")
  expect_equal(length(pfg), 3)

  expect_true(is.data.frame(pfg[[1]]), "First element must be a data frame")
  expect_true(is.data.frame(pfg[[2]]), "Second element must be a data frame")

  expect_true(inherits(pfg[[3]], "ggplot"), "Third element should be a ggplot object")
})
