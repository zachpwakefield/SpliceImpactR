test_that("getForeground returns the expected list structure", {

  pdir <- system.file("extdata", package="SpliceImpactR")
  dataDirectory <- paste0(pdir, "/")
  skip_if_not(dir.exists(dataDirectory), "Data directory not found.")

  test_group <- paste0(dataDirectory, "rawData/", c("test1","test2","test3"))
  control_group <- paste0(dataDirectory, "rawData/", c("control1","control2","control3"))

  gtf_file <- paste0(dataDirectory, "gtf_limited.csv")
  transcript_file <- paste0(dataDirectory, "transcript_gtf_limited.csv")
  translations_file <- paste0(dataDirectory, "translations_limited.csv")

  skip_if_not(file.exists(gtf_file), "Missing gtf_limited.csv")
  skip_if_not(file.exists(transcript_file), "Missing transcript_gtf_limited.csv")
  skip_if_not(file.exists(translations_file), "Missing translations_limited.csv")

  result <- differential_inclusion_HITindex(
    test_names = test_group,
    control_names = control_group,
    et = "AFE",
    outlier_threshold = "Inf",
    minReads = 10,
    min_prop_samples = 0,
    chosen_method = "qbGLM"
  )

  gtf_sample <- list(
    gtf = readr::read_csv(gtf_file),
    transcript_gtf = readr::read_csv(transcript_file)
  )
  translations_sample <- readr::read_lines(translations_file)

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

  expect_true(is.list(fg), "Output should be a list.")
  expect_equal(length(fg), 4)

  expect_true(is.data.frame(fg[[1]]), "First element must be a data frame.")
  expect_true(is.data.frame(fg[[4]]), "Fourth element must be a data frame.")

  expect_true(is.character(fg[[2]]), "Second element must be a character vector.")
  expect_true(is.character(fg[[3]]), "Third element must be a character vector.")
})
