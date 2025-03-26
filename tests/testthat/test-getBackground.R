test_that("getBackground returns the expected list structure", {

  pdir <- system.file("extdata", package="SpliceImpactR")
  dataDirectory <- paste0(pdir, "/")
  skip_if_not(dir.exists(dataDirectory), "Data directory not found.")

  test_group <- paste0(dataDirectory, "rawData/", c("test1","test2", "test3"))
  control_group <- paste0(dataDirectory, "rawData/", c("control1","control2", "control3"))

  gtf_file <- file.path(dataDirectory, "gtf_limited.csv")
  transcript_file <- file.path(dataDirectory, "transcript_gtf_limited.csv")
  translations_file <- file.path(dataDirectory, "translations_limited.csv")

  skip_if_not(file.exists(gtf_file), "Missing gtf_limited.csv")
  skip_if_not(file.exists(transcript_file), "Missing transcript_gtf_limited.csv")
  skip_if_not(file.exists(translations_file), "Missing translations_limited.csv")

  # Load the example data
  gtf_sample <- list(
    gtf = readr::read_csv(gtf_file),
    transcript_gtf = readr::read_csv(transcript_file)
  )
  translations_sample <- readr::read_lines(translations_file)

  # Run the function
  bg <- getBackground(
    input = c(test_group, control_group),
    mOverlap = 0.1,
    cores = 1,
    exon_type = "AFE",
    output_location = NULL,
    gtf = gtf_sample,
    translations = translations_sample
  )

  # Basic checks
  expect_true(is.list(bg), "getBackground should return a list.")
  expect_equal(length(bg), 4)

  # First three should be data frames
  for (i in 1:3) {
    expect_true(is.data.frame(bg[[i]]), paste("Element", i, "should be a data frame."))
  }

  # Fourth should be a vector
  expect_true(is.vector(bg[[4]]), "Fourth element should be a vector.")
})
