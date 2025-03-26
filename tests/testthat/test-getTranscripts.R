test_that("getTranscripts works with sample data", {
  pdir <- system.file("extdata", package = "SpliceImpactR")
  dataDirectory <- paste0(pdir, "/")

  transcripts <- getTranscripts(transcripts_location = dataDirectory)


  expect_type(transcripts, "list")
  expect_length(transcripts, 2)

  expect_s3_class(transcripts$transDF, "data.frame")

  expect_type(transcripts$c_trans, "character")
})
