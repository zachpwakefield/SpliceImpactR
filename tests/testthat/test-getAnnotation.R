test_that("getAnnotation produces the expected structure", {
  pdir <- system.file("extdata", package = "SpliceImpactR")
  dataDirectory <- paste0(pdir, "/")

  biomart_data_sample <- list(
    setup_gtf_exon_data = readr::read_csv(paste0(dataDirectory, "biomart_setup_gtf_exon_data.csv")),
    ptg_init            = readr::read_csv(paste0(dataDirectory, "biomart_ptg_init.csv")),
    transcript_data     = readr::read_csv(paste0(dataDirectory, "biomart_transcript_data.csv"))
  )

  gtf_list <- getAnnotation(biomart_data = biomart_data_sample)

  expect_type(gtf_list, "list")

  expected_names <- c(
    "gtf",
    "transcript_gtf",
    "hybrid_last_extract",
    "hybrid_first_extract",
    "tgp_biomart",
    "hybrid_first_extract_transcripts",
    "hybrid_last_extract_transcripts"
  )

  expect_true(all(expected_names %in% names(gtf_list)))

  expect_s3_class(gtf_list$gtf,                       "data.frame")
  expect_s3_class(gtf_list$transcript_gtf,            "data.frame")
  expect_s3_class(gtf_list$hybrid_last_extract,       "data.frame")
  expect_s3_class(gtf_list$hybrid_first_extract,      "data.frame")
  expect_s3_class(gtf_list$tgp_biomart,               "data.frame")

  expect_type(gtf_list$hybrid_first_extract_transcripts, "character")
  expect_type(gtf_list$hybrid_last_extract_transcripts,  "character")
})
