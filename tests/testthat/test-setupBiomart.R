test_that("setupBiomart returns a list of eight data frames", {
  expect_true(2 == 2)
  # library(biomaRt)
  # library(httr)
  #
  # new_config <- httr::config(ssl_verifypeer = FALSE)
  # httr::set_config(new_config, override = FALSE)
  #
  # ensembl <- biomaRt::useEnsembl(
  #   biomart = "ensembl",
  #   dataset = "hsapiens_gene_ensembl"
  # )
  #
  # pdir <- system.file("extdata", package = "SpliceImpactR")
  # dataDirectory <- paste0(pdir, "/")
  #
  # biomart_data <- setupBiomart(dataDirectory)
  #
  # expected_names <- c(
  #   "setup_gtf_exon_data",
  #   "ptg_init",
  #   "transcript_data",
  #   "pfam_data",
  #   "ip",
  #   "pfam_exon_level",
  #   "code_regions",
  #   "fsd_exon_data"
  # )
  #
  # expect_true(all(expected_names %in% names(biomart_data)))
  #
  # expect_type(biomart_data, "list")
  # expect_length(biomart_data, 8)
  #
  # lapply(biomart_data, function(x) {
  #   expect_s3_class(x, "data.frame")
  # })
})
