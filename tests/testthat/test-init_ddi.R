test_that("init_ddi works with sample data", {

  pdir <- system.file("extdata", package="SpliceImpactR")
  dataDirectory <- paste0(pdir, "/")

  biomart_data_sample <- list(
    ip               = readr::read_csv(paste0(dataDirectory, "biomart_ip.csv")),
    code_regions     = readr::read_csv(paste0(dataDirectory, "biomart_code_regions.csv")),
    pfam_exon_level  = readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon_level.csv")),
    fsd_exon_data    = readr::read_csv(paste0(dataDirectory, "biomart_data_sample.csv")),
    pfam_data        = readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon.csv"))
  )

  initDDI <- init_ddi(
    pdir            = dataDirectory,
    output_location = NULL,
    ppidm_class     = c("Gold_Standard", "Gold", "Silver", "Bronze")[1],
    removeDups      = TRUE,
    cores           = 1,
    pfam_data       = biomart_data_sample$pfam_data
  )

  expect_type(initDDI, "list")
  expect_length(initDDI, 2)

  expect_s3_class(initDDI[[1]], "igraph")

  expect_true(is.matrix(initDDI[[2]]))
})
