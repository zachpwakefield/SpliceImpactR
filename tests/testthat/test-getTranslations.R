test_that("getTranslations works with sample data", {
  pdir <- system.file("extdata", package = "SpliceImpactR")
  dataDirectory <- paste0(pdir, "/")

  translations <- getTranslations(translations_location = dataDirectory)

  expect_type(translations, "character")
  expect_gt(length(translations), 0)
})
