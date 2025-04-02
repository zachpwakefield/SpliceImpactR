test_that("getOverviewComparison works with example data", {

  pdir <- system.file("extdata", package="SpliceImpactR")
  dataDirectory <- paste0(pdir, "/rawData/")

  skip_if_not(dir.exists(dataDirectory), "Data directory not found.")

  test_group <- paste0(dataDirectory, c("test1", "test2", "test3"))
  control_group <- paste0(dataDirectory, c("control1", "control2", "control3"))

  data_df <- data.frame(
    sample_names = c(control_group, test_group),
    phenotype_names = c(
      rep("control", length(control_group)),
      rep("test", length(test_group))
    ),
    stringsAsFactors = FALSE
  )

  data_df$utc <- "control"
  data_df$utc[data_df$phenotype_names == unique(data_df$phenotype_names)[2]] <- "test"

  overview <- getOverviewComparison(
    data_df,
    "AFE",
    output_location = NULL,
    minReads = 10
  )

  expect_true(is.list(overview) || inherits(overview, "list"),
              info = "Output should be a list of ggplot objects.")
  expect_equal(length(overview), 4,
               info = "Output should contain 4 separate ggplot plots.")

  all_ggplots <- vapply(overview, function(x) inherits(x, "ggplot"), logical(1))
  expect_true(all(all_ggplots),
              info = "All returned objects should be ggplot plots.")
})
