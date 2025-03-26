test_that("getHitCompare returns two ggplot objects and one data frame", {

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

  compareHIT <- getHitCompare(
    data_df,
    output_location = NULL,
    threshold = 0.1
  )

  # Basic checks on the returned list
  expect_true(is.list(compareHIT), "Expected a list as output.")
  expect_equal(length(compareHIT), 3)

  expect_true(is.list(compareHIT[[1]]),
              "First element should be a list of ggplot objects.")
  expect_equal(length(compareHIT[[1]]), 2)

    for (i in seq_along(compareHIT[[1]])) {
    expect_true(inherits(compareHIT[[1]][[i]], "ggplot"),
                paste("Sub-element", i, "of compareHIT[[1]] must be a ggplot object."))
  }
  expect_true(inherits(compareHIT[[2]], "pheatmap"),
              "Second element should be a ggplot object.")

  expect_true(is.data.frame(compareHIT[[3]]),
              "Third element should be a data frame.")
})
