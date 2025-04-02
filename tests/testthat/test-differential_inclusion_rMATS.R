test_that("differential_inclusion_rMATS works", {

  pdir <- system.file("extdata", package="SpliceImpactR")
  dataDirectory <- paste0(pdir, "/rawData/")

  skip_if_not(dir.exists(dataDirectory), "Data directory not found.")

  test_group <- paste0(dataDirectory, c("test1", "test2", "test3"))
  control_group <- paste0(dataDirectory, c("control1", "control2", "control3"))

  result <- differential_inclusion_rMATS(
    test_names = test_group,
    control_names = control_group,
    et = "SE",
    outlier_threshold = "Inf",
    minReads = 10,
    min_prop_samples = 0,
    chosen_method = "qbGLM"
  )
  # Basic checks
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) >= 0)
  expect_equal(ncol(result), 13)

  expected_cols <- c("id", "p.val", "delta.psi",
                     "test_average_psi", "control_average_psi",
                     "count_test", "count_control", "zero_count",
                     "p.adj", "type", "gene", "exon", "add_inf")
  expect_identical(colnames(result), expected_cols)
})
