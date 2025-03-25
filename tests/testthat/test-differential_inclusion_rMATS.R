test_that("differential_inclusion_rMATS works", {
  test_group <- c(
    "/projectnb2/evolution/zwakefield/SpliceImpactR/tests/testdata/rawData/test1",
    "/projectnb2/evolution/zwakefield/SpliceImpactR/tests/testdata/rawData/test2",
    "/projectnb2/evolution/zwakefield/SpliceImpactR/tests/testdata/rawData/test3"
  )
  control_group <- c(
    "/projectnb2/evolution/zwakefield/SpliceImpactR/tests/testdata/rawData/control1",
    "/projectnb2/evolution/zwakefield/SpliceImpactR/tests/testdata/rawData/control2",
    "/projectnb2/evolution/zwakefield/SpliceImpactR/tests/testdata/rawData/control3"
  )

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
