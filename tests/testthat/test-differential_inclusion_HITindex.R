test_that("differential_inclusion_HITindex works", {
  # Paths to test and control sample directories (already provided)
  test_group2 <- c(
    "/projectnb2/evolution/zwakefield/SpliceImpactR/tests/testdata/rawData/test1",
    "/projectnb2/evolution/zwakefield/SpliceImpactR/tests/testdata/rawData/test2",
    "/projectnb2/evolution/zwakefield/SpliceImpactR/tests/testdata/rawData/test3"
  )
  control_group2 <- c(
    "/projectnb2/evolution/zwakefield/SpliceImpactR/tests/testdata/rawData/control1",
    "/projectnb2/evolution/zwakefield/SpliceImpactR/tests/testdata/rawData/control2",
    "/projectnb2/evolution/zwakefield/SpliceImpactR/tests/testdata/rawData/control3"
  )

  result <- differential_inclusion_HITindex(
    test_names = test_group2,
    control_names = control_group2,
    et = "AFE",
    outlier_threshold = "Inf",
    minReads = 10,
    min_prop_samples = 0,    # as indicated in your example call
    chosen_method = "qbGLM"
  )
  # Basic checks
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) >= 0)
  expect_equal(ncol(result), 12)

  expected_cols <- c("gene", "exon", "p.val", "delta.psi",
                     "test_average_psi", "control_average_psi",
                     "count_test", "count_control", "zero_count",
                     "p.adj", "add_inf", "type")
  expect_identical(colnames(result), expected_cols)
})
