test_that("getPfam works with sample data", {
  pdir <- system.file("extdata", package="SpliceImpactR")
  dataDirectory <- paste0(pdir, "/")
  test_group <- paste0(dataDirectory, "rawData/", c("test1","test2", "test3"))
  control_group <- paste0(dataDirectory, "rawData/", c("control1", "control2", "control3"))

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

  transcripts_sample <- list(
    transDF = readr::read_csv(paste0(dataDirectory, "transcripts_limited_transDF.csv")),
    c_trans = readr::read_lines(paste0(dataDirectory, "transcripts_limited_c_trans.csv"))
  )

  gtf_sample <- list(
    gtf           = readr::read_csv(paste0(dataDirectory, "gtf_limited.csv")),
    transcript_gtf = readr::read_csv(paste0(dataDirectory, "transcript_gtf_limited.csv"))
  )

  translations_sample <- readr::read_lines(paste0(dataDirectory, "translations_limited.csv"))

  biomart_data <- list(
    ip               = readr::read_csv(paste0(dataDirectory, "biomart_ip.csv")),
    code_regions     = readr::read_csv(paste0(dataDirectory, "biomart_code_regions.csv")),
    pfam_exon_level  = readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon_level.csv")),
    fsd_exon_data    = readr::read_csv(paste0(dataDirectory, "biomart_data_sample.csv"))
  )

  result <- differential_inclusion_HITindex(
    test_names        = test_group,
    control_names     = control_group,
    et                = "AFE",
    outlier_threshold = "Inf",
    minReads          = 10,
    min_prop_samples  = 0,
    chosen_method     = "qbGLM"
  )

  fg <- getForeground(
    input            = result,
    test_names       = test_group,
    control_names    = control_group,
    thresh           = 0.1,
    fdr              = 0.05,
    mOverlap         = 0.1,
    exon_type        = "AFE",
    output_location  = NULL,
    cores            = 1,
    gtf              = gtf_sample,
    max_zero_prop    = 1,
    min_prop_samples = 0,
    translations     = translations_sample
  )

  bg <- getBackground(
    input           = c(test_group, control_group),
    mOverlap        = 0.1,
    cores           = 1,
    exon_type       = "AFE",
    output_location = NULL,
    gtf_sample,
    translations_sample
  )

  pfg <- getPaired(
    foreground      = fg$proBed,
    et              = "AFE",
    nucleotides     = transcripts_sample,
    newGTF          = gtf_sample,
    cores           = 1,
    output_location = NULL,
    saveAlignments  = FALSE,
    exon_data       = biomart_data$fsd_exon_data
  )

  pfamData <- getPfam(
    background      = bg,
    foreground      = fg,
    pdir            = pdir,
    output_location = NULL,
    cores           = 1,
    biomart_data    = biomart_data
  )

  expect_type(pfamData, "list")
  expect_length(pfamData, 2)

  expect_s3_class(pfamData[[1]], "data.frame")
  expect_s3_class(pfamData[[2]], "data.frame")
})
