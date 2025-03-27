test_that("getfxnlASoutcome returns a list of length 10 with specific sub-structures", {
  pdir <- system.file("extdata", package = "SpliceImpactR")
  dataDirectory <- paste0(pdir, "/")
  test_group <- paste0(dataDirectory, "rawData/", c("test1", "test2", "test3"))
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
    transDF  = readr::read_csv(paste0(dataDirectory, "transcripts_limited_transDF.csv")),
    c_trans = readr::read_lines(paste0(dataDirectory, "transcripts_limited_c_trans.csv"))
  )

  gtf_sample <- list(
    gtf            = readr::read_csv(paste0(dataDirectory, "gtf_limited.csv")),
    transcript_gtf = readr::read_csv(paste0(dataDirectory, "transcript_gtf_limited.csv")),
    tgp_biomart    = readr::read_csv(paste0(dataDirectory, "tgp_biomart_limited"))
  )

  translations_sample <- readr::read_lines(paste0(dataDirectory, "translations_limited.csv"))

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

  oneASrun <- getfxnlASoutcome(
    output_location = NULL,
    test_group      = test_group,
    control_group   = control_group,
    data_df         = data_df,
    exon_type       = "AFE",
    cutoff          = 0.1,
    outlier_handle  = "Inf",
    cores           = 1,
    tti_location    = NULL,
    full_pipe       = FALSE,
    bg              = NA,
    mOverlap        = 0.05,
    gtf      = gtf_sample,
    plotAlignments  = FALSE,
    transcripts     = transcripts_sample,
    translations    = translations_sample,
    biomart_data    = biomart_data_sample,
    max_zero_prop   = 1,
    min_prop_samples = 0,
    chosen_method   = "qbGLM",
    initTTI         = initDDI
  )
  expect_type(oneASrun, "list")
  expect_length(oneASrun, 10)

  expected_top_names <- c(
    "diAS", "fg", "pfg", "bg", "pfam",
    "gD", "tti", "proxPlot", "length_comparison", "initial_comparison"
  )
  expect_true(all(expected_top_names %in% names(oneASrun)))

  expect_s3_class(oneASrun$diAS, "data.frame")

  fg <- oneASrun$fg
  expect_type(fg, "list")
  expect_length(fg, 4)
  expect_s3_class(fg$proBed,   "data.frame")
  expect_type(fg$proFast,     "character")
  expect_type(fg$protCode,    "character")
  expect_s3_class(fg$matched, "data.frame")

  pfg <- oneASrun$pfg
  expect_type(pfg, "list")
  expect_length(pfg, 3)
  expect_s3_class(pfg$exon_pairs,    "data.frame")
  expect_s3_class(pfg$paired_proBed, "data.frame")
  expect_s3_class(pfg$gdf,          "ggplot")

  bg <- oneASrun$bg
  expect_type(bg, "list")
  expect_length(bg, 4)
  expect_s3_class(bg$matched, "data.frame")
  expect_s3_class(bg$bed,     "data.frame")
  expect_s3_class(bg$proBed,  "data.frame")
  expect_type(bg$proFast,     "character")

  pfam <- oneASrun$pfam
  expect_type(pfam, "list")
  expect_length(pfam, 2)
  expect_s3_class(pfam$fg_out, "data.frame")
  expect_s3_class(pfam$bg_out, "data.frame")

  gD <- oneASrun$gD
  expect_type(gD, "list")
  expect_length(gD, 2)
  expect_s3_class(gD$data, "data.frame")
  expect_type(gD$enrichmentPlots, "logical")

  tti <- oneASrun$tti
  expect_type(tti, "list")
  expect_length(tti, 2)
  expect_true(is.atomic(tti$differences) || is.list(tti$differences))
  expect_s3_class(tti$results, "data.frame")

  proxPlot <- oneASrun$proxPlot
  expect_type(proxPlot, "list")
  expect_length(proxPlot, 2)
  expect_s3_class(proxPlot$proxPlot,      "ggplot")
  expect_s3_class(proxPlot$proximalShift, "table")

  expect_s3_class(oneASrun$length_comparison, "ggplot")

  init_comp <- oneASrun$initial_comparison
  expect_type(init_comp, "list")
  expect_length(init_comp, 4)
  expect_s3_class(init_comp$p1,        "ggplot")
  expect_s3_class(init_comp$p2,        "ggplot")
  expect_s3_class(init_comp$p3,        "ggplot")
  expect_s3_class(init_comp$comb_plot, "ggplot")
})
