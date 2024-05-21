fullASoutcome <- function(as_types = c("AFE", "ALE", "HFE", "HLE", "SE", "MXE", "RI", "A5SS", "A3SS"),
                          output_directory, data_directory,
                          data_df, outlier_handle,
                          cutoff = .2, cores = 6, bg_pre = NA,
                          tti_location = "/projectnb/evolution/zwakefield/allison_mettl/analysis/sir/", mOverlap = .01) {
  system(paste0("mkdir ",  output_directory))
  pdir <- system.file(package="SpliceImpactR")
  ##get bg for all classes

  # Annotate control / test groups
  data_df$utc <- "control"
  data_df$utc[data_df$phenotype_names == unique(data_df$phenotype_names)[2]] <- "test"
  print(paste0(unique(data_df$phenotype_names)[1], ": control group"))
  print(paste0(unique(data_df$phenotype_names)[2], ": test group"))


  control_group <- data_df$sample_names[data_df$utc == "control"]
  test_group <- data_df$sample_names[data_df$utc == "test"]

  if (is.na(bg_pre[[1]])) {
    bg_input <- gsub("[^/]*$", "", c(control_group, test_group))
    bg <- getBackground(input=bg_input,
                        mOverlap = .01,
                        cores = cores,
                        nC = length(control_group),
                        nE = length(test_group),
                        exon_type = as_types[1],
                        pdir = pdir,
                        output_location = output_directory)
  } else {
    bg <- bg_pre
  }
  lapply(as_types, function(x) {
    print(paste0(x, " analysis..."))
    system(paste0("mkdir ",  paste0(output_directory, x, "/")))

    fAS <- getfxnlASoutcome(output_location = paste0(output_directory, x, "/"),
                             test_group = test_group,control_group = control_group, data_df = data_df,
                             exon_type = x, cutoff = cutoff, outlier_handle = outlier_handle, cores = cores,
                             tti_location = tti_location, full_pipe = T, mOverlap = mOverlap, bg = bg)
  })

}
