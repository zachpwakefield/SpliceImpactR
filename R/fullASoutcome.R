fullASoutcome <- function(as_types = c("AFE", "ALE", "HFE", "HLE", "SE", "MXE", "RI", "A5SS", "A3SS"),
                          output_directory, data_directory,
                          test_group, control_group,
                          cutoff = .2, cores = 6,
                          tti_location = "/projectnb/evolution/zwakefield/allison_mettl/analysis/sir/") {
  system(paste0("mkdir ",  output_directory))
  pdir <- system.file(package="SpliceImpactR")
  fullOutcome <- lapply(as_types, function(x) {
    print(paste0(x, " analysis..."))
    system(paste0("mkdir ",  paste0(output_directory, x, "/")))

    fAS <- getfxnlASoutcome2(output_location = paste0(output_directory, x, "/"),
                             test_group = test_group,control_group = control_group,
                             exon_type = x, cutoff = .2, cores = 6,
                             tti_location = tti_location)
    fAS
  })

  return(fullOutcome)
}
