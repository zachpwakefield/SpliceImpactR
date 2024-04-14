fullASoutcome <- function(as_types = c("AFE", "ALE", "HFE", "HLE", "SE", "MXE", "RI", "A5SS", "A3SS"),
                          output_directory, data_directory,
                          test_group, control_group, outlier_handle,
                          cutoff = .2, cores = 6,
                          tti_location = "/projectnb/evolution/zwakefield/allison_mettl/analysis/sir/") {
  system(paste0("mkdir ",  output_directory))
  pdir <- system.file(package="SpliceImpactR")
  ##get bg for all classes
  bg_input <- gsub("[^/]*$", "", c(control_group, test_group))
  bg <- getBackground(input=bg_input,
                      mOverlap = .5,
                      cores = cores,
                      nC = length(control_group),
                      nE = length(test_group),
                      exon_type = as_types[1],
                      pdir = pdir,
                      output_location = output_location)
  lapply(as_types, function(x) {
    print(paste0(x, " analysis..."))
    system(paste0("mkdir ",  paste0(output_directory, x, "/")))

    fAS <- getfxnlASoutcome(output_location = paste0(output_directory, x, "/"),
                             test_group = test_group,control_group = control_group,
                             exon_type = x, cutoff = .2, outlier_handle = outlier_handle, cores = 6,
                             tti_location = tti_location, full_pipe = T)
  })

}
