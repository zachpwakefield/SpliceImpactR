organizeSamples <- function(samples, cores = 1) {
  system(paste0("mkdir ", paste0(output_location, 'data/')))
  sample_names <- unlist(parallel::mclapply(samples, mc.cores = cores, function(y) {
    sample <- strsplit(y, '/')[[1]][length(strsplit(y, '/')[[1]])]
    system(paste0("mkdir ", paste0(output_location, 'data/', sample)))
    se_dir <- paste0(y, 'rmats/', "SE.MATS.JC.txt")
    system(
      paste0("cp ", se_dir, " ", output_location, 'data/', sample, '/',
             sample, '.SEPSI')
    )
    hit_dir <- paste0(y, 'hit/')
    system(
      paste0("cp ", hit_dir, list.files(hit_dir)[grep("exon", list.files(hit_dir))], " ", output_location, 'data/', sample, '/')
    )
    system(
      paste0("cp ", hit_dir, list.files(hit_dir)[grep("ALE", list.files(hit_dir))], " ", output_location, 'data/', sample, '/')
    )
    system(
      paste0("cp ", hit_dir, list.files(hit_dir)[grep("AFE", list.files(hit_dir))], " ", output_location, 'data/', sample, '/')
    )


  }))
}
