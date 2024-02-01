## if in structure:
# sample
# # rmats
# # hit

organize_samples <- function(samples, cores) {
  system(paste0("mkdir ", paste0(output_location, 'data/')))
  sample_names <- unlist(parallel::mclapply(samples, mc.cores=cores, function(y) {
    sample <- strsplit(y, '/')[[1]][length(strsplit(y, '/')[[1]])]
    system(paste0("mkdir ", paste0(output_location, 'data/', sample)))
    se_dir <- paste0(y, 'rmats/', "SE.MATS.JC.txt")
    lapply(se_dir, function(x) {
      system(
        paste0("cp ", x, " ", output_location, 'data/', unlist(lapply(strsplit(x, '/'), "[[", 8)), '/',
               unlist(lapply(strsplit(y, '/'), "[[", 8)), '.SEPSI')
      )
    })
    hit_dir <- paste0(y, 'hit/')
    lapply(hit_dir, function(x) {
      system(
        paste0("cp ", x, list.files(x)[grep("exon", list.files(x))], " ", output_location, 'data/', unlist(lapply(strsplit(x, '/'), "[[", 8)), '/')
      )
      system(
        paste0("cp ", x, list.files(x)[grep("ALE", list.files(x))], " ", output_location, 'data/', unlist(lapply(strsplit(x, '/'), "[[", 8)), '/')
      )
      system(
        paste0("cp ", x, list.files(x)[grep("AFE", list.files(x))], " ", output_location, 'data/', unlist(lapply(strsplit(x, '/'), "[[", 8)), '/')
      )
    })

  }))
}
