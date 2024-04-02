## This function organizes sample data by creating directories for each sample and copying specific files into them.
## It is designed to work within a specific directory structure that includes 'rmats' and 'hit' folders.


## if in structure:
# sample
# # rmats
# # hit

organizeSamples <- function(samples, output_location, cores = 1) {

  # Create a main directory to store all the data
  system(paste0("mkdir ", paste0(output_location, 'data/')))

  # Process each sample to organize its files
  sample_names <- unlist(parallel::mclapply(samples, mc.cores = cores, function(y) {

    # Extract the sample name from the path
    sample <- strsplit(y, '/')[[1]][length(strsplit(y, '/')[[1]])]

    # Create a directory for the sample within the main data directory
    system(paste0("mkdir ", paste0(output_location, 'data/', sample)))

    # Define the directory containing splicing event files (SE.MATS.JC.txt)
    se_dir <- paste0(y, 'rmats/', "SE.MATS.JC.txt")

    # Copy the splicing event file to the sample directory, renaming it with a '.SEPSI' suffix
    system(
      paste0("cp ", se_dir, " ", output_location, 'data/', sample, '/',
             sample, '.SEPSI')
    )

    # Define the directory containing splicing event files (SE.MATS.JC.txt)
    a5ss_dir <- paste0(y, 'rmats/', "A5SS.MATS.JC.txt")

    # Copy the splicing event file to the sample directory, renaming it with a '.SEPSI' suffix
    system(
      paste0("cp ", a5ss_dir, " ", output_location, 'data/', sample, '/',
             sample, '.A5SSPSI')
    )

    # Define the directory containing splicing event files (SE.MATS.JC.txt)
    a3ss_dir <- paste0(y, 'rmats/', "A3SS.MATS.JC.txt")

    # Copy the splicing event file to the sample directory, renaming it with a '.SEPSI' suffix
    system(
      paste0("cp ", a3ss_dir, " ", output_location, 'data/', sample, '/',
             sample, '.A3SSPSI')
    )

    # Define the directory containing splicing event files (SE.MATS.JC.txt)
    mxe_dir <- paste0(y, 'rmats/', "MXE.MATS.JC.txt")

    # Copy the splicing event file to the sample directory, renaming it with a '.SEPSI' suffix
    system(
      paste0("cp ", mxe_dir, " ", output_location, 'data/', sample, '/',
             sample, '.MXEPSI')
    )

    # Define the directory containing splicing event files (SE.MATS.JC.txt)
    ri_dir <- paste0(y, 'rmats/', "RI.MATS.JC.txt")

    # Copy the splicing event file to the sample directory, renaming it with a '.SEPSI' suffix
    system(
      paste0("cp ", ri_dir, " ", output_location, 'data/', sample, '/',
             sample, '.RIPSI')
    )

    # Define the directory containing hit files and copy specific types of hit files to the sample directory
    hit_dir <- paste0(y, 'hit/')

    # Copy exon hit files to the sample directory
    system(
      paste0("cp ", hit_dir, list.files(hit_dir)[grep("exon", list.files(hit_dir))], " ", output_location, 'data/', sample, '/')
    )

    # Copy ALE hit files to the sample directory
    system(
      paste0("cp ", hit_dir, list.files(hit_dir)[grep("ALE", list.files(hit_dir))], " ", output_location, 'data/', sample, '/')
    )

    # Copy AFE hit files to the sample directory
    system(
      paste0("cp ", hit_dir, list.files(hit_dir)[grep("AFE", list.files(hit_dir))], " ", output_location, 'data/', sample, '/')
    )


  }))
}
