getfxnlASoutcome <- function(output_location,
                             test_group,control_group,
                             exon_type, cutoff = .25,
                             cores = 4,
                             tti_location = "") {
  system(paste0("mkdir ",  output_location))
  pdir <- system.file(package="SpliceImpactR")
  if (length(test_group) > 3 & length(control_group) > 3 ) {
    boolUse <- T
  } else {boolUse <- F}

  diHIT <- differential_inclusion_HITindex(test_names = test_group, control_names = control_group,
                                                          cores = cores,stat_model_bool = boolUse, outlier_bool = boolUse,
                                                          outlier_threshold = "4/n", min_proportion_samples_per_phenotype = .15)
  if (exon_type %in% c("AFE", "HFE")) {
      diAS <- diHIT[diHIT$type == "AFE",]
  } else if (exon_type %in% c("AFE", "HFE")) {
      diAS <- diHIT[diHIT$type == "ALE",]
    } else {
    diAS <- differential_inclusion_rMATS(test_names = test_group, control_names = control_group,
                                                        et = exon_type, cores = cores, outlier_threshold = "4/n",
                                                        min_proportion_samples_per_phenotype = .15)
  }


  fg <- getForeground(input=diAS,
                     test_names = test_group,
                     control_names = control_group,
                     thresh = cutoff,
                     fdr=.05,
                     mOverlap=.5,
                     cores=cores,
                     nC=length(control_group),
                     nE=length(test_group),
                     exon_type=exon_type,
                     pdir=pdir,
                     output_location=output_location)

  pfg <- getPaired(foreground = fg$proBed, et = exon_type, nucleotides = c_nucs)

  if (tti_location == "") {
    iDDI <- init_ddi(pdir = pdir, output_location = output_location, ppidm_class = "Gold_Standard", removeDups = T)
    tti_location <- output_location
  }


  bg_input <- gsub("[^/]*$", "", c(control_group, test_group))
  bg <- getBackground(input=bg_input,
                     mOverlap = .5,
                     cores = cores,
                     nC = length(control_group),
                     nE = length(test_group),
                     exon_type = exon_type,
                     pdir = pdir,
                     output_location = output_location)

  #####
  pfam <- getPfam(foreground = fg, background = bg, pdir = pdir, cores = cores, output_location = output_location)


  gD <- getData(fg = fg, bg = bg, fg_out = pfam$fg_out, bg_out = pfam$bg_out, output_location = output_location,
                fdr_use = .05, min_sample_success = 3, engine = "Pfam")


  tti <- getTTI(paired_foreground = pfg$paired_proBed, background = bg$proBed,
                               pdir = pdir,
                               steps=1,
                               max_vertices_for_viz = 5000,
                               fdr = .05,
                               plot_bool = T,
                               ppidm_class = "Gold_Standard",
                               write_igraphs = T,
                               ddi = "Gold",
                               output_location = output_location,
                               tti_location = tti_location)

  return(list(fg = fg,
              pfg = pfg,
              bg = bg,
              tti = tti,
              pfam = pfam,
              gD = gD,
              diHIT = diHIT,
              diAS = diAS))
}
