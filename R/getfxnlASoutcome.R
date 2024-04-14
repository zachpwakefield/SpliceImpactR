getfxnlASoutcome <- function(output_location,
                             test_group,control_group,
                             exon_type, cutoff = .25, outlier_handle = "4/n",
                             cores = 4,
                             tti_location = "", full_pipe = T, boolUse = T, bg = NA, mOverlap = .5) {
  system(paste0("mkdir ",  output_location))
  pdir <- system.file(package="SpliceImpactR")

  if (exon_type %in% c("AFE", "HFE")) {
      diHIT <- differential_inclusion_HITindex(test_names = test_group, control_names = control_group, et = "AFE",
                                               cores = cores,stat_model_bool = boolUse, outlier_bool = boolUse,
                                               outlier_threshold = outlier_handle, min_proportion_samples_per_phenotype = .15)

      diAS <- diHIT[diHIT$type == "AFE",]
  } else if (exon_type %in% c("ALE", "HLE")) {
      diHIT <- differential_inclusion_HITindex(test_names = test_group, control_names = control_group, et = "ALE",
                                               cores = cores,stat_model_bool = boolUse, outlier_bool = boolUse,
                                               outlier_threshold = outlier_handle, min_proportion_samples_per_phenotype = .15)

      diAS <- diHIT[diHIT$type == "ALE",]
    } else {
    diAS <- differential_inclusion_rMATS(test_names = test_group, control_names = control_group,
                                         stat_model_bool = boolUse, outlier_bool = boolUse,
                                         et = exon_type, cores = cores, outlier_threshold = outlier_handle,
                                         min_proportion_samples_per_phenotype = .15)
  }


  fg <- getForeground(input=diAS,
                     test_names = test_group,
                     control_names = control_group,
                     thresh = cutoff,
                     fdr=.05,
                     mOverlap=mOverlap,
                     cores=cores,
                     nC=length(control_group),
                     nE=length(test_group),
                     exon_type=exon_type,
                     pdir=pdir,
                     output_location=output_location)

  initial_comparison <- getOverviewComparison(control_group, test_group, exon_type, output_location)

  pfg <- getPaired(foreground = fg$proBed, et = exon_type, nucleotides = c_nucs, output_location = output_location, newGTF = newGTF, saveAlignments = F)

  length_comparison <- getLengthComparison(pfg$paired_proBed, output_location)

  if (tti_location == "") {
    iDDI <- init_ddi(pdir = pdir, output_location = output_location, ppidm_class = "Gold_Standard", removeDups = T)
    tti_location <- output_location
  }

  if (is.na(bg)) {
    bg_input <- gsub("[^/]*$", "", c(control_group, test_group))
    bg <- getBackground(input=bg_input,
                        mOverlap = .5,
                        cores = cores,
                        nC = length(control_group),
                        nE = length(test_group),
                        exon_type = exon_type,
                        pdir = pdir,
                        output_location = output_location)
  }


  #####
  pfam <- getPfam(foreground = fg, background = bg, pdir = pdir, cores = cores, output_location = output_location)


  gD <- getData(fg = fg, bg = bg, cores = cores, pfam = pfam, output_location = output_location,
                fdr_use = .05, min_sample_success = 5, engine = "Pfam")

  if (nrow(pfg$paired_proBed) > 1) {
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
  } else {
    tti <- NA
  }

  if (full_pipe) {
    return(NA)
  } else {
    return(list(diAS = diAS,
                fg = fg,
                pfg = pfg,
                bg = bg,
                pfam = pfam,
                gD = gD,
                tti = tti))
  }

}
