getfxnlASoutcome <- funciton(output_location,
                             test_group,control_group, exon_type, cutoff = .25, cores = 4, translations, transcripts, gtf, tti_location = "") {
  system(paste0("mkdir ",  output_location))
  pdir <- system.file(package="SpliceImpactR")
  if (length(test_group) > 3 & length(control_group) > 3 ) {
    boolUse <- T
  } else {boolUse <- F}

  diHIT <- SpliceImpactR::differential_inclusion_HITindex(test_names = test_group, control_names = control_group, cores = cores,stat_model_bool = boolUse, outlier_bool = boolUse, outlier_threshold = "4/n", min_proportion_samples_per_phenotype = .15)
  if (exon_type %in% c("AFE", "ALE", "HFE", "HLE")) {
    if (exon_type %in% c("AFE", "HFE")) {
      diAS <- diHIT[diHIT$type == "AFE",]
    } else {
      diAS <- diHIT[diHIT$type == "ALE",]
    }
  } else {
    diAS <- SpliceImpactR::differential_inclusion_rMATS(test_names = test_group, control_names = control_group, et = exon_type, cores = cores, outlier_threshold = "4/n", min_proportion_samples_per_phenotype = .15)
  }


  c_trans <- SpliceImpactR::get_c_trans(translations)

  c_nucs <- SpliceImpactR::get_c_nucs(transcripts)

  newGTF <- SpliceImpactR::setup_gtf(gtf)
  gtf <- newGTF$gtf

  fg <- SpliceImpactR::getForeground(input=diAS,
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

  pfg <- SpliceImpactR::getPaired(foreground = fg$proBed, et = exon_type, nucleotides = c_nucs)

  if (tti_location != "") {
    iDDI <- init_ddi(pdir = pdir, output_location = output_location, ppidm_class = "Gold_Standard", removeDups = T)
    tti_location <- output_location
  }


  bg_input <- gsub("[^/]*$", "", c(control_group, test_group))
  bg <- SpliceImpactR::getBackground(input=bg_input,
                                     mOverlap = .5,
                                     cores = cores,
                                     nC = length(control_group),
                                     nE = length(test_group),
                                     exon_type = exon_type,
                                     pdir = pdir,
                                     output_location = output_location)

  #####
  pfam <- SpliceImpactR::getPfam(foreground = fg, background = bg, pdir = pdir, cores = cpres, output_location = output_location)


  gD <- SpliceImpactR::getData(output_location = output_location, fdr_use = .05, min_sample_success = 3, engine = "Pfam")


  tti <- SpliceImpactR::getTTI(paired_foreground = pfg$paired_proBed, background = bg$proBed,
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
              diAS = diAS,
              newGTF = newGTF,
              gtf = gtf,
              c_trans = c_trans,
              c_nucs = c_nucs))
}
