# proteinExtract_pipe
# inputs:
### files_dir : either a list of hit exon output directories for background domain generation or a differential inclusion output file for foreground domains
### differential inclusion final output : dataframe output from differential inclusion function
### background : T for generating background domains or F for generating foreground (differentially included) domains
#                defauls to T
### updown : whether looking at exons differentially included more (up) or less (down)
#            defaults to up
### thresh : delta psi threshold absolute value
#            defaults to 0.4
### fdr : false discovery rate to count as significant in differential inclusion of exons
#         defaults to 0.05
### mOverlap : proportion of overlap of exons to count as intersection
#              defaults to 0.5
### saveOutput : whether to write output or not
#                defaults to F
### inCores : number of cores for use
#             defaults to 8
### nC : number of control samples
#        defaults to 0
### nE : number of experimental samples
#        defaults to 0
### exon_type : AFE, ALE, or SE
#               defaults to AFE
### location : package location, defaults to system.file(package="domainEnrichment")
### output_location : directory to write output

# outputs:
### matched : transcripts matched to input exons
### bed : bed file format for exons of transcripts matched to input exons
### proBed : file of transcripts matched to input exons with protein code
### proFast : fasta file for proteins matched to input exons
### paired_bed : bed file for exons of transcripts matched to input paired exons (background = F)
### paired_proBed : bed file of transcripts matched to input exons (background = F)
### paired_proFast : fasta file for proteins matched to input exons (background = F)
### df.l : differential exon inclusion output with log2 fold change (background = F)
### gdf : alignment plot (background = F)
### deExons : volcano plot (background = F)

proteinExtract_pipe <- function(files_dir, diff_incl_df, background = T, updown = c('up', 'down')[1], thresh = .3, fdr = .05, mOverlap = .5,
                                saveOutput = F, inCores = 8, nC = 0, nE = 0, exon_type = "AFE",
                                location = system.file(package="domainEnrichment"), output_location) {

  if (background == T) {

    ## If using background set, extract all first exons and create combined data.frame with gene, location
    files <- paste(files_dir, list.files(files_dir)[grep('[.]exon', list.files(files_dir))], sep = "")
    cat(files)
    if (exon_type == "AFE") {
      lim <- c("first")
    } else if (exon_type == "ALE") {
      lim <- c("last")
    } else if (exon_type == "SE") {
      lim <- c("internal")
    }
    first_exons <- unique(unlist(lapply(files, function(x) {
      in_file <- read.delim(x)
      in_file <- in_file[in_file$ID == lim,]
      paste(in_file$gene, ';', in_file$exon, ';',  in_file$strand, sep = "")})))
    redExon <- data.frame(geneR = unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 1)), split = '[.]'), "[[", 1)),
                          chr = unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 1)), split = ":"), "[[", 1)),
                          start = unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 1)), split = ":"), "[[", 2)),
                          stop = unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 2))
    )

    ## Remove duplicate rows
    redExon <- redExon[!duplicated(redExon),]

  } else {

    ## If using foreground set, read in diExon file and extract differentially included exons using lfc()
    df <- diff_incl_df
    df.l <- lfc(df, numCont = nC, numExp = nE, exon_type = exon_type, cores = inCores)

    ## Filter df.l for paired analysis
    bdf.l <- df.l[abs(df.l$delta_PSI) >= thresh & df.l$p_value <= fdr,]

    ## Make volcano plot with make_lfcPlot()
    lfcPlot <- make_lfcPlot(df.l)


    ## Filter based on fdr and thresh inputs for all analysis
    if (updown == "up") {
      cdf.l <- df.l[df.l$delta_PSI >= thresh & df.l$p_value <= fdr,]
    } else {
      cdf.l <- df.l[df.l$delta_PSI <= -(thresh) & df.l$p_value <= fdr,]
    }

    ## Make data.frame with gene, location of each exon on up/down filtered redExon
    redExon_filt <- data.frame(geneR = unlist(lapply(strsplit(cdf.l$gene, split = "[.]"), "[[", 1)),
                          chr = sapply(strsplit(cdf.l$exon, split = ":"), "[[", 1),
                          start = sapply(strsplit(sapply(strsplit(cdf.l$exon, split = ":"), "[[", 2), split = "[-]"), "[[", 1),
                          stop = sapply(strsplit(sapply(strsplit(cdf.l$exon, split = ":"), "[[", 2), split = "[-]"), "[[", 2)
    )
    colnames(redExon_filt) <- c("geneR", "chr", "start", "stop")
    redExon_filt$start <- as.numeric(redExon_filt$start)
    redExon_filt$stop <- as.numeric(redExon_filt$stop)

    ## Make data.frame with gene, location of each exon on total foreground
    redExon <- data.frame(geneR = unlist(lapply(strsplit(bdf.l$gene, split = "[.]"), "[[", 1)),
                               chr = sapply(strsplit(bdf.l$exon, split = ":"), "[[", 1),
                               start = sapply(strsplit(sapply(strsplit(bdf.l$exon, split = ":"), "[[", 2), split = "[-]"), "[[", 1),
                               stop = sapply(strsplit(sapply(strsplit(bdf.l$exon, split = ":"), "[[", 2), split = "[-]"), "[[", 2)
    )
  }

  ## Standardize redExon column names and column types
  colnames(redExon) <- c("geneR", "chr", "start", "stop")
  redExon$start <- as.numeric(redExon$start)
  redExon$stop <- as.numeric(redExon$stop)
  print("exon loaded...")

  ## Use getTranscript() to extract total and (if background = F) paired transcripts matched to input exons
  matched <- getTranscript(gtf = gtf, redExon = redExon, ex_type = exon_type, minOverlap = mOverlap, swaps = !(background), cores = inCores)
  print("exons matched, bed-ifying...")

  ## Use bedify() to extract the total or matched (if background = F) bed file
  ## num = 1 is paired in foreground, and total in background

  bed <- bedify(matched, num = 1, saveBED=F, outname = outname, cores = inCores)
  print("done bed-ifying...")

  ## extract unqiue transcript names as trans and all trancript names as possT
  trans <- unlist(lapply(strsplit(unique(bed$name), "#"), "[[", 1))
  possT <- unlist(lapply(strsplit(bed$name, "#"), "[[", 1))


  ## Find annotated proteins for transcripts if possible
  protCode <- unlist(parallel::mclapply(trans, mc.cores = 8, function(x) {
    rc <- c_trans[which(c_trans == x)+1]
    if (length(rc) > 0) {
      rc[1]
    } else {"none"}
  }))

  if (background) {print("making background output...")} else {print("making paired output...")}

  ## Make dataframe proBed for output of matched transcripts with protein code
  proBed <- data.frame(id = unique(bed$name), strand = unlist(lapply(unique(bed$name), function(x) unique(bed$strand[bed$name == x][1]))), prot = protCode) %>%
    tidyr::separate(id, c("transcript", "id"), "#") %>%
    tidyr::separate("id", c("gene", "chr"), ";") %>%
    tidyr::separate('chr', c('chr', 'coords'), ':') %>%
    tidyr::separate('coords', c('start', 'stop'), '-')

  ## Make fasta file with id & strand in first line and protein code
  proFast <- c()
  for (i in 1:length(proBed[,1])) {
    proFast <- c(proFast, paste(">", proBed$transcript[i], "#", proBed$gene[i], ";", proBed$chr[i], ":", proBed$start[i], "-", proBed$stop[i], ";", proBed$strand[i], sep = ""),
                 proBed$prot[i])
  }


  if (background == T) {
    if (saveOutput == T) {
      write_csv(proBed, paste0(output_location, "bgoutBed.csv"))
      write_lines(proFast, paste0(output_location, "bgoutFast.fa"))
      write_csv(matched$out_matched,  paste0(output_location, "bgmatched.csv"))
      write_csv(bed,  paste0(output_location, "bgbed.csv"))
    }
    return(list(matched = matched,
                bed = bed,
                proBed = proBed,
                proFast = proFast))
  } else {

    ## In paired bed file, mark each pair as either onePC, nonPC, PC and as either Match, Different, Frameshift, or Partial Match
    ## Possible output of aligned protein sequence pdf
    protAlign <- list()
    protC <- c()
    pMatch <- c()
    alignType <- c()
    cate <- c()
    for (i in seq(from=1,to=(length(protCode)-1), by=2)) {
      if (protCode[i] == "none" | protCode[i+1] == "none") {
        if (protCode[i] == "none" & protCode[i+1] != "none") {
          protC <- c(protC, "nonPC", "PC")
          protAlign[[i]] <- "none"
          protAlign[[i]] <- "onePC"
          alignType <- c(alignType, "onePC")
        } else if (protCode[i] != "none" & protCode[i+1] == "none") {
          protC <- c(protC, "PC", "nonPC")
          protAlign[[i]] <- "onePC"
          alignType <- c(alignType, "onePC")
        } else {
          protC <- c(protC, "nonPC", "nonPC")
          protAlign[[i]] <- "none"
          alignType <- c(alignType, "noPC")
        }
        pMatch <- c(pMatch, 0)
      } else if (protCode[i] == protCode[i+1]) {
        protC <- c(protC, "Same", "Same")
        protAlign[[i]] <- msa::msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])))
        pMatch <- c(pMatch, 1.04)
        alignType <- c(alignType, "Match")
      } else {
        protC <- c(protC, "Different", "Different")
        protAlign[[i]] <- msa::msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE)

        minPc <- min(nchar(protCode[i]), nchar(protCode[i+1]))
        pMatch <- c(pMatch, table(unlist(lapply(strsplit(msa::msaConsensusSequence(protAlign[[i]]), split = ""), function(x) x == "?")))[1]/min(nchar(protCode[i]), nchar(protCode[i+1])))
        if (nchar(paste(strsplit(msa::msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]][nchar(strsplit(msa::msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]]) > (.1*minPc)], collapse = "")) > .2*minPc)  {
          alignType <- c(alignType, "PartialMatch")
          # msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE
          # , file = paste(out_dir, "prettyAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_pm_prettyAlignment.pdf", sep = ""), output = "pdf")
        } else {
          alignType <- c(alignType, "FrameShift")
          # msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE
          # , file = paste(out_dir, "prettyAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_fs_prettyAlignment.pdf", sep = ""), output = "pdf")
        }
      }
    }
    print(table(protC))
    proBed$match <- protC
    proBed$prop <- rep(pMatch, each = 2)

    gdf_df <- data.frame(dens = as.numeric(pMatch), type = alignType)
    gdf_df2 <- gdf_df[gdf_df$type != "noPC",]
    # Alignment plot showing distribution of different type of exon swapping
    (gdf <- ggplot2::ggplot(gdf_df, ggplot2::aes(x = dens, fill = type)) +
        ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(count)/sum(ggplot2::after_stat(count))), colour = 1,
                                bins = 20) + ggplot2::geom_density(ggplot2::aes(y=.0005*ggplot2::after_stat(count)), color = 'black', fill = "coral2", bw = .1, alpha = .3) +
        ggplot2::scale_fill_manual(values=c('noPC' = "azure4", 'Match' = "#E69F00", 'onePC' = "#56B4E9", 'FrameShift' = "pink", 'PartialMatch' = "deeppink4")) +
        ggplot2::theme_classic() + ggplot2::xlab("Alignment Score") + ggplot2::ylab("Fraction"))

    # Alignment plot showing distribution of different type of exon swapping
    (gdf2 <- ggplot2::ggplot(gdf_df2, ggplot2::aes(x = dens, fill = type)) +
        ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(count)/sum(ggplot2::after_stat(count))), colour = 1,
                                bins = 20) + ggplot2::geom_density(ggplot2::aes(y=.0005*ggplot2::after_stat(count)), color = 'black', fill = "coral2", bw = .1, alpha = .3) +
        ggplot2::scale_fill_manual(values=c('noPC' = "azure4", 'Match' = "#E69F00", 'onePC' = "#56B4E9", 'FrameShift' = "pink", 'PartialMatch' = "deeppink4")) +
        ggplot2::theme_classic() + ggplot2::xlab("Alignment Score") + ggplot2::ylab("Fraction"))


    proBed$matchType <- rep(alignType, each = 2)

    print("making all output...")
    ## Reproduce previous output but on entire differentially includeded foreground set, not just matched transcripts for domain enrichment analysis
    matched_all <- getTranscript(gtf = gtf, redExon = redExon_filt, ex_type = exon_type, minOverlap = mOverlap, swaps = !(background), cores = inCores)
    bed_all <- bedify(matched_all, num = 3, saveBED=F, outname = outname, cores = inCores)
    trans_all <- unlist(lapply(strsplit(unique(bed_all$name), "#"), "[[", 1))
    possT_all <- unlist(lapply(strsplit(bed_all$name, "#"), "[[", 1))

    ## Find annotated proteins for transcripts if possible
    protCode_all <- unlist(mclapply(trans_all, mc.cores = 8, function(x) {
      rc <- c_trans[which(c_trans == x)+1]
      if (length(rc) > 0) {
        rc[1]
      } else {"none"}
    }))

    ## Make dataframe proBed for output of matched transcripts withprotein code
    proBed_all <- data.frame(id = unique(bed_all$name), strand = unlist(lapply(unique(bed_all$name), function(x) unique(bed_all$strand[bed_all$name == x][1])[1])), prot = protCode_all) %>%
      tidyr::separate(id, c("transcript", "id"), "#") %>%
      tidyr::separate("id", c("gene", "chr"), ";") %>%
      tidyr::separate('chr', c('chr', 'coords'), ':') %>%
      tidyr::separate('coords', c('start', 'stop'), '-')

    ## Make fasta file with id & strand in first line and protein code
    proFast_all <- c()
    for (i in 1:length(proBed_all[,1])) {
      proFast_all <- c(proFast_all, paste(">", proBed_all$transcript[i], "#", proBed_all$gene[i], ";", proBed_all$chr[i], ":", proBed_all$start[i], "-", proBed_all$stop[i], ";", proBed_all$strand[i], sep = ""),
                       proBed_all$prot[i])
    }
    if (saveOutput == T) {
      write_csv(proBed_all, paste0(output_location, "fgoutBed.csv"))
      write_lines(proFast_all, paste0(output_location, "fgoutFast.fa"))
      write_csv(matched$tot_matched,  paste0(output_location, "fgmatched.csv"))
      write_csv(bed_all,  paste0(output_location, "fgbed.csv"))
      write_csv(proBed, paste0(output_location, "paired_fgoutBed.csv"))
      write_lines(proFast, paste0(output_location, "paired_fgoutFast.fa"))
      write_csv(matched_all$out_matched,  paste0(output_location, "paired_fgmatched.csv"))
      write_csv(bed,  paste0(output_location, "paired_fgbed.csv"))
      write_csv(df.l,  paste0(output_location, "fglfc.csv"))

      pdf(file = paste0(output_location, "alignPlot.pdf"))
      print(gdf)
      dev.off()

      pdf(file = paste0(output_location, "alignPlot_pc.pdf"))
      print(gdf2)
      dev.off()

      pdf(file = paste0(output_location, "volcano.pdf"))
      print(lfcPlot)
      dev.off()
    }
    return(list(matched = matched,
                bed = bed_all,
                proBed = proBed_all,
                proFast = proFast_all,
                paired_bed = bed,
                paired_proBed = proBed,
                paired_proFast = proFast,
                df.l = df.l,
                bdf.l = bdf.l,
                cdf.l = cdf.l,
                gdf = gdf,
                gdf2 = gdf2,
                deExons = lfcPlot))
  }

}
