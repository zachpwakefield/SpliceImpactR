## This function retrieves and processes data for domain enrichment analysis using various protein domain databases.
## It compares domain occurrences in foreground (differentially expressed) and background datasets to identify enriched domains.

getData <- function(fg, bg, pfg, pfam, cores = cores,
                    output_location, fdr_use, min_sample_success,
                    engine = c("FunFam","Gene3D","CDD","PANTHER","SMART","ProSiteProfiles","Pfam","SUPERFAMILY","MobiDBLite","Coils","PRINTS","ProSitePatterns","PIRSF","NCBIfam","Hamap")[7],
                    repeatingDomains = F, topViz = 15) {

  ## extract bed, fasta, and interproscan files from output_location
  tf <- list.files(output_location)
  ipscan <- list(pfam$bg_out, pfam$fg_out) #fg_out, bg_out
  outFast<- list(bg$proFast, fg$proFast) #fg$proFast, fg$proFast
  outBed <- list(bg$proBed, pfg$paired_proBed) #bg$proBed, bg$proBed

  system(paste0("mkdir ", output_location, "DomainEnrichment/"))
  ## Process interproscan results for background and foreground datasets
  interproscan_results <- lapply(1:2, function(o) {
    # Read domain scan results and sequence information
    ip <- ipscan[[o]]
    s <- outBed[[o]]
    s$id <- paste(paste(paste(paste(s$transcript, s$gene, sep = "#"), s$chr, sep = ";"), paste(s$start, s$stop, sep = "-"), sep = ":"), s$strand, sep = ";")
    fa <- outFast[[o]]
    if (o == 2) {
      s <- s[s$alignType %in% c("TooLong", "PartialMatch", "FrameShift"),]
      fa <- fa[sort(c(which(fa %in% paste0(">", s$id)), (which(fa %in% paste0(">", s$id))+1)))]
      nFa <- fa[grep(">", fa)]
      gN <- gsub(">", "", nFa)
      gN <- unlist(lapply(s$id, function(faIter) {gN[gN %in% faIter]}))
    } else {
      nFa <- fa[grep(">", fa)]
      gN <- gsub(">", "", nFa)
    }



    ## filter by whichever domain identifying engine(s) selected
    ips <- ip
    # ips <- ip %>% dplyr::filter(X4 %in% engine)

    # Process and aggregate information for each protein domain
    protInf <- mclapply(1:length(gN), mc.cores = cores, function(j) {
      paste(ips$X6[ips$X1 == gN[j]], collapse = ';')
    })
    protInf.o <- unlist(protInf)
    protInf.o[protInf.o == ""] <- "none"

    # Construct final output data frames for background and foreground
    if (o == 2) {
      finOut <- data.frame(exon = gN,
                           protein = s$prot[match(gN, s$id)],
                           protInfor = protInf.o,
                           delta.psi = s$delta.psi[match(gN, s$id)],
                           p.adj = s$p.adj[match(gN, s$id)])
    } else {
      finOut <- data.frame(exon = gN,
                           protein = s$prot[match(gN, s$id)],
                           protInfor = protInf.o)

    }
    finOut
  })

  # Perform domain enrichment analysis for upregulated and downregulated genes
  data <- lapply(c("up", "down"), function(x) {
    # Separate upregulated and downregulated genes based on delta.psi
    if (x == "up") {
      fg_ip <- interproscan_results[[2]][interproscan_results[[2]]$delta.psi > 0,]
      fg_sc <- interproscan_results[[2]][interproscan_results[[2]]$delta.psi < 0,]
    } else {
      fg_ip <- interproscan_results[[2]][interproscan_results[[2]]$delta.psi < 0,]
      fg_sc <- interproscan_results[[2]][interproscan_results[[2]]$delta.psi > 0,]
    }


    # Extract key values for phyper() hypergeometric
    if (repeatingDomains == F) {
      bg_dom <- unlist(lapply(interproscan_results[[1]]$protInfor, function(dom) {
        unique(unlist(strsplit(dom, split = ";")))
      }))
      fg_dom <- unlist(lapply(1:nrow(fg_ip), function(ind) {
        setdiff(unique(unlist(strsplit(fg_ip$protInfor[ind], split = ";"))), unique(unlist(strsplit(fg_sc$protInfor[ind], split = ";"))))
      }))
    } else if (repeatingDomains) {
      bg_dom <- unlist(strsplit(interproscan_results[[1]]$protInfor, split = ';'))
      fg_dom <- unlist(lapply(1:nrow(fg_ip), function(ind) {
        setdiff(unlist(strsplit(fg_ip$protInfor[ind], split = ";")), unlist(strsplit(fg_sc$protInfor[ind], split = ";")))
      }))
    }

    fg_dom_li <- lapply(strsplit(fg_ip$protInfor, split = ';'), unique)
    searcher <- unique(fg_dom)

    successes <- lapply(searcher, function(x) c(sum(fg_dom == x), sum(bg_dom == x)))
    pop_size <- length(bg_dom)
    sample_size <- length(fg_dom)

    # Compute hypergeometric test for domain enrichment
    pvals <- suppressWarnings(signif(stats::phyper(sapply(successes, "[[", 1)-1,
                                                   sapply(successes, "[[", 2),
                                                   pop_size-sapply(successes, "[[", 2),
                                                   sample_size, lower.tail=FALSE), 5))


    # Adjust p-values for multiple hypothesis testing (FDR correction)
    pvals.adj <- signif(stats::p.adjust(pvals, method="fdr"), 5)

    ## Make summary data frame
    data <- data.frame(domain = searcher,
                       pval = pvals,
                       fdr = pvals.adj,
                       sample_size = rep(length(fg_dom), length(pvals)),
                       sample_successes = sapply(successes, "[[", 1),
                       sample_prop = sapply(successes, "[[", 1) / length(fg_dom),

                       pop_size = rep(length(bg_dom), length(pvals)),
                       pop_successes = sapply(successes, "[[", 2),
                       pop_prop = sapply(successes, "[[", 2) / length(bg_dom),
                       protein_contain = unlist(lapply(searcher, function(y) {paste(unique(fg_ip$exon[which(unlist(lapply(fg_dom_li, function(uu) {y %in% uu})))]), collapse = "%")}))

    )

    data$rel_prop <- (data$sample_prop + .01) / (data$pop_prop + .01)
    data <- data[data$domain != "none" & data$domain != "" & data$domain != "-",] %>% dplyr::arrange(fdr)

    # Write domain enrichment results to CSV files
    write.csv(data, paste0(output_location, "DomainEnrichment/", ifelse(x == "up", "(+)", "(-)"), "domainEnrichment.csv"))
    data


  })

  # Generate enrichment plots for upregulated and downregulated genes
  eP <-lapply(1:2, function(x) {
    # Prepare data for plotting
    sp <- data[[x]][data[[x]]$fdr <= fdr_use & data[[x]]$sample_successes >= min_sample_success,c(1, 6, 11)]
    if (nrow(sp) == 0) {
      print("Lower min_sample_success and/or increase fdr, none meet criteria")
      return(NA)
    }
    sp$category <- "foreground"
    colnames(sp)[2] <- c("proportion")
    pp <- data[[x]][data[[x]]$fdr <= fdr_use & data[[x]]$sample_successes >= min_sample_success,c(1, 9, 11)]
    pp$category <- "background"
    colnames(pp)[2] <- c("proportion")
    df <- rbind(sp, pp)
    df <- df[df$rel_prop >= ifelse(nrow(df) > topViz, sort(df$rel_prop, decreasing = T)[topViz], max(df$rel_prop)),]
    # Create and save enrichment plots using ggplot2
    enrichmentPlot <- ggplot2::ggplot(data=df, ggplot2::aes(x=domain, y=proportion, fill=category)) +
      ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge())+
      ggplot2::coord_flip() +
      ggplot2::theme_bw()+
      ggplot2::theme(axis.ticks.y=ggplot2::element_blank()  #remove y axis ticks
      )+
      ggplot2::scale_fill_manual(values=c('brown','chartreuse4')) + ggplot2::ggtitle(paste0("Domain Enrichment of ", ifelse(x == "1", "(+)", "(-)"),
                                                                                            " PSI Exons")) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black"))

    pdf(paste0(output_location, "DomainEnrichment/", ifelse(x == "1", "(+)", "(-)"), 'enrichmentPlot.pdf'))
    print(enrichmentPlot)
    dev.off()
    enrichmentPlot

  })

  # Return the enrichment plot object for display
  # print(eP)

  # Return a list containing the domain enrichment data and plots
  return(list(data=data,
              enrichmentPlots = eP))
}
