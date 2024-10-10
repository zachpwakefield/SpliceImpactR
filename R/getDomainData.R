#' idenitify enriched domains across phenotype
#'
#' @param fg foreground from getForeground
#' @param bg background from getBackground
#' @param pfam pfam from getPfam
#' @param cores the number of cores requested
#' @param fdr_use the fdr to set as a threshold
#' @param min_sample_success the number of appearances of a domain in the sample set to visualize
#' @param engine pfam only for now
#' @param repeatingDomains whether to identify repeating domains being enriched or not
#' @param topViz the max number of domains to put in each visualization
#' @param output_location location to make background directory
#' @return the domain enrichment data and the enrichment plots
#' @importFrom dplyr arrange filter first group_by
#' @importFrom stats phyper p.adjust
#' @importFrom ggplot2 scale_fill_manual theme_classic ggplot aes xlab ylab theme geom_bar geom_boxplot theme_classic coord_flip theme_bw ggtitle element_blank element_line
#' @importFrom ggpubr ggarrange
#' @importFrom data.table as.data.table
#' @export
getDomainData <- function(fg, bg, pfg, pfam, cores = 1,
                    output_location, fdr_use = .05, min_sample_success = 3,
                    engine = c("FunFam","Gene3D","CDD","PANTHER","SMART","ProSiteProfiles","Pfam","SUPERFAMILY","MobiDBLite","Coils","PRINTS","ProSitePatterns","PIRSF","NCBIfam","Hamap")[7],
                    repeatingDomains = FALSE, topViz = 15) {

  ## extract bed, fasta, and interproscan files from output_location
  tf <- list.files(output_location)
  ipscan <- list(pfam$bg_out, pfam$fg_out) #fg_out, bg_out
  outFast<- list(bg$proFast, fg$proFast) #fg$proFast, fg$proFast
  outBed <- list(bg$proBed, pfg$paired_proBed) #bg$proBed, bg$proBed

  ## remove dup'd transcripts from bg_out to make background correctly
  ipscan[[1]] <- ipscan[[1]] %>%
    dplyr::group_by(transcript) %>%
    dplyr::filter(X1 == dplyr::first(X1))

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
    ips_dt <- data.table::as.data.table(ip)

    result <- ips_dt[, .(protInf = paste(domains, collapse = ";")), by = X1]

    protInf_df <- data.frame(result[gN, on = .(X1), nomatch = NA])
    protInf_df$protInf[is.na(protInf_df$protInf)] <- "none"
    order_index <- match(protInf_df$X1, gN)

    protInf_df_ordered <- protInf_df[order(order_index), , drop = FALSE]

    protInf.o <- protInf_df_ordered$protInf

    if (o == 2) {
      finOut <- data.frame(exon = gN,
                           protein = s$prot,
                           protInfor = protInf.o,
                           delta.psi = s$delta.psi,
                           p.adj = s$p.adj)
    } else {
      finOut <- data.frame(exon = gN,
                           protein = s$prot,
                           protInfor = protInf.o)

    }
    finOut
  })

  # Perform domain enrichment analysis for upregulated and downregulated genes
  dataList <- lapply(c("up", "down"), function(x) {
    # Separate upregulated and downregulated genes based on delta.psi
    if (x == "up") {
      fg_ip <- interproscan_results[[2]][interproscan_results[[2]]$delta.psi > 0,]
      fg_sc <- interproscan_results[[2]][interproscan_results[[2]]$delta.psi < 0,]
    } else {
      fg_ip <- interproscan_results[[2]][interproscan_results[[2]]$delta.psi < 0,]
      fg_sc <- interproscan_results[[2]][interproscan_results[[2]]$delta.psi > 0,]
    }


    # Extract key values for phyper() hypergeometric
    if (repeatingDomains == FALSE) {
      bg_dom <- unlist(lapply(interproscan_results[[1]]$protInfor, function(dom) {
        unique(unlist(strsplit(dom, split = ";")))
      }))
      fg_dom_ul <- lapply(1:nrow(fg_ip), function(ind) {
        setdiff(c(unique(unlist(strsplit(fg_ip$protInfor[ind], split = ";"))), "none"), c(unique(unlist(strsplit(fg_sc$protInfor[ind], split = ";"))), "none"))
      })
      fg_dom <- unlist(fg_dom_ul)
    } else if (repeatingDomains) {
      bg_dom <- unlist(strsplit(interproscan_results[[1]]$protInfor, split = ';'))
      fg_dom_ul <- lapply(1:nrow(fg_ip), function(ind) {
        setdiff(c(unlist(strsplit(fg_ip$protInfor[ind], split = ";")), "none"), c(unlist(strsplit(fg_sc$protInfor[ind], split = ";")), "none"))
      })
      fg_dom <- unlist(fg_dom_ul)
    }
    bg_dom <- bg_dom[bg_dom != "none"]

    if (sum(lengths(fg_dom_ul)) == 0) {
      noneFound <- paste0("No domains enriched in ", x, " set")
      message(noneFound)
      return(list("none", data.frame(vals = 0, types = x), data.frame(vals = 0, types = x)))
    }


    fg_dom_li <- lapply(strsplit(fg_ip$protInfor, split = ';'), unique)
    searcher <- unique(fg_dom)

    l_dd <- lengths(fg_dom_ul)[lengths(fg_dom_ul) != 0]


    lengthsDistributionDF <- data.frame(vals = l_dd, types = rep(x, length(l_dd)))
    barPlotDF <- data.frame(vals = length(l_dd), types = x)

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

    list(data, lengthsDistributionDF, barPlotDF)


  })


  lengthDist <- do.call(rbind, lapply(dataList, "[[", 2))
  changeNum <- do.call(rbind, lapply(dataList, "[[", 3))
  lengthDist$types2 <- "(-)"
  lengthDist$types2[lengthDist$types == "up"] <- "(+)"
  changeNum$types2 <- "(-)"
  changeNum$types2[changeNum$types == "up"] <- "(+)"

  domainChangesNums <- ggplot2::ggplot(lengthDist, ggplot2::aes(x = types2, y = vals, fill = types)) + ggplot2::geom_boxplot() +
    ggplot2::scale_fill_manual(values=c("brown", "chartreuse4"), breaks = c("down", "up")) +
    ggplot2::theme_classic() + ggplot2::ylab("Number of domain changes per swap") +
    ggplot2::xlab("") +
    ggplot2::theme(legend.position = "none")

  domainChanges <- ggplot2::ggplot(changeNum, ggplot2::aes(x = types2, y = vals, fill = types)) + ggplot2::geom_bar(stat = 'identity') +
    ggplot2::scale_fill_manual(values=c("brown", "chartreuse4"), breaks = c("down", "up")) +
    ggplot2::theme_classic() + ggplot2::xlab("") + ggplot2::ylab("Count of swaps with domain changes")+
    ggplot2::theme(legend.position = "none")

  comb_plot <- ggpubr::ggarrange(domainChanges, domainChangesNums, nrow = 1, widths = c(1, 2.5))

  pdf(paste0(output_location, "DomainEnrichment/", "domainStats.pdf"), height = 8, width = 6)
  print(comb_plot)
  dev.off()

  dataList[[1]][[1]]$reg <- "Upregulated"
  dataList[[2]][[1]]$reg <- "Downregulated"

  dataFinal <- rbind(dataList[[1]][[1]], dataList[[2]][[1]])
  dataFinal$reg <- factor(dataFinal$reg, levels = unique(dataFinal$reg))
  dataFinal$domain2 <- 1:length(dataFinal$domain)
  # dataFinal <- transform(dataFinal, variable=reorder(domain, domain2) )
  enrichmentPlot <- ggplot2::ggplot(data=dataFinal, ggplot2::aes(x=reorder(domain, domain2), y=-log10(fdr), fill = reg)) +
    ggplot2::geom_bar(stat="identity")+
    ggplot2::coord_flip() +
    ggplot2::theme_bw()+
    ggplot2::theme(axis.ticks.y=ggplot2::element_blank()  #remove y axis ticks
    )+
    ggplot2::scale_fill_manual(values=c('brown','chartreuse4'), breaks = c("Downregulated", "Upregulated")) +
    ggplot2::ggtitle("Global Domain Enrichment") +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black"))

  pdf(paste0(output_location, "DomainEnrichment/GlobalEnrichmentPlot.pdf"))
  print(enrichmentPlot)
  dev.off()
  enrichmentPlot



  # Return a list containing the domain enrichment data and plots
  return(list(data=dataFinal,
              enrichmentPlots = enrichmentPlot))
}
