getData <- function(cuD, fdr_use, min_sample_success, engine = c("FunFam","Gene3D","CDD","PANTHER","SMART","ProSiteProfiles","Pfam","SUPERFAMILY","MobiDBLite","Coils","PRINTS","ProSitePatterns","PIRSF","NCBIfam","Hamap")[7]) {

  ## extract bed, fasta, and interproscan files from output_location
  tf <- list.files(cuD)
  ipscan <- c("bgoutFast.fa.tsv", "fgoutFast.fa.tsv") #tf[grep("fa.tsv", tf)]
  outFast<- c("bgoutFast.fa", "fgoutFast.fa") #tf[(!(grepl(".tsv", tf)) & !(grepl("xml", tf)) & !(grepl("gff3", tf)) & !(grepl("json", tf)) & grepl("outFast.fa", tf) & !grepl("paired", tf))]
  outBed <- c("bgoutBed.csv", "fgoutBed.csv") #tf[grepl("outBed.csv", tf) & !grepl("paired", tf)]

  ## collate interproscan results
  interproscan_results <- lapply(c("bg", "fg"), function(o) {
    ip <- readr::read_delim(paste(cuD, ipscan[grep(o, ipscan)], sep = ""), col_names = T, delim = '\t')
    s <- readr::read_csv(paste(cuD, outBed[grep(o, outBed)], sep =""))
    s$id <- paste(paste(paste(paste(s$transcript, s$gene, sep = "#"), s$chr, sep = ";"), paste(s$start, s$stop, sep = "-"), sep = ":"), s$strand, sep = ";")
    fa <- readr::read_lines(paste(cuD, outFast[grep(o, outFast)], sep = ""))
    nFa <- fa[grep(">", fa)]
    gN <- gsub(">", "", nFa)
    geneL <- unlist(lapply(strsplit(unlist(lapply(strsplit(fa[grep(">", fa)], split = ";"), "[[", 1)), split = "#"), "[[", 2))

    ## filter by whichever domain identifying engine(s) selected
    ips <- ip %>% dplyr::filter(X4 %in% engine)
    protInf <- list()
    for (j in 1:length(gN)) {
      protInf[[j]] <- ips$X6[ips$X1 == gN[j]]
      protInf[[j]] <- paste(protInf[[j]], collapse = ';')
    }
    protInf.o <- unlist(protInf)
    protInf.o[protInf.o == ""] <- "none"
    if (o == "fg") {
      finOut <- data.frame(exon = gN,
                           protein = unlist(lapply(gN, function(x) s$prot[s$id == x])),
                           protInfor = protInf.o,
                           delta.psi = unlist(lapply(gN, function(x) s$delta.psi[s$id == x])),
                           p.adj = unlist(lapply(gN, function(x) s$p.adj[s$id == x])))
    } else {
      finOut <- data.frame(exon = gN,
                           protein = unlist(lapply(gN, function(x) s$prot[s$id == x])),
                           protInfor = protInf.o)
    }
  })
  data <- lapply(c("up", "down"), function(x) {
    if (x == "up") {
      fg_ip <- interproscan_results[[2]][interproscan_results[[2]]$delta.psi > 0,]

    } else {
      fg_ip <- interproscan_results[[2]][interproscan_results[[2]]$delta.psi < 0,]
    }


    bg_dom <- unlist(strsplit(interproscan_results[[1]]$protInfor, split = ';'))
    fg_dom <- unlist(strsplit(fg_ip$protInfor, split = ';'))
    fg_dom_li <- lapply(strsplit(fg_ip$protInfor, split = ';'), unique)
    searcher <- unique(fg_dom)

    successes <- lapply(searcher, function(x) c(sum(fg_dom == x), sum(bg_dom == x)))
    pop_size <- length(bg_dom)
    sample_size <- length(searcher)

    ## compute hypergeometric test
    pvals <- suppressWarnings(signif(stats::phyper(sapply(successes, "[[", 1)-1,
                                                   sapply(successes, "[[", 2),
                                                   pop_size-sapply(successes, "[[", 2),
                                                   sample_size, lower.tail=FALSE), 5))


    ## multiple hypothesis correction, fdr
    pvals.adj <- signif(stats::p.adjust(pvals, method="fdr"), 5)

    ## output
    data <- data.frame(domain = searcher,
                       pval = pvals,
                       fdr = pvals.adj,
                       sample_size = rep(length(fg_dom), length(pvals)),
                       sample_successes = sapply(successes, "[[", 1),
                       sample_prop = sapply(successes, "[[", 1) / length(fg_dom),

                       pop_size = rep(length(bg_dom), length(pvals)),
                       pop_successes = sapply(successes, "[[", 2),
                       pop_prop = sapply(successes, "[[", 2) / length(bg_dom),
                       protein_contain = unlist(lapply(searcher, function(y) {paste(unique(interproscan_results[[2]]$exon[which(unlist(lapply(fg_dom_li, function(uu) {y %in% uu})))]), collapse = "%")}))

    )

    data$rel_prop <- (data$sample_prop + .01) / (data$pop_prop + .01)
    data <- data[data$domain != "none" & data$domain != "" & data$domain != "-",] %>% dplyr::arrange(fdr)
    write.csv(data, paste0(output_location, ifelse(x == "1", "(+)", "(-)"), "domainEnrichment.csv"))
    data


  })
  eP <-lapply(1:2, function(x) {
    sp <- data[[x]][data[[x]]$fdr <= fdr_use & data[[x]]$sample_successes >= min_sample_success,c(1, 6)]
    sp$category <- "foreground"
    colnames(sp)[2] <- c("proportion")
    pp <- data[[x]][data[[x]]$fdr <= fdr_use & data[[x]]$sample_successes >= min_sample_success,c(1, 9)]
    pp$category <- "background"
    colnames(pp)[2] <- c("proportion")
    df <- rbind(sp, pp)

    enrichmentPlot <- ggplot2::ggplot(data=df, ggplot2::aes(x=domain, y=proportion, fill=category)) +
      ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge())+
      ggplot2::coord_flip() +
      ggplot2::theme_bw()+
      ggplot2::theme(axis.ticks.y=ggplot2::element_blank()  #remove y axis ticks
      )+
      ggplot2::scale_fill_manual(values=c('brown','chartreuse4')) + ggplot2::ggtitle(paste0("Domain Enrichment of ", ifelse(x == "1", "(+)", "(-)"),
                                                                                   " PSI Exons"))

    pdf(paste0(output_location, ifelse(x == "1", "(+)", "(-)"), 'enrichmentPlot.pdf'))
    print(enrichmentPlot)
    dev.off()
    enrichmentPlot

  })
  print(eP)

  return(list(data=data,
              enrichmentPlots = eP))
}
