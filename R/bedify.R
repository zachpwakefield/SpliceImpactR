bedifyForeground <- function(matched = matched, outname = outname, cores = 8) {
  tID <- matched$transcriptID
  eiID <- matched$input_id
  p.adj <- matched$p.adj
  delta.psi <- matched$delta.psi
  toBed <- list()
  toBed <- parallel::mclapply(1:length(tID), mc.cores = cores, function(i) {
    bed <- gtf[gtf$transcriptID == tID[i] & gtf$type == "exon",] %>% dplyr::select(chr, start, stop, transcriptID, geneID, strand)
    bed$eiID <- rep(eiID[i], length(gtf$geneID[gtf$transcriptID == tID[i] & gtf$type == "exon"]))
    bed$score <- 0
    bed$squish <- paste(bed$transcriptID, "#", bed$eiID, sep = "")
    bed$delta.psi <- rep(delta.psi[i], length(gtf$geneID[gtf$transcriptID == tID[i] & gtf$type == "exon"]))
    bed$p.adj <- rep(p.adj[i], length(gtf$geneID[gtf$transcriptID == tID[i] & gtf$type == "exon"]))

    colnames(bed) <- c('chrom','chromStart', 'chromEnd', 'transcriptID', "geneID", "strand", "eiID", "score", "name", "delta.psi", "p.adj")
    bed <- bed %>% dplyr::select(chrom, chromStart, chromEnd, name, score, strand, delta.psi, p.adj)
    if (unique(bed$strand) == "+") {
      bed$chromStart <- as.integer(bed$chromStart) - 1
    } else {
      bed$chromStart <- as.integer(bed$chromEnd) + 1
    }
    bed
  })
  toBed <- do.call(rbind, toBed)
  return(toBed)
}


bedifyBackground <- function(matched = matched, outname = outname, cores = 8) {
  tID <- matched$transcriptID
  eiID <- matched$input_id
  toBed <- list()
  toBed <- parallel::mclapply(1:length(tID), mc.cores = cores, function(i) {
    bed <- gtf[gtf$transcriptID == tID[i] & gtf$type == "exon",] %>% dplyr::select(chr, start, stop, transcriptID, geneID, strand)
    bed$eiID <- rep(eiID[i], length(gtf$geneID[gtf$transcriptID == tID[i] & gtf$type == "exon"]))
    bed$score <- 0
    bed$squish <- paste(bed$transcriptID, "#", bed$eiID, sep = "")

    colnames(bed) <- c('chrom','chromStart', 'chromEnd', 'transcriptID', "geneID", "strand", "eiID", "score", "name")
    bed <- bed %>% dplyr::select(chrom, chromStart, chromEnd, name, score, strand)
    if (unique(bed$strand) == "+") {
      bed$chromStart <- as.integer(bed$chromStart) - 1
    } else {
      bed$chromStart <- as.integer(bed$chromEnd) + 1
    }
    bed
  })
  toBed <- do.call(rbind, toBed)
  return(toBed)
}
