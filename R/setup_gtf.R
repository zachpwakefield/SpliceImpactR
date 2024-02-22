setup_gtf <- function(gtf_location) {
  gtf <- rtracklayer::import(gtf_location)
  gtf_df=as.data.frame(gtf)
  pcgtf <- gtf_df[gtf_df$gene_type == "protein_coding",]
  pcgtf$classification <- ""
  pcgtf$rowname <- 1:nrow(pcgtf)

  transcript_indices <- c(which(pcgtf$type == 'transcript'), length(pcgtf$type))
  gene_indices <- c(which(pcgtf$type == 'gene'), length(pcgtf$type))
  strands <- pcgtf$strand[transcript_indices]
  pcgtf$type[is.na(pcgtf$type)] <- ""

  overlap <- function(start1, end1, start2, end2) {
    start_min <- min(start1, end1)
    end_max <- max(start1, end1)
    temp_start2 <- pmin(start2, end2)
    temp_end2 <- pmax(start2, end2)
    return(pmax(start_min, start2) <= pmin(end_max, end2))
  }


  fl_exons <- parallel::mclapply(1:(length(transcript_indices)-1), mc.cores = 8, function(x) {
    min_gtf <- pcgtf[((transcript_indices[x]+1):(transcript_indices[x+1]-1)),]
    s <- min_gtf$rowname[ifelse("start_codon" %in% min_gtf$type,
                                which(overlap(min_gtf$start[min_gtf$type == "start_codon"], min_gtf$end[min_gtf$type == "start_codon"],
                                              min_gtf$start, min_gtf$end) & min_gtf$type == "exon"),
                                ifelse(strands[x] == "+", which(min_gtf$type == "exon")[1],
                                       rev(which(min_gtf$type == "exon"))[1]))]
    e <- min_gtf$rowname[ifelse("stop_codon" %in% min_gtf$type,
                                which(overlap(min_gtf$start[min_gtf$type == "stop_codon"], min_gtf$end[min_gtf$type == "stop_codon"],
                                              min_gtf$start, min_gtf$end) & min_gtf$type == "exon"),
                                ifelse(strands[x] == "+", rev(which(min_gtf$type == "exon"))[1],
                                       which(min_gtf$type == "exon")[1]))]
    i <- min_gtf$rowname[!min_gtf$rowname %in% c(s, e) & min_gtf$type == "exon"]
    list(s, e, i)
  })
  pcgtf$classification[pcgtf$rowname %in% unlist(lapply(fl_exons, "[[", 3))] <- "internal"
  pcgtf$classification[pcgtf$rowname %in% unlist(lapply(fl_exons, "[[", 1))] <- "first"
  pcgtf$classification[pcgtf$rowname %in% unlist(lapply(fl_exons, "[[", 2))] <- "last"
  pcgtf$classification[pcgtf$rowname %in% intersect(unlist(lapply(fl_exons, "[[", 1)), unlist(lapply(fl_exons, "[[", 2)))] <- "first/last"

  hybrid_first_extract <- unlist(unlist(parallel::mclapply(1:(length(gene_indices)-1), mc.cores = 8, function(x) {
    tr_in_gene <- pcgtf[(gene_indices[x]+1):(gene_indices[x+1]-1),]
    if (sum(tr_in_gene$type == 'transcript') > 1) {
      lapply(unique(tr_in_gene$transcript_id), function(y) {
        tr_start <- tr_in_gene[(tr_in_gene$transcript_id %in% y) & tr_in_gene$classification %in% c("first/last", "first"),]
        tr_internal <- tr_in_gene[!(tr_in_gene$transcript_id %in% y) & tr_in_gene$classification == "internal",]
        over <- overlap(tr_start$start, tr_start$end, tr_internal$start, tr_internal$end)
        if (sum(over) >= 1) {
          return(lapply(tr_internal$rowname[over], function(z) c(unique(tr_start$rowname), z)))
        } else {return(NA)}
      })
    } else {return(NA)}
  }), recursive = F), recursive = F)

  hybrid_last_extract <- unlist(unlist(parallel::mclapply(1:(length(gene_indices)-1), mc.cores = 8, function(x) {
    tr_in_gene <- pcgtf[(gene_indices[x]+1):(gene_indices[x+1]-1),]
    if (sum(tr_in_gene$type == 'transcript') > 1) {
      lapply(unique(tr_in_gene$transcript_id), function(y) {
        tr_start <- tr_in_gene[(tr_in_gene$transcript_id %in% y) & tr_in_gene$classification %in% c("first/last", "last"),]
        tr_internal <- tr_in_gene[!(tr_in_gene$transcript_id %in% y) & tr_in_gene$classification == "internal",]
        over <- overlap(tr_start$start, tr_start$end, tr_internal$start, tr_internal$end)
        if (sum(over) >= 1) {
          return(lapply(tr_internal$rowname[over], function(z) c(unique(tr_start$rowname), z)))
        } else {return(NA)}
      })
    } else {return(NA)}
  }), recursive = F), recursive = F)

  hybrid_last_extract <- hybrid_last_extract[!is.na(hybrid_last_extract)]
  hybrid_first_extract <- hybrid_first_extract[!is.na(hybrid_first_extract)]


  pcgtf$HFE <- "non-hybrid"
  pcgtf$HLE <- "non-hybrid"

  pcgtf$HFE[pcgtf$rowname %in% unique(unlist(lapply(hybrid_first_extract, "[[", 1)))] <- "HFE_F"
  pcgtf$HFE[pcgtf$rowname %in% unique(unlist(lapply(hybrid_first_extract, "[[", 2)))] <- "HFE_I"

  pcgtf$HLE[pcgtf$rowname %in% unique(unlist(lapply(hybrid_last_extract, "[[", 1)))] <- "HLE_F"
  pcgtf$HLE[pcgtf$rowname %in% unique(unlist(lapply(hybrid_last_extract, "[[", 2)))] <- "HLE_I"

  ngtf <- data.frame(geneID = unlist(lapply(strsplit(pcgtf$gene_id, split = "[.]"), "[[", 1)),
                     transcriptID = unlist(lapply(strsplit(pcgtf$transcript_id, split = "[.]"), "[[", 1)),
                     transcriptName = pcgtf$transcript_name,
                     geneName = pcgtf$gene_name,
                     classification = pcgtf$classification,
                     rownum = pcgtf$rowname,
                     chr = pcgtf$seqnames,
                     start = pcgtf$start,
                     stop = pcgtf$end,
                     strand = pcgtf$strand,
                     method = pcgtf$source,
                     type = as.character(pcgtf$type),
                     gpc = as.character(pcgtf$gene_type),
                     tpc = as.character(pcgtf$transcript_type))


  ngtf$transcriptID[is.na(ngtf$transcriptID)] <- "gene"
  ngtf$transcriptName[is.na(ngtf$transcriptName)] <- "gene"
  ngtf$classification[ngtf$classification == ""] <- ngtf$type[ngtf$classification == ""]
  gtf <- ngtf[ngtf$classification %in% c("gene", "transcript", "first", "internal", "last"),]
  return(gtf)
}
