#' setting up gtf with exon labels, first, internal, last, hybrid, etc
#' @importFrom dplyr filter group_by mutate ungroup inner_join select case_when
#' @importFrom biomaRt useEnsembl getBM
#' @return annotated gtf and various hybrid info. also gene to transcript to protein naming dataframe
#' @keywords internal

setup_gtf <- function() {
  ensembl <- biomaRt::useEnsembl(biomart = "ensembl",
                                 dataset = "hsapiens_gene_ensembl",
                                 mirror = "useast")
  attributes <- c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand", "ensembl_gene_id", "ensembl_transcript_id", "gene_biotype",
                  "transcript_biotype", "external_gene_name", "rank", "ensembl_exon_id")
  exon_data <- biomaRt::getBM(attributes = attributes, mart = ensembl)

  exon_data_exon_label <- exon_data %>% dplyr::filter(chromosome_name %in% c(1:23, "X", "Y")) %>%
    dplyr::group_by(ensembl_gene_id, ensembl_transcript_id) %>%
    dplyr::mutate(exon_label = dplyr::case_when(
      rank == min(rank) ~ "first",
      rank == max(rank) ~ "last",
      TRUE ~ "internal"
    )) %>%
    dplyr::ungroup()


  first_exons <- exon_data_exon_label %>% dplyr::filter(exon_label == "first")
  internal_exons <- exon_data_exon_label %>% dplyr::filter(exon_label == "internal")
  last_exons <- exon_data_exon_label %>% dplyr::filter(exon_label == "last")

  # Join the first and internal exons to find overlaps within the same gene
  hybrid_first_exons <- dplyr::inner_join(first_exons, internal_exons, by = "ensembl_gene_id", suffix = c("_first", "_internal")) %>%
    dplyr::filter(chromosome_name_first == chromosome_name_internal &
             ensembl_transcript_id_first != ensembl_transcript_id_internal &
             ((exon_chrom_start_first <= exon_chrom_end_internal & exon_chrom_end_first >= exon_chrom_start_internal) |
                (exon_chrom_start_internal <= exon_chrom_end_first & exon_chrom_end_internal >= exon_chrom_start_first)))

  # Select relevant columns to display the hybrid exons
  hybrid_first_exons <- hybrid_first_exons %>% dplyr::select(ensembl_gene_id, ensembl_transcript_id_first, exon_chrom_start_first, exon_chrom_end_first, ensembl_exon_id_first,
                                                             ensembl_transcript_id_internal, exon_chrom_start_internal, exon_chrom_end_internal, ensembl_exon_id_internal)

  # Join the first and internal exons to find overlaps within the same gene
  hybrid_last_exons <- dplyr::inner_join(last_exons, internal_exons, by = "ensembl_gene_id", suffix = c("_last", "_internal")) %>%
    dplyr::filter(chromosome_name_last == chromosome_name_internal &
             ensembl_transcript_id_last != ensembl_transcript_id_internal &
             ((exon_chrom_start_last <= exon_chrom_end_internal & exon_chrom_end_last >= exon_chrom_start_internal) |
                (exon_chrom_start_internal <= exon_chrom_end_last & exon_chrom_end_internal >= exon_chrom_start_last)))

  # Select relevant columns to display the hybrid exons
  hybrid_last_exons <- hybrid_last_exons %>% dplyr::select(ensembl_gene_id, ensembl_transcript_id_last, exon_chrom_start_last, exon_chrom_end_last, ensembl_exon_id_last,
                                                           ensembl_transcript_id_internal, exon_chrom_start_internal, exon_chrom_end_internal, ensembl_exon_id_internal)



  attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name", "external_transcript_name", "ensembl_peptide_id", "chromosome_name")
  ptg_init <- biomaRt::getBM(attributes = attributes, mart = ensembl, uniqueRows = TRUE)
  ptg_conversion <- ptg_init[ptg_init$chromosome_name %in% c(1:23, "X", "Y"),]
  ptg_conversion$chromosome_name <- NULL
  ptg_conversion$external_gene_name[ptg_conversion$external_gene_name == ""] <- ptg_conversion$ensembl_gene_id[ptg_conversion$external_gene_name == ""]
  ptg_conversion$external_transcript_name[ptg_conversion$external_transcript_name == ""] <- ptg_conversion$ensembl_gene_id[ptg_conversion$external_transcript_name == ""]
  ptg_conversion$ensembl_peptide_id[ptg_conversion$ensembl_peptide_id == ""] <- "none"

  tgp_biomart <- data.frame(
    gene_id = ptg_conversion$ensembl_gene_id,
    transcript_id = ptg_conversion$ensembl_transcript_id,
    gene_name = ptg_conversion$external_gene_name,
    transcript_name = ptg_conversion$external_transcript_name,
    protein_id = ptg_conversion$ensembl_peptide_id)

  gtf <- data.frame(
    geneID = exon_data_exon_label$ensembl_gene_id,
    transcriptID = exon_data_exon_label$ensembl_transcript_id,
    geneName = exon_data_exon_label$external_gene_name,
    exonID = exon_data_exon_label$ensembl_exon_id,
    classification = exon_data_exon_label$exon_label,
    rownum = 1:nrow(exon_data_exon_label),
    chr = paste0("chr", exon_data_exon_label$chromosome_name),
    start = exon_data_exon_label$exon_chrom_start,
    stop = exon_data_exon_label$exon_chrom_end,
    strand = exon_data_exon_label$strand,
    type = rep("exon", nrow(exon_data_exon_label)),
    gpc = exon_data_exon_label$gene_biotype,
    tpc = exon_data_exon_label$transcript_biotype
  )
  gtf <- gtf[gtf$gpc == "protein_coding",]


  transcript_attributes <- c("ensembl_transcript_id", "chromosome_name", "strand", "transcript_start", "transcript_end")
  transcript_data <- biomaRt::getBM(attributes = transcript_attributes, mart = ensembl, values = c(1:23, "X", "Y"), filters = 'chromosome_name')
  merge_transcript_exon <- gtf[gtf$classification == "first",c("transcriptID", "rownum")]
  merge_transcript_exon <- merge_transcript_exon[!duplicated(merge_transcript_exon),]
  transcript_gtf <- dplyr::left_join(transcript_data, merge_transcript_exon, by = join_by('ensembl_transcript_id' == 'transcriptID'))
  transcript_gtf <- transcript_gtf[transcript_gtf$ensembl_transcript_id %in% gtf$transcriptID,]
  transcript_gtf$chromosome_name <- paste0('chr', transcript_gtf$chromosome_name)
  colnames(transcript_gtf) <- c('transcriptID', 'chr', 'strand', 'start', 'stop', 'rownum')
  transcript_gtf <- transcript_gtf %>% dplyr::mutate(strand = ifelse(strand == -1, "-", "+"))
  gtf <- gtf %>% dplyr::mutate(strand = ifelse(strand == -1, "-", "+"))



  hleList <- unique(c(paste0(hybrid_last_exons$ensembl_transcript_id_last, ';', hybrid_last_exons$ensembl_transcript_id_internal),
                      paste0(hybrid_last_exons$ensembl_transcript_id_internal, ';', hybrid_last_exons$ensembl_transcript_id_last)))
  hfeList <- unique(c(paste0(hybrid_first_exons$ensembl_transcript_id_first, ';', hybrid_first_exons$ensembl_transcript_id_internal),
                      paste0(hybrid_first_exons$ensembl_transcript_id_internal, ';', hybrid_first_exons$ensembl_transcript_id_first)))


  return(list(gtf = gtf,
              transcript_gtf = transcript_gtf,
              hybrid_last_extract = hybrid_last_exons,
              hybrid_first_extract = hybrid_first_exons,
              hybrid_first_extract_transcripts = hfeList,
              hybrid_last_extract_transcripts = hleList,
              tgp_biomart = tgp_biomart))
}
# setup_gtf <- function(gtf_location, cores = 1) {
#   gtf <- rtracklayer::import(gtf_location)
#   gtf_df <- as.data.frame(gtf)
#   pcgtf <- gtf_df[gtf_df$gene_type == "protein_coding",]
#   pcgtf$classification <- ""
#   pcgtf$rowname <- 1:nrow(pcgtf)
#
#   transcript_indices <- c(which(pcgtf$type == 'transcript'), length(pcgtf$type))
#   gene_indices <- c(which(pcgtf$type == 'gene'), length(pcgtf$type))
#   strands <- pcgtf$strand[transcript_indices]
#   pcgtf$type[is.na(pcgtf$type)] <- ""
#
#   overlap <- function(start1, end1, start2, end2) {
#     start_min <- min(start1, end1)
#     end_max <- max(start1, end1)
#     temp_start2 <- pmin(start2, end2)
#     temp_end2 <- pmax(start2, end2)
#     return(pmax(start_min, start2) <= pmin(end_max, end2))
#   }
#
#   partial_overlap <- function(start1, end1, start2, end2) {
#     start_min <- min(start1, end1)
#     end_max <- max(start1, end1)
#     temp_start2 <- pmin(start2, end2)
#     temp_end2 <- pmax(start2, end2)
#     return(start_min < temp_start2 | end_max < temp_end2)
#   }
#   firstExonID <- function(min_gtf) {
#     if ("start_codon" %in% min_gtf$type) {
#       rn <- which(overlap(min_gtf$start[min_gtf$type == "start_codon"], min_gtf$end[min_gtf$type == "start_codon"],
#                           min_gtf$start, min_gtf$end) & min_gtf$type == "exon")
#     } else if ("CDS" %in% min_gtf$type) {
#       rn <- which(min_gtf$type == "exon" & overlap(min_gtf$start[min_gtf$type == "CDS"], min_gtf$end[min_gtf$type == "CDS"],
#                                                    min_gtf$start, min_gtf$end))[1]
#     } else if ("UTR" %in% min_gtf$type) {
#       rn <- which(min_gtf$type == "exon" & partial_overlap(min_gtf$start[min_gtf$type == "UTR"], min_gtf$end[min_gtf$type == "UTR"],
#                                                            min_gtf$start, min_gtf$end))[1]
#     } else {
#       rn <- which(min_gtf$type == "exon")[1]
#     }
#     return(rn)
#   }
#
#   lastExonID <- function(min_gtf) {
#     if ("stop_codon" %in% min_gtf$type) {
#       rn <- which(overlap(min_gtf$start[min_gtf$type == "stop_codon"], min_gtf$end[min_gtf$type == "stop_codon"],
#                           min_gtf$start, min_gtf$end) & min_gtf$type == "exon")
#     } else if ("CDS" %in% min_gtf$type) {
#       rn <- rev(which(min_gtf$type == "exon" & overlap(min_gtf$start[min_gtf$type == "CDS"], min_gtf$end[min_gtf$type == "CDS"],
#                                                        min_gtf$start, min_gtf$end)))[1]
#     } else if ("UTR" %in% min_gtf$type) {
#       rn <- rev(which(min_gtf$type == "exon" & partial_overlap(min_gtf$start[min_gtf$type == "UTR"], min_gtf$end[min_gtf$type == "UTR"],
#                                                                min_gtf$start, min_gtf$end)))[1]
#     } else {
#       rn <- rev(which(min_gtf$type == "exon"))[1]
#     }
#     return(rn)
#   }
#
#   fl_exons <- parallel::mclapply(1:(length(transcript_indices)-1), mc.cores = cores, function(x) {
#     min_gtf <- pcgtf[((transcript_indices[x]+1):(transcript_indices[x+1]-2)),]
#     if (sum(min_gtf$type == "exon") == 1) {
#       list(min_gtf$rowname[which(min_gtf$type == "exon")])
#     } else {
#       s <- sort(min_gtf$rowname[firstExonID(min_gtf)], decreasing = ifelse(unique(min_gtf$strand)[1] == "+", TRUE, FALSE))[1]
#       e <- sort(min_gtf$rowname[lastExonID(min_gtf)], decreasing = ifelse(unique(min_gtf$strand)[1] == "+", FALSE, TRUE))[1]
#       i <- min_gtf$rowname[!min_gtf$rowname %in% c(s, e) & min_gtf$type == "exon" & min_gtf$rowname %in% seq(s+1, e-1)]
#       utr <- min_gtf$rowname[!(min_gtf$rowname %in% c(i, s, e)) & min_gtf$type == "exon"]
#       list(s, e, i, utr)
#     }
#   })
#   single_exon_transcripts <- unlist(lapply(fl_exons, function(y) {
#     if (length(y) == 1) {
#       y
#     }
#   }))
#   multi_exon_transcripts <- lapply(fl_exons, function(y) {
#     if (length(y) != 1) {
#       y
#     }
#   })
#   pcgtf$classification[pcgtf$rowname %in% single_exon_transcripts] <- "single_exon"
#   pcgtf$classification[pcgtf$rowname %in% unlist(lapply(multi_exon_transcripts, "[[", 3))] <- "internal"
#   pcgtf$classification[pcgtf$rowname %in% unlist(lapply(multi_exon_transcripts, "[[", 1))] <- "first"
#   pcgtf$classification[pcgtf$rowname %in% unlist(lapply(multi_exon_transcripts, "[[", 2))] <- "last"
#   pcgtf$classification[pcgtf$rowname %in% unlist(lapply(multi_exon_transcripts, "[[", 4))] <- "UTR"
#   pcgtf$classification[pcgtf$rowname %in% intersect(unlist(lapply(multi_exon_transcripts, "[[", 1)), unlist(lapply(multi_exon_transcripts, "[[", 2)))] <- "single_exon"
#   pcgtf$transcriptID <- unlist(lapply(strsplit(pcgtf$transcript_id, split = "[.]"), "[[", 1))
#
#   # hybrid_first_extract <- unlist(unlist(parallel::mclapply(1:(length(gene_indices)-1), mc.cores = 8, function(x) {
#   hybrid_first_extract_duo <- unlist(unlist(parallel::mclapply(1:(length(gene_indices)-1), mc.cores = cores, function(x) {
#     tr_in_gene <- pcgtf[(gene_indices[x]+1):(gene_indices[x+1]-1),]
#     if (sum(tr_in_gene$type == 'transcript') > 1) {
#       lapply(unique(tr_in_gene$transcript_id), function(y) {
#         tr_start <- tr_in_gene[tr_in_gene$transcript_id == y & tr_in_gene$classification == "first",]
#         tr_internal <- tr_in_gene[tr_in_gene$transcript_id != y & tr_in_gene$classification == "internal",]
#         if (nrow(tr_start) > 0 & nrow(tr_internal) > 0) {
#           over <- overlap(tr_start$start, tr_start$end, tr_internal$start, tr_internal$end)
#         } else {return(NA)}
#         if (sum(over) >= 1) {
#           return(list(lapply(tr_internal$rowname[over], function(z) c(unique(tr_start$rowname), z)),
#                       lapply(tr_internal$transcriptID[over], function(z) c(paste0(unique(tr_start$transcriptID), ';', z),
#                                                                            paste0(z, ";", unique(tr_start$transcriptID)))))
#           )
#         } else {return(NA)}
#       })
#     } else {return(NA)}
#   }), recursive = FALSE), recursive = FALSE)
#
#   interM_first <- hybrid_first_extract_duo[!is.na(hybrid_first_extract_duo)]
#   hybrid_first_extract <- unlist(interM_first[seq(1, length(interM_first), by = 2)], recursive = FALSE)
#   hybrid_first_extract_transcripts <- unlist(interM_first[seq(2, length(interM_first), by = 2)], recursive = FALSE)
#
#   hybrid_last_extract_duo <- unlist(unlist(parallel::mclapply(1:(length(gene_indices)-1), mc.cores = cores, function(x) {
#     tr_in_gene <- pcgtf[(gene_indices[x]+1):(gene_indices[x+1]-1),]
#     if (sum(tr_in_gene$type == 'transcript') > 1) {
#       lapply(unique(tr_in_gene$transcript_id), function(y) {
#         tr_stop <- tr_in_gene[tr_in_gene$transcript_id == y & tr_in_gene$classification == "last",]
#         tr_internal <- tr_in_gene[tr_in_gene$transcript_id != y & tr_in_gene$classification == "internal",]
#         if (nrow(tr_stop) > 0 & nrow(tr_internal) > 0) {
#           over <- overlap(tr_stop$start, tr_stop$end, tr_internal$start, tr_internal$end)
#         } else {return(NA)}
#         if (sum(over) >= 1) {
#           return(list(lapply(tr_internal$rowname[over], function(z) c(unique(tr_stop$rowname), z)),
#                       lapply(tr_internal$transcriptID[over], function(z) c(paste0(unique(tr_stop$transcriptID), ';', z),
#                                                                            paste0(z, ";", unique(tr_stop$transcriptID)))))
#           )
#         } else {return(NA)}
#       })
#     } else {return(NA)}
#   }), recursive = FALSE), recursive = FALSE)
#
#   interM_last <- hybrid_last_extract_duo[!is.na(hybrid_last_extract_duo)]
#   hybrid_last_extract <- unlist(interM_last[seq(1, length(interM_last), by = 2)], recursive = FALSE)
#   hybrid_last_extract_transcripts <- unlist(interM_last[seq(2, length(interM_last), by = 2)], recursive = FALSE)
#
#
#   pcgtf$HFE <- "non-hybrid"
#   pcgtf$HLE <- "non-hybrid"
#
#   pcgtf$HFE[pcgtf$rowname %in% unique(unlist(lapply(hybrid_first_extract, "[[", 1)))] <- "HFE_F"
#   pcgtf$HFE[pcgtf$rowname %in% unique(unlist(lapply(hybrid_first_extract, "[[", 2)))] <- "HFE_I"
#
#   pcgtf$HLE[pcgtf$rowname %in% unique(unlist(lapply(hybrid_last_extract, "[[", 1)))] <- "HLE_F"
#   pcgtf$HLE[pcgtf$rowname %in% unique(unlist(lapply(hybrid_last_extract, "[[", 2)))] <- "HLE_I"
#
#   ngtf <- data.frame(geneID = unlist(lapply(strsplit(pcgtf$gene_id, split = "[.]"), "[[", 1)),
#                      transcriptID = pcgtf$transcriptID,
#                      transcriptName = pcgtf$transcript_name,
#                      geneName = pcgtf$gene_name,
#                      exonID = unlist(lapply(strsplit(pcgtf$exon_id, split = "[.]"), "[[", 1)),
#                      classification = pcgtf$classification,
#                      rownum = pcgtf$rowname,
#                      chr = pcgtf$seqnames,
#                      start = pcgtf$start,
#                      stop = pcgtf$end,
#                      strand = pcgtf$strand,
#                      method = pcgtf$source,
#                      type = as.character(pcgtf$type),
#                      gpc = as.character(pcgtf$gene_type),
#                      tpc = as.character(pcgtf$transcript_type))
#
#   ngtf$transcriptID[is.na(ngtf$transcriptID)] <- "gene"
#   ngtf$transcriptName[is.na(ngtf$transcriptName)] <- "gene"
#   ngtf$classification[ngtf$classification == ""] <- ngtf$type[ngtf$classification == ""]
#   gtf <- ngtf[ngtf$classification %in% c("gene", "transcript", "first", "internal", "last"),]
#   return(list(gtf = gtf,
#               hybrid_first_extract = hybrid_first_extract,
#               hybrid_last_extract = hybrid_last_extract,
#               hybrid_first_extract_transcripts = hybrid_first_extract_transcripts,
#               hybrid_last_extract_transcripts = hybrid_last_extract_transcripts))
# }
