get_pfam <- function(background, foreground, pdir, output_location) {
  pfam_hg38 <- read_tsv(paste0(pdir, '/protein_code_from_gencodev43_headerFix.txt.tsv'), col_names = F)
  pfam_hg38$transcriptID <- unlist(lapply(strsplit(pfam_hg38$X1, split = "#"),
                                          "[[", 1))
  pfam_hg38 <- pfam_hg38 %>% dplyr::relocate(transcriptID)
  pfam_hg38$geneID <- unlist(lapply(strsplit(unlist(lapply(strsplit(pfam_hg38$X1,
                                                                    split = "#"), "[[", 2)), split = ";"), "[[", 1))
  pfam_hg38 <- pfam_hg38 %>% dplyr::relocate(geneID, .after = transcriptID)


  fg_out <- do.call(rbind, mclapply(1:length(foreground$matched$tot_matched$input_id), mc.cores = 8, function(c) {
    cx <- foreground$matched$tot_matched[c,]
    cxdf <- pfam_hg38[pfam_hg38$transcriptID %in% cx$transcriptID,]
    loc <- unlist(lapply(strsplit(cx$input_id, split = ";"), "[[", 2))
    id <- unlist(lapply(strsplit(cxdf$X1, split = ";"), "[[", 1))
    strand <- unlist(lapply(strsplit(cxdf$X1, split = ";"), "[[", 2))
    cxdf$X1 <- paste0(id, ";", loc, ";", strand)
    cxdf
  }))

  bg_out <- do.call(rbind, mclapply(1:length(background$matched$tot_matched$input_id), mc.cores = 8, function(c) {
    cx <- background$matched$tot_matched[c,]
    cxdf <- pfam_hg38[pfam_hg38$transcriptID %in% cx$transcriptID,]
    loc <- unlist(lapply(strsplit(cx$input_id, split = ";"), "[[", 2))
    id <- unlist(lapply(strsplit(cxdf$X1, split = ";"), "[[", 1))
    strand <- unlist(lapply(strsplit(cxdf$X1, split = ";"), "[[", 2))
    cxdf$X1 <- paste0(id, ";", loc, ";", strand)
    cxdf
  }))

  write_tsv(fg_out, paste0(output_location, "fgoutFast.fa.tsv"))
  write_tsv(bg_out, paste0(output_location, "bgoutFast.fa.tsv"))
  return(list(fg_out = fg_out,
              bg_out = bg_out))
}
