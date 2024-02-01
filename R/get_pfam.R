get_pfam <- function(background, foreground, pdir, output_location, cores) {
  pfam_hg38 <- read_tsv(paste0(pdir, '/protein_code_from_gencodev43_headerFix.txt.tsv'), col_names = F)
  pfam_hg38$transcriptID <- unlist(lapply(strsplit(pfam_hg38$X1, split = "#"),
                                          "[[", 1))
  pfam_hg38 <- pfam_hg38 %>% dplyr::relocate(transcriptID)
  pfam_hg38$geneID <- unlist(lapply(strsplit(unlist(lapply(strsplit(pfam_hg38$X1,
                                                                    split = "#"), "[[", 2)), split = ";"), "[[", 1))
  pfam_hg38 <- pfam_hg38 %>% dplyr::relocate(geneID, .after = transcriptID)

  fg_out <- do.call(rbind, parallel::mclapply(1:length(foreground$proBed$transcript), mc.cores = cores, function(i) {
    if (foreground$proBed$transcript[i] %in% pfam_hg38$transcriptID) {
      df <- pfam_hg38[pfam_hg38$transcriptID == foreground$proBed$transcript[i],]
      id <- paste0(foreground$proBed$transcript[i], "#",
                   foreground$proBed$gene[i], ";",
                   foreground$proBed$chr[i], ":",
                   foreground$proBed$start[i], "-",
                   foreground$proBed$stop[i], ";",
                   foreground$proBed$strand[i])
      df$X1 <- rep(id, length(df$X1))
      df
    }
  }))
  bg_out <- do.call(rbind, paralell::mclapply(1:length(background$proBed$transcript), mc.cores = cores, function(i) {
    if (background$proBed$transcript[i] %in% pfam_hg38$transcriptID) {
      df <- pfam_hg38[pfam_hg38$transcriptID == background$proBed$transcript[i],]
      id <- paste0(background$proBed$transcript[i], "#",
                   background$proBed$gene[i], ";",
                   background$proBed$chr[i], ":",
                   background$proBed$start[i], "-",
                   background$proBed$stop[i], ";",
                   background$proBed$strand[i])
      df$X1 <- rep(id, length(df$X1))
      df
    }
  }))

  write_tsv(fg_out, paste0(output_location, "fgoutFast.fa.tsv"))
  write_tsv(bg_out, paste0(output_location, "bgoutFast.fa.tsv"))
  return(list(fg_out = fg_out,
              bg_out = bg_out))
}
