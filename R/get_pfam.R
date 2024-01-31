get_pfam <- function(background, foreground, pdir, output_location) {
  pfam_hg38 <- read_tsv(paste0(pdir, '/protein_code_from_gencodev43_headerFix.txt.tsv'), col_names = F)
  pfam_hg38$transcriptID <- unlist(lapply(strsplit(pfam_hg38$X1, split = "#"),
                                          "[[", 1))
  pfam_hg38 <- pfam_hg38 %>% dplyr::relocate(transcriptID)
  pfam_hg38$geneID <- unlist(lapply(strsplit(unlist(lapply(strsplit(pfam_hg38$X1,
                                                                    split = "#"), "[[", 2)), split = ";"), "[[", 1))
  pfam_hg38 <- pfam_hg38 %>% dplyr::relocate(geneID, .after = transcriptID)


  fg_out <- pfam_hg38[pfam_hg38$transcriptID %in% foreground$matched$tot_matched$transcriptID,]
  bg_out <- pfam_hg38[pfam_hg38$transcriptID %in% background$matched$tot_matched$transcriptID,]

  write_tsv(fg_out, paste0(output_location, "fgoutFast.fa.tsv"))
  write_tsv(bg_out, paste0(output_location, "bgoutFast.fa.tsv"))
  return(list(fg_out = fg_out,
              bg_out = bg_out))
}
