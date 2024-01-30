get_pfam <- function(background, output_location) {
  pfam_hg38 <- read_tsv('/projectnb2/evolution/zwakefield/proteinImpacts/protein_code_from_gencodev43_headerFix.txt.tsv', col_names = F)
  pfam_hg38$transcriptID <- unlist(lapply(strsplit(pfam_hg38$X1, split = "#"),
                                          "[[", 1))
  pfam_hg38 <- pfam_hg38 %>% dplyr::relocate(transcriptID)
  pfam_hg38$geneID <- unlist(lapply(strsplit(unlist(lapply(strsplit(pfam_hg38$X1,
                                                                    split = "#"), "[[", 2)), split = ";"), "[[", 1))
  pfam_hg38 <- pfam_hg38 %>% dplyr::relocate(geneID, .after = transcriptID)
  ip_out <- pfam_hg38[pfam_hg38$transcriptID %in% background$matched$tot_matched$transcriptID,]
  write_tsv(ip_out, paste0(output_location, "bgoutFast.fa.tsv"))
  return(ip_out)
}
