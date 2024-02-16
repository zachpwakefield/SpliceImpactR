# get_c_trans
# input: directory containing f_protein_code_from_gencodev43.txt, default should be package directory or: system.file(package="domainEnrichment")
# output: f_protein_code_from_gencodev43.txt lines

get_c_trans <- function(location) {
  return(readr::read_lines(paste(location,"/f_protein_code_from_gencodev43.txt",sep="")))
}

# get_gtf
# input: directory containing all the gtf files, default should be package directory or: system.file(package="domainEnrichment")
# output: cg csv

get_gtf <- function(location) {
  cg <- do.call(rbind, lapply(paste0(location, paste0(c("/gene", "/transcript", "/single_exon", "/first", paste0("/internal", 1:5), "/last"), "_gtf.csv")), function(x) {
    read.csv(x, header = T)
  })) %>% dplyr::arrange(rownum)
  return(cg)
}
