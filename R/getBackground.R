#' getBackground set of transcripts in samples
#'
#' @param input the path containing the .exon file
#' @param mOverlap overlap to identify a match to annotation
#' @param exon_type placeholder for other functions
#' @param pdir location of the package
#' @param gtf full output from setup_gtf
#' @param output_location location to make background directory
#' @param translations from getTranslations
#' @return matched : matched transcripts dataframe, bed : bed file of the matched transcripts
#' proBed : output for further functions with protein code and protein info,
#' proFast : fasta file of proteins identified in proBed
#' @importFrom dplyr arrange first left_join group_by summarise
#' @importFrom tidyr separate
#' @importFrom stringr str_extract
#' @export
getBackground <- function(input, mOverlap, cores, exon_type, pdir, output_location, gtf, translations) {
    ## extract all first exons and create combined data.frame with gene, location
    files <- paste(input, unlist(lapply(input, function(x) list.files(x)[grep('[.]exon', list.files(x))])), sep = "")

    first_exons <- unique(unlist(lapply(files, function(x) {
        in_file <- read.delim(x)
        paste(in_file$gene, ';', in_file$exon, ';',  in_file$strand, sep = "")})))
        first_exons <- first_exons[grepl('[-]', first_exons) & grepl(';', first_exons)]
        redExon <- data.frame(geneR = unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 1)), split = '[.]'), "[[", 1)),
            chr = unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 1)), split = ":"), "[[", 1)),
            start = as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 1)), split = ":"), "[[", 2))),
            stop = as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 2))
                 )
    )

    ## Remove duplicate rows
    redExon <- redExon[!duplicated(redExon),]

    message("background exons loaded, matching now, this may take a little bit...")
    matched <- getTranscriptBackground(gtf$gtf, redExon, exon_type, mOverlap, cores, gtf$transcript_gtf)

    message("Bedifying and adding protein info to background...")
    bed <- bedifyBackground(matched, outname = output_location, cores, gtf$gtf)

    ## extract unique transcript names as trans and all transcript names as possT
    trans <- stringr::str_extract(unique(bed$name), "^[^#]*")

    # For all names
    possT <- stringr::str_extract(bed$name, "^[^#]*")

    ## Find annotated proteins for transcripts if possible
    protCode <- translations[match(trans, translations) + 1]
    protCode[is.na(protCode)]<- "none"

    merger <- data.frame(name = unique(bed$name),
                       prot = protCode)

    ## Make dataframe proBed for output of matched transcripts with protein code
    bed_summary <- bed %>%
        dplyr::group_by(.data$name) %>%
        dplyr::summarise(strand = dplyr::first(.data$strand),
                        .groups = 'drop')

    proBed <- bed_summary %>% dplyr::left_join(merger, by='name') %>%
        tidyr::separate(.data$name, c("transcript", "id"), "#") %>%
        tidyr::separate(.data$id, c("gene", "chr"), ";") %>%
        tidyr::separate(.data$chr, c("chr", "coords"), ':') %>%
        tidyr::separate(.data$coords, c("start", "stop"), '-')

    ## Create the FASTA headers using vectorized paste function
    fasta_headers <- paste0(">", proBed$transcript, "#", proBed$gene, ";", proBed$chr, ":", proBed$start, "-", proBed$stop, ";", proBed$strand)

    ## Assuming 'proBed$prot' contains the protein sequences
    fasta_sequences <- proBed$prot

    ## Interleave headers and sequences
    proFast <- paste(rbind(fasta_headers, fasta_sequences))

    system(paste0("mkdir ", output_location, "Background/"))
    write_csv(proBed, paste0(output_location, "Background/", "bgoutBed.csv"))
    write_lines(proFast, paste0(output_location, "Background/", "bgoutFast.fa"))
    write_csv(matched,  paste0(output_location, "Background/", "bgmatched.csv"))
    write_csv(bed,  paste0(output_location, "Background/", "bgbed.csv"))

    message("Background complete...")
    return(list(matched = matched,
                    bed = bed,
                proBed = proBed,
                proFast = proFast))
}
