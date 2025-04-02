#' get paired AS events / exons and identify similarities/differences through sequence alignment + identify frame shifts
#'
#' @param foreground proBed from getForeground
#' @param nucleotides made from setup get_c_nucs
#' @param et exon_type being investigated
#' @param output_location location to make background directory
#' @param cores number of requested cores
#' @param saveAlignments bool to save alignment pdfs
#' @param newGTF from setup get_gtf
#' @param exon_data from biomart_data
#' @return figures and dataframes with paired data
#' @importFrom dplyr %>% mutate filter rename left_join
#' @importFrom ggplot2 ggplot aes geom_histogram after_stat geom_density scale_fill_manual theme_classic xlab ylab
#' @importFrom ggpubr ggarrange
#'
#' @examples
#'
#' pdir <- system.file("extdata", package="SpliceImpactR")
#' dataDirectory <- paste0(pdir, "/")
#' test_group <- paste0(dataDirectory, "rawData/", c("test1","test2", "test3"))
#' control_group <- paste0(dataDirectory, "rawData/", c("control1", "control2", "control3"))
#' transDF <- readr::read_csv(paste0(dataDirectory, "transcripts_limited_transDF.csv"))
#' c_trans <- readr::read_lines(paste0(dataDirectory, "transcripts_limited_c_trans.csv"))
#'
#' transcripts_sample <- list(transDF = transDF,
#'                            c_trans = c_trans)
#'
#' gtf_sample <- list(gtf = readr::read_csv(paste0(dataDirectory, "gtf_limited.csv")),
#'             transcript_gtf = readr::read_csv(paste0(dataDirectory, "transcript_gtf_limited.csv")))
#' translations_sample <- readr::read_lines(paste0(dataDirectory, "translations_limited.csv"))
#' biomart_data_sample <- readr::read_csv(paste0(dataDirectory, "biomart_data_sample.csv"))
#'
#'
#' result <- differential_inclusion_HITindex(test_names = test_group,
#'                                           control_names = control_group,
#'                                           et = "AFE",
#'                                           outlier_threshold = "Inf",
#'                                           minReads = 10,
#'                                           min_prop_samples = 0,
#'                                           chosen_method = "qbGLM"
#'                                           )
#'
#' fg <- getForeground(input = result,
#'                             test_names = test_group,
#'                             control_names = control_group,
#'                             thresh = .1,
#'                             fdr = .05,
#'                             mOverlap = .1,
#'                             exon_type = "AFE",
#'                             output_location = NULL,
#'                             cores = 1,
#'                             gtf = gtf_sample,
#'                             max_zero_prop = 1,
#'                             min_prop_samples = 0,
#'                             translations = translations_sample)
#' library(msa)
#' pfg <- getPaired(foreground = fg$proBed,
#'           et = "AFE",
#'           nucleotides = transcripts_sample,
#'           newGTF = gtf_sample,
#'           cores = 1,
#'           output_location = NULL,
#'           saveAlignments = FALSE,
#'           exon_data = biomart_data_sample)
#'
#' @export
getPaired <- function(foreground, et, nucleotides, newGTF, cores = 1, output_location = NULL, saveAlignments = TRUE, exon_data) {
  # Create a unique identifier for each exon
  foreground <- foreground %>%
    dplyr::mutate(exon_id = paste0(chr, ":", start, "-", stop))

  # incorporate pair number into gene name for the pairing process
  foreground$gene <- paste0(foreground$gene,
                            ifelse(unique(foreground$add_inf) == "none", "#",
                               paste0("#", sub(".*;", "", foreground$add_inf))))

  # Split the data into positive and negative delta.psi, and filter out genes with only positive or only negative delta.psi
  pos_exons <- foreground %>% dplyr::filter(delta.psi > 0)
  neg_exons <- foreground %>% dplyr::filter(delta.psi < 0)

  valid_genes <- intersect(pos_exons$gene, neg_exons$gene)  # Genes with both positive and negative delta.psi

  # Filter out genes with only one type of delta.psi (edge case handled here)
  pos_exons <- pos_exons %>% dplyr::filter(gene %in% valid_genes)
  neg_exons <- neg_exons %>% dplyr::filter(gene %in% valid_genes)

  # Initialize an empty list to store exon pairs for each valid gene
  exon_pairs_list <- list()
  combined_rows <- list()

  # Loop through each valid gene to find all possible exon pairs
  for (gene_select in valid_genes) {
    pos <- pos_exons %>% dplyr::filter(gene == gene_select)
    neg <- neg_exons %>% dplyr::filter(gene == gene_select)


    # Create all possible combinations of positive and negative exons for the gene
    combos <- expand.grid(pos_exon_id = pos$exon_id, neg_exon_id = neg$exon_id)
    combos$gene <- gene_select  # Add a column for the gene

    # Join back to get the delta.psi values for the positive and negative exons
    combos <- combos %>%
      dplyr::left_join(pos %>% dplyr::select(exon_id, delta.psi), by = c("pos_exon_id" = "exon_id")) %>%
      dplyr::left_join(neg %>% dplyr::select(exon_id, delta.psi), by = c("neg_exon_id" = "exon_id"))

    # Rename delta.psi columns
    combos <- combos %>%
      dplyr::rename(pos_delta.psi = delta.psi.x, neg_delta.psi = delta.psi.y)



    # Store the combinations in the list, using the gene as the list name
    exon_pairs_list[[gene_select]] <- combos

    # Join combos back to the original df to get the full rows for each exon pair
    combined <- combos %>%
      dplyr::left_join(pos, by = c("pos_exon_id" = "exon_id")) %>%
      dplyr::left_join(neg, by = c("neg_exon_id" = "exon_id"), suffix = c(".pos", ".neg"))

    # Store combined rows for each gene
    combined_rows[[gene_select]] <- combined
  }

    if (length(exon_pairs_list) == 0) {
      message("No paired exons found")
      return(NA)
    } else {

    # Combine all combinations into a single data frame
    exon_pairs_df <- do.call(rbind, exon_pairs_list) %>% dplyr::relocate(gene)
    # Combine all rows into a single data frame
    combined_rows_df <- do.call(rbind, combined_rows)
    combined_rows_df_expanded <- do.call(rbind, lapply(seq_len(nrow(combined_rows_df)), function(x) {
      rbind(foreground[foreground$exon_id == combined_rows_df$pos_exon_id[x] &
                         foreground$add_inf == combined_rows_df$add_inf.pos[x] &
                         foreground$gene == combined_rows_df$gene[x],],
            foreground[foreground$exon_id == combined_rows_df$neg_exon_id[x] &
                         foreground$add_inf == combined_rows_df$add_inf.neg[x] &
                         foreground$gene == combined_rows_df$gene[x],])
    }))

    # Use matchAlignType to identify protein alignment score and type
    placeholder_tri <- matchAlignType(proBed = combined_rows_df_expanded, protCode = combined_rows_df_expanded$prot, nucleotides = nucleotides, output_location = output_location,
                                                   saveAlignments = saveAlignments, exon_type = et, newgtf = newGTF$gtf, exon_data = exon_data)
    proBed <- placeholder_tri[[1]]
    pMatch <- placeholder_tri[[2]]
    alignType <- placeholder_tri[[3]]

    combined_rows_df_expanded$pMatch <- as.numeric(pMatch)
    combined_rows_df_expanded$alignType <- alignType

    exon_pairs_df$pMatch <- as.numeric(pMatch)[seq(1, length(pMatch), by=2)]
    exon_pairs_df$alignType <- alignType[seq(1, length(pMatch), by=2)]

    if (et == "HFE") {
      pair_trans <- unlist(lapply(seq(1, nrow(combined_rows_df_expanded), by=2), function(x) {
        paste0(combined_rows_df_expanded$transcript[x], ';', combined_rows_df_expanded$transcript[x+1])
      }))
      in_hfe_pair <- pair_trans %in% unlist(newGTF$hybrid_first_extract_transcripts)
      combined_rows_df_expanded <- combined_rows_df_expanded[rep(in_hfe_pair, each = 2),]
      exon_pairs_df <- exon_pairs_df[in_hfe_pair,]
    } else if (et == "HLE") {
      pair_trans <- unlist(lapply(seq(1, nrow(combined_rows_df_expanded), by=2), function(x) {
        paste0(combined_rows_df_expanded$transcript[x], ';', combined_rows_df_expanded$transcript[x+1])
      }))
      in_hle_pair <- pair_trans %in% unlist(newGTF$hybrid_last_extract_transcripts)
      combined_rows_df_expanded <- combined_rows_df_expanded[rep(in_hle_pair, each = 2),]
      exon_pairs_df <- exon_pairs_df[in_hle_pair,]
    }


    # Make dataframe for plotting in ggplot2
    gdf_df <- data.frame(dens = as.numeric(exon_pairs_df$pMatch), type = exon_pairs_df$alignType)
    gdf_df <- gdf_df[gdf_df$type != "noPC" & gdf_df$type != "onePC",]


    dfProp <- data.frame(vals = c(sum(gdf_df == 'noPC'), sum(gdf_df == 'onePC'), sum(gdf_df == 'Match'), sum(gdf_df == 'FrameShift'), sum(gdf_df == 'PartialMatch')),
                         type = c("noPC", "onePC", "Match", "FrameShift", "PartialMatch"),
                         col = c("Pair Type", "Pair Type", "Pair Type", "Pair Type", "Pair Type"))

    propCoding <- ggplot2::ggplot(dfProp, ggplot2::aes(fill=.data$type, x = .data$col, y = .data$vals)) +
        ggplot2::geom_bar(position="stack", stat="identity") +
      ggplot2::scale_fill_manual(values=c('noPC' = "azure4", 'Match' = "#E69F00", 'onePC' = "#56B4E9", 'FrameShift' = "pink", 'PartialMatch' = "deeppink4")) +
      ggplot2::theme_classic() + ggplot2::xlab("") + ggplot2::ylab("Count")


    gdf <- ggplot2::ggplot(gdf_df, ggplot2::aes(x = dens, fill = type)) +
        ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(count)/sum(ggplot2::after_stat(count))), colour = 1,
                                bins = 20) +
        ggplot2::scale_fill_manual(values=c('noPC' = "azure4", 'Match' = "#E69F00", 'onePC' = "#56B4E9", 'FrameShift' = "pink", 'PartialMatch' = "deeppink4")) +
        ggplot2::theme_classic() + ggplot2::xlab("Alignment Score") + ggplot2::ylab("Fraction")

    gdf_comp <- ggpubr::ggarrange(propCoding, gdf, nrow = 1, widths = c(1, 2),
                                   common.legend = TRUE)
    # remove pair number from gene name after the pairing process
    exon_pairs_df$gene <- unlist(lapply(strsplit(exon_pairs_df$gene, split = "#"), "[[", 1))
    combined_rows_df_expanded$gene <- unlist(lapply(strsplit(combined_rows_df_expanded$gene, split = "#"), "[[", 1))

    if (!is.null(output_location)) {
    system(paste0("mkdir ", output_location, "pairedOutput/"))
    readr::write_csv(exon_pairs_df, paste0(output_location, "pairedOutput/", "exon_pairs.csv"))
    readr::write_csv(combined_rows_df_expanded, paste0(output_location, "pairedOutput/", "paired_combined_rows.csv"))

    pdf(paste0(output_location, "pairedOutput/", "primarySequencePlot.pdf"))
    print(gdf_comp)
    dev.off()
    }

    return(list(exon_pairs = exon_pairs_df,
                paired_proBed = combined_rows_df_expanded,
                gdf=gdf_comp))
    # Return the combined dataframe of exon pairs with their corresponding rows from the input dataframe
    # exon_pairs is minimal information
    # combined_rows is maximal information
  }
}
