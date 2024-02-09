getPaired <- function(foreground) {
  # Create a unique identifier for each exon
  foreground <- foreground %>%
    dplyr::mutate(exon_id = paste0(chr, ":", start, "-", stop))

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
      print("No paired exons found")
      return(NA)
    } else {

    # Combine all combinations into a single data frame
    exon_pairs_df <- do.call(rbind, exon_pairs_list) %>% dplyr::relocate(gene)
    # Combine all rows into a single data frame
    combined_rows_df <- do.call(rbind, combined_rows)
    combined_rows_df_expanded <- do.call(rbind, lapply(1:nrow(combined_rows_df), function(x) {
      rbind(foreground[foreground$exon_id == combined_rows_df$pos_exon_id[x],],
            foreground[foreground$exon_id == combined_rows_df$neg_exon_id[x],])
    }))

    # Use matchAlignType to identify protein alignment score and type
    c(proBed, pMatch, alignType) := matchAlignType(proBed = combined_rows_df_expanded, protCode = combined_rows_df_expanded$prot)

    # Make dataframe for plotting in ggplot2
    gdf_df <- data.frame(dens = as.numeric(pMatch), type = alignType)
    gdf_df2 <- gdf_df[gdf_df$type != "noPC",]


    # Alignment plot showing distribution of different type of exon swapping
    (gdf <- ggplot2::ggplot(gdf_df, ggplot2::aes(x = dens, fill = type)) +
        ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(count)/sum(ggplot2::after_stat(count))), colour = 1,
                                bins = 20) + ggplot2::geom_density(ggplot2::aes(y=.0005*ggplot2::after_stat(count)), color = 'black', fill = "coral2", bw = .1, alpha = .3) +
        ggplot2::scale_fill_manual(values=c('noPC' = "azure4", 'Match' = "#E69F00", 'onePC' = "#56B4E9", 'FrameShift' = "pink", 'PartialMatch' = "deeppink4")) +
        ggplot2::theme_classic() + ggplot2::xlab("Alignment Score") + ggplot2::ylab("Fraction"))

    # Alignment plot showing distribution of different type of exon swapping
    (gdf2 <- ggplot2::ggplot(gdf_df2, ggplot2::aes(x = dens, fill = type)) +
        ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(count)/sum(ggplot2::after_stat(count))), colour = 1,
                                bins = 20) + ggplot2::geom_density(ggplot2::aes(y=.0005*ggplot2::after_stat(count)), color = 'black', fill = "coral2", bw = .1, alpha = .3) +
        ggplot2::scale_fill_manual(values=c('noPC' = "azure4", 'Match' = "#E69F00", 'onePC' = "#56B4E9", 'FrameShift' = "pink", 'PartialMatch' = "deeppink4")) +
        ggplot2::theme_classic() + ggplot2::xlab("Alignment Score") + ggplot2::ylab("Fraction"))

    # Return the combined dataframe of exon pairs with their corresponding rows from the input dataframe
    # exon_pairs is minimal information
    # combined_rows is maximal information
    return(list(exon_pairs = exon_pairs_df, paired_proBed = combined_rows_df_expanded, gdf=gdf, gdf2=gdf2))
  }
}
}
