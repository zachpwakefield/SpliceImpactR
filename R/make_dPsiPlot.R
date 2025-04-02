#' gets the tti interactions with various helper functions to do so
#'
#' @param dpsi_df pre filtering differential inclusion analysis df
#' @param gtf gtf from setupGTF
#' @param thresh threshold for differential inclusion
#'
#' @return differences between each tti pair and the overall results
#' @importFrom igraph graph_from_edgelist V make_ego_graph write_graph simplify E layout.fruchterman.reingold
#' @importFrom tidyr crossing
#' @importFrom ggplot2 ggplot aes geom_point theme_classic ylab xlab theme geom_text scale_fill_manual geom_bar
#' @importFrom ggpubr ggarrange
#' @importFrom dplyr select relocate
#' @keywords internal
make_dPsiPlot <- function(dpsi_df, gtf, thresh = .1) {

  # Map Ensembl IDs in the lfc_df to HGNC symbols using the conversion table
  gene_id_to_name <- setNames(gtf$geneName, gtf$geneID)

  gene_ids <- unlist(lapply(strsplit(dpsi_df$gene, "[.]"), '[[', 1))

  # If gene_id is in hg38.conv$gene_id, get the gene name, else keep the original gene_id
  dpsi_df$hgnc <- ifelse(gene_ids %in% names(gene_id_to_name),
                         gene_id_to_name[gene_ids],
                         dpsi_df$gene)
  # Prepare the data for labeling in the plot, sorting by absolute log fold change and adjusted p-value
  lab_thresh <- dpsi_df %>% dplyr::arrange(desc(abs(.data$delta.psi)), .data$p.adj)

  # Create the dpsi change plot using ggplot2
  deExons <- ggplot2::ggplot(dpsi_df, ggplot2::aes(x = .data$delta.psi, y = -log10(.data$p.adj), color = .data$col, label = .data$hgnc)) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$type), size = 2, color = dpsi_df$col) +
    ggplot2::theme_classic() + ggplot2::ylab("-Log2(FDR)") +
    ggplot2::xlab("Delta Psi") +
    ggplot2::theme(legend.position = "none")


  diE <- data.frame(val = c(sum(dpsi_df$delta.psi < -1*(thresh) & dpsi_df$p.adj < .05), sum(dpsi_df$delta.psi > thresh & dpsi_df$p.adj < .05)),
                    type = c("-", "+"))

  deExons_chart <- ggplot2::ggplot(diE, ggplot2::aes(x = .data$type, y = .data$val, fill = .data$type)) +
    ggplot2::geom_bar(stat="identity") + ggplot2::ylab("Count") + ggplot2::xlab("diExons") +
    ggplot2::theme_classic() + ggplot2::scale_fill_manual(values=c("brown", "chartreuse4"), breaks = c("-", "+")) +
    ggplot2::geom_text(aes(label = .data$val), vjust = 1.5, colour = "white", size = 5) +
    ggplot2::theme(legend.position = "none")


  comb_plot <- ggpubr::ggarrange(deExons_chart, deExons, widths = c(1, 2.5))
  return(comb_plot) # Return the ggplot2 object for the plot
}
