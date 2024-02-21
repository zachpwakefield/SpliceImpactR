make_dPsiPlot <- function(dpsi_df, pdir, num_thresh = 30) {
  # Load the gene reference conversion table
  hg38.conv <- readr::read_csv('/projectnb2/evolution/zwakefield/proteinImpacts/gencodev42_transcriptGeneProtein.csv')

  # Map Ensembl IDs in the lfc_df to HGNC symbols using the conversion table
  gene_id_to_name <- setNames(hg38.conv$gene_name, hg38.conv$gene_id)

  # Assuming the gene format is "GENEID.something"
  gene_ids <- sapply(strsplit(dpsi_df$gene, "\\."), `[`, 1)

  # If gene_id is in hg38.conv$gene_id, get the gene name, else keep the original gene_id
  dpsi_df$hgnc <- ifelse(gene_ids %in% names(gene_id_to_name),
                         gene_id_to_name[gene_ids],
                         dpsi_df$gene)
  # Prepare the data for labeling in the plot, sorting by absolute log fold change and adjusted p-value
  lab_thresh <- dpsi_df %>% dplyr::arrange(desc(abs(delta.psi)), p.adj)

  # Remove gene labels for points below the threshold (either by p-value or log fold change)
  dpsi_df$hgnc[-log(dpsi_df$p.adj) <= -log(lab_thresh$p.adj[num_thresh]) | abs(dpsi_df$delta.psi) < abs(lab_thresh$delta.psi[num_thresh])]  <- ""


  # Create the log fold change plot using ggplot2
  (deExons <- ggplot2::ggplot(dpsi_df,ggplot2:: aes(x = delta.psi, y = -log(p.adj), color = col, label = hgnc)) + ggplot2::geom_point(ggplot2::aes(shape = type), size = 2, color = dpsi_df$col) +
      ggplot2::theme_classic() + ggplot2::ylab("-Log2(FDR)") + ggplot2::xlab("Delta Psi")
    # +coord_cartesian(xlim = c(-100, 20))
    + ggplot2::geom_text(hjust=.2, vjust=0, size = 3)
  )
  print("delta psi plot done!") # Print completion message
  return(deExons) # Return the ggplot2 object for the plot
}
