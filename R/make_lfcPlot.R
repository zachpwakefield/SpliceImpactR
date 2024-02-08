make_lfcPlot <- function(lfc_df, pdir, num_thresh = 30) {
  # Load the gene reference conversion table
  hg38.conv <- readRDS(paste0(pdir, "/hg38_geneRef_conv.RDS"))

  # Map Ensembl IDs in the lfc_df to HGNC symbols using the conversion table

  # Extract Ensembl ID before the dot
  # Return HGNC symbol if found
  # Return original Ensembl ID if HGNC symbol not found
  lfc_df$hgnc <- unlist(lapply(lfc_df$gene, function(x) {ifelse(unlist(lapply(strsplit(x, split = "[.]"), "[[", 1)) %in% hg38.conv$ens,
                                                                unique(hg38.conv$hgnc[hg38.conv$ens == unlist(lapply(strsplit(x, split = "[.]"), "[[", 1))]), x)

  }))

  # Prepare the data for labeling in the plot, sorting by absolute log fold change and adjusted p-value
  lab_thresh <- lfc_df %>% dplyr::arrange(desc(abs(lfc)), p.adj)

  # Remove gene labels for points below the threshold (either by p-value or log fold change)
  lfc_df$hgnc[-log(lfc_df$p.adj) <= -log(lab_thresh$p.adj[num_thresh]) | abs(lfc_df$lfc) < abs(lab_thresh$lfc[num_thresh])]  <- ""


  # Create the log fold change plot using ggplot2
  (deExons <- ggplot2::ggplot(lfc_df,ggplot2:: aes(x = lfc, y = -log(p.adj), color = col, label = hgnc)) + ggplot2::geom_point(ggplot2::aes(shape = type), size = 2, color = lfc_df$col) +
      ggplot2::theme_classic() + ggplot2::ylab("-Log Adj P Value") + ggplot2::xlab("Log2FC")
    # +coord_cartesian(xlim = c(-100, 20))
    + ggplot2::geom_text(hjust=.2, vjust=0, size = 4)
  )
  print("lfc plot done!") # Print completion message
  return(deExons) # Return the ggplot2 object for the plot
}
