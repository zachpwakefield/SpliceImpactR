make_lfcPlot <- function(lfc_df, location = system.file(package="domainEnrichment"), num_thresh = 10) {
  hg38.conv <- readRDS(paste0(location, "/hg38_geneRef_conv.RDS"))
  lfc_df$hgnc <- unlist(lapply(lfc_df$gene, function(x) {ifelse(unlist(lapply(strsplit(x, split = "[.]"), "[[", 1)) %in% hg38.conv$ens,
                                                                unique(hg38.conv$hgnc[hg38.conv$ens == unlist(lapply(strsplit(x, split = "[.]"), "[[", 1))]), x)

  }))
  hg38.conv$hgnc[hg38.conv$ens %in% unlist(lapply(strsplit(lfc_df$gene, split = "[.]"), "[[", 1))]


  lab_thresh <- lfc_df %>% dplyr::arrange(desc(abs(lfc)), p.adj)
  lfc_df$hgnc[-log(lfc_df$p.adj) <= -log(lab_thresh$p.adj[num_thresh]) | abs(lfc_df$lfc) < abs(lab_thresh$lfc[num_thresh])]  <- ""

  (deExons <- ggplot2::ggplot(lfc_df,ggplot2:: aes(x = lfc, y = -log(p.adj), color = col, label = hgnc)) + ggplot2::geom_point(ggplot2::aes(shape = type), size = 2, color = lfc_df$col) +
      ggplot2::theme_classic() + ggplot2::ylab("-Log Adj P Value") + ggplot2::xlab("Log2FC")
    # +coord_cartesian(xlim = c(-100, 20))
    + ggplot2::geom_text(hjust=.2, vjust=0, size = 4)
  )
  print("lfc done!")
  return(deExons)
}
