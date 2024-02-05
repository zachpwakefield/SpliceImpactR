getPaired <- function(foreground) {
  proBed <- foreground$proBed[order(foreground$proBed$gene),]
  paired_proBed <- proBed[proBed$gene %in% names(table(proBed$gene))[as.numeric(table(proBed$gene)) > 1],]
  moreThanTwo <- names(table(paired_proBed$gene))[table(paired_proBed$gene) > 2]
  if (length(moreThanTwo) > 0) {
    lapply(moreThanTwo, function(x) {
      sg <- paired_proBed[paired_proBed$gene %in% x,]
      if (sum(sg$prot != "none") == 2) {
        paired_proBed <- rbind(paired_proBed[paired_proBed$gene!= x,], sg[sg$prot != "none",])
      } else {paired_proBed <- rbind(paired_proBed[paired_proBed$gene!= x,], sg[sample(length(sg[,1]), 2),])}
    })
  }
  c(proBed, pMatch, alignType) := matchAlignType(proBed = paired_proBed, protCode = paired_proBed$prot)
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
  return(list(paired_proBed = proBed,
              gdf = gdf,
              gdf2 = gdf2))
}
