diff_info <- function(de_df, numCont, numExp, exon_type, cores = 8, test_names, control_names, color_thresh = .2) {

  de_df$p.adj[is.na(de_df$p.adj)] <- 1
  # Filter out non-significant and NA adjusted p-values
  de_df <- de_df[de_df$p.adj >= 0 & !is.na(de_df$p.adj),]

  # Assign colors based on log fold change and significance for visualization
  col <- lapply(1:length(de_df$gene), function(i) {
    if (de_df$delta.psi[i] <= -(color_thresh) & de_df$p.adj[i] < .05) {
      'brown'  # Color for significant negative log fold change
    } else if (de_df$delta.psi[i] >= (color_thresh) & de_df$p.adj[i] < .05) {
      'chartreuse4'  # Color for significant positive log fold change
    } else {
      "#A7A9AC"  # Default color for non-significant changes
    }
  })


  de_df$col <- unlist(col) # Add color information to the dataframe

  # Return the processed dataframe with differential expression analysis results
  return (de = de_df)
}
