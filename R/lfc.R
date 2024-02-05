## This function computes the log fold change (LFC) for differential expression analysis between control and experimental groups.
## It filters based on exon type, excludes outliers, and categorizes results based on significance and fold change direction.


lfc <- function(de_df, numCont, numExp, exon_type, cores = 8) {
  # Filter the dataframe for the specified exon type and relevant columns
  de_df <- de_df[de_df$type == exon_type,c(1:11, seq(12, ncol(de_df), by = 3), seq(13, ncol(de_df), by = 3), seq(14, ncol(de_df), by = 3))]
  samps <- colnames(de_df)

  lfc <- list()  # Initialize list to store log fold changes
  col <- list()  # Initialize list to store colors for visualization

  # Calculate log fold change for each gene
  lfc <- mclapply(1:length(de_df$gene), mc.cores = cores, function(i) {
    # Identify outliers in the dataset
    outlier <- which(unlist(lapply(samps, function(x) grepl(x, de_df$outlier[i]))))+11

    # Define control and experimental groups, excluding outliers
    cont <- c(12:(12-1+numCont))
    cont <- cont[!(cont %in% outlier)]
    exp <- c((numCont+11+1):(12-1+numCont+numExp))
    exp <- exp[!(exp %in% outlier)]


    # Compute log fold change based on available data, handling edge cases
    if (length(exp) == 0) {
      lfc_t <- log2((.01)/as.numeric(rowMeans(de_df[i,cont, drop = FALSE])+.01))
    }
    if (length(cont) == 0) {
      lfc_t <- log2((as.numeric(rowMeans(de_df[i,exp, drop = FALSE]))+.01)/.01)
    }
    if (length(exp) == 0 & length(cont) == 0) {
      lfc_t <- 0
    }
    if (length(exp) != 0 & length(cont) != 0) {
      lfc_t <- log2((as.numeric(rowMeans(de_df[i,exp, drop = FALSE]))+.01)/as.numeric(rowMeans(de_df[i,cont, drop = FALSE])+.01))
    }


    lfc_t # Return the computed log fold change
  })

  # Assign computed log fold changes to the dataframe
  de_df$lfc <- unlist(lfc)
  # Filter out non-significant and NA adjusted p-values
  de_df <- de_df[de_df$p.adj >= 0 & !is.na(de_df$p.adj),]

  # Assign colors based on log fold change and significance for visualization
  col <- lapply(1:length(de_df$gene), function(i) {
    if (de_df$lfc[i] <= -1.0 & de_df$p.adj[i] < .05) {
      'brown'  # Color for significant negative log fold change
    } else if (de_df$lfc[i] >= 1.0 & de_df$p.adj[i] < .05) {
      'chartreuse4'  # Color for significant positive log fold change
    } else {
      "#A7A9AC"  # Default color for non-significant changes
    }
  })



  de_df$col <- unlist(col) # Add color information to the dataframe

  # Count the number of outliers for each gene
  p.scDE <- de_df
  p.scDE$numOutliers <- str_count(p.scDE$outlier, "#")

  # Return the processed dataframe with differential expression analysis results
  return (de = p.scDE)
}
