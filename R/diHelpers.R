#' add colors to differential inclusion analysis for plotting
#'
#' @param de_df from differential inclusion functions
#' @param color_thresh the threshold to color by, default of .2
#' @return the input dataframe with color info added
#' @keywords internal
diColor <- function(de_df, color_thresh = .2) {
  de_df$col <- "#A7A9AC"
  de_df$col[de_df$delta.psi <= -(color_thresh) & de_df$p.adj < .05] <- 'brown'
  de_df$col[de_df$delta.psi >= color_thresh & de_df$p.adj < .05] <- 'chartreuse4'
  return(de = de_df)
}

#' filter using minimums for identifying events in control/test and max proportion of zeros allowed
#'
#' @param final_data from diHIT or diRMATS passed through diColor
#' @param nT number of test samples
#' @param nC number of control samples
#' @param min_prop_samples min proportion of samples needed to not be counted as outliers
#' @param max_zero_prop max number of zeros in to count an exon
#' @return a filtered di inclusion dataframe
#' @keywords internal
qualityFilter <- function(df, nT, nC, min_prop_samples = 0.5, max_zero_prop = 0.5) {
  df_filtered <- df[df$count_control >= min_prop_samples * nC &
                      df$count_test >= min_prop_samples * nT &
                      df$zero_count <= max_zero_prop * (nC+nT),]
  return(df_filtered)
}

#' filter using delta psi and fdr
#'
#' @param df from diHIT or diRMATS (or post filtering)
#' @param fdr fdr threshold
#' @param d.psi delta psi threshold
#' @return a filtered di inclusion dataframe
#' @keywords internal
significanceFilter <- function(df, fdr = .05, d.psi = .1) {
  df_filtered <- df[abs(df$delta.psi) >= d.psi & df$p.adj <= fdr,]
  return(df_filtered)
}
