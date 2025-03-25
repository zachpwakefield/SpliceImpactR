#' calculates jaccard between sequences
#' @return jaccard index and alt jacc length
#' @keywords internal
calculate_jaccard <- function(start1, stop1, start2, stop2) {
  intersection_length <- pmax(0, pmin(stop1, stop2) - pmax(start1, start2) + 1)
  union_length <- (stop1 - start1 + 1) + (stop2 - start2 + 1) - intersection_length
  jaccard_index <- intersection_length / union_length
  length_jacc <- intersection_length/(abs(start2-stop2)+1)
  return(list(jaccard_index = jaccard_index,
              length_jacc = length_jacc))
}
