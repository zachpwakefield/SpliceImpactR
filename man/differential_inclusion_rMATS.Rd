% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_inclusion_rMATS.R
\name{differential_inclusion_rMATS}
\alias{differential_inclusion_rMATS}
\title{Extract differentially included AS events from HIT Index data}
\usage{
differential_inclusion_rMATS(
  control_names,
  test_names,
  et,
  cores = 1,
  outlier_threshold = c("4/n", "1", "Inf")[1],
  minReads = 10
)
}
\arguments{
\item{control_names}{a vector of control_names}

\item{test_names}{a vector of test_names}

\item{et}{string of the exon_type}

\item{cores}{the number of cores requested}

\item{outlier_threshold}{the thresholding of the cooks distance, no outlier removal is "Inf"}
}
\value{
a dataframe with differential inclusion information
}
\description{
Extract differentially included AS events from HIT Index data
}
