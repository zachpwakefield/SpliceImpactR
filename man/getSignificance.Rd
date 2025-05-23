% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats_models.R
\name{getSignificance}
\alias{getSignificance}
\title{Extract significance stats for each aRNAp}
\usage{
getSignificance(psi_data_sf, chosen_method)
}
\arguments{
\item{psi_data_sf}{a df with exons for signifcance calc}

\item{chosen_method}{stats method selected}
}
\value{
a dataframe with psi values, LR, cooks d
}
\description{
Extract significance stats for each aRNAp
}
\keyword{internal}
