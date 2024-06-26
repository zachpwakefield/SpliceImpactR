% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getPfam.R
\name{getPfam}
\alias{getPfam}
\title{get the pfam domains for each protein}
\usage{
getPfam(background, foreground, pdir, output_location, cores = 1)
}
\arguments{
\item{background}{output from getBackground}

\item{foreground}{output from getForeground}

\item{pdir}{directory of package}

\item{output_location}{location to make background directory}

\item{cores}{number of requested cores}
}
\value{
figures and dataframes with paired data
}
\description{
get the pfam domains for each protein
}
