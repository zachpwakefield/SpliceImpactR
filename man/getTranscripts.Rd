% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setup.R
\name{getTranscripts}
\alias{getTranscripts}
\title{setting up transcriptions}
\usage{
getTranscripts(transcripts_location)
}
\arguments{
\item{transcripts_location}{location of transcripts file / where to save the transcript file if not already imported}
}
\value{
modified nucleotide code fasta
}
\description{
setting up transcriptions
}
\examples{

pdir <- system.file("extdata", package = "SpliceImpactR")
dataDirectory <- paste0(pdir, "/")
transcripts <- getTranscripts(transcripts_location = dataDirectory)

}
