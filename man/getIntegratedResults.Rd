% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getIntegrative.R
\name{getIntegratedResults}
\alias{getIntegratedResults}
\title{Perform various integrative analyses on multiple aRNAp results}
\usage{
getIntegratedResults(fg_list, pfg_list, domain_list)
}
\arguments{
\item{fg_list}{a list of the fglfc.csv output from getForeground, names are the event types}

\item{pfg_list}{a list of the paired_output from getPaired, names are the event types}

\item{domain_list}{a list of the data output from getDomainData, names are the event types}
}
\value{
a dataframe with summary stats for each event type supplied along with various plots
}
\description{
Perform various integrative analyses on multiple aRNAp results
}
\examples{

pdir <- system.file("extdata", package="SpliceImpactR")
dataDirectory <- paste0(pdir, "/")
test_group <- paste0(dataDirectory, "rawData/", c("test1","test2", "test3"))
control_group <- paste0(dataDirectory, "rawData/", c("control1", "control2", "control3"))
data_df <- data.frame(
    sample_names = c(control_group, test_group),
    phenotype_names = c(
      rep("control", length(control_group)),
      rep("test", length(test_group))
     ),
   stringsAsFactors = FALSE
  )

transDF <- readr::read_csv(paste0(dataDirectory, "transcripts_limited_transDF.csv"))
c_trans <- readr::read_lines(paste0(dataDirectory, "transcripts_limited_c_trans.csv"))

transcripts_sample <- list(transDF = transDF,
                           c_trans = c_trans)

gtf_sample <- list(gtf = readr::read_csv(paste0(dataDirectory, "gtf_limited.csv")),
            transcript_gtf = readr::read_csv(paste0(dataDirectory, "transcript_gtf_limited.csv")),
            tgp_biomart = readr::read_csv(paste0(dataDirectory, "tgp_biomart_limited"))
            )
translations_sample <- readr::read_lines(paste0(dataDirectory, "translations_limited.csv"))

ip <- readr::read_csv(paste0(dataDirectory, "biomart_ip.csv"))
code_regions <- readr::read_csv(paste0(dataDirectory, "biomart_code_regions.csv"))
pfam_exon_level <- readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon_level.csv"))
fsd_exon_data <- readr::read_csv(paste0(dataDirectory, "biomart_data_sample.csv"))
pfam_data = readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon.csv"))
biomart_data_sample <- list(ip = ip,
                     code_regions = code_regions,
                     fsd_exon_data = fsd_exon_data,
                     pfam_exon_level = pfam_exon_level,
                     pfam_data = pfam_data)

initTTI <- init_ddi(pdir = dataDirectory,
                    output_location = NULL,
                    ppidm_class = c("Gold_Standard", "Gold", "Silver", "Bronze")[1],
                    removeDups = TRUE,
                    cores = 1,
                    pfam_data = biomart_data_sample$pfam_data)
twoASfullRun <- fullASoutcome(as_types = c("AFE", "SE", "HIT"),
                              output_directory = NULL,
                              data_directory = dataDirectory,
                              data_df,
                              outlier_handle = "Inf",
                              cutoff = .1,
                              cores = 1,
                              bg_pre = NA,
                              tti_location = NULL,
                              initTTI = initTTI,
                              mOverlap = .05,
                              s_gtf = gtf_sample,
                              plotAlignments = FALSE,
                              transcripts = transcripts_sample,
                              translations = translations_sample,
                              biomart_data = biomart_data_sample,
                              max_zero_prop = 1,
                              min_prop_samples = 0,
                              chosen_method = 'qbGLM')

fg_list <- list("AFE" = twoASfullRun$AFE$fg$proBed)
pfg_list <- list("AFE" = twoASfullRun$AFE$pfg$paired_proBed)
domain_list <- list("AFE" = twoASfullRun$AFE$gD$data)

integrated <- getIntegratedResults(fg_list,
                                   pfg_list,
                                   domain_list)
}
