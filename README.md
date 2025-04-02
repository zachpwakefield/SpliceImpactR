# SpliceImpactR <img src="./inst/screenshot1.png" alt="Fiszbein Lab Logo" width="110" align="right"/> <img src="./inst/screenshot2.png" width="200" align="right"/>

SpliceImpactR is an R package designed for studying the impact of alternative splicing on protein structure and function. It provides tools for analyzing RNA-seq data to identify differentially included splicing events and predict their consequences on the resulting protein products. SpliceImpactR output involves identifying key changes in proteins at various levels: primary sequence, domain content, and transcript-transcript interactions.

The suite of funcitons is designed to anaylyze the consequences of AFE, ALE, SE, MXE, A5SS, RI, and A3SS, along with hybrid exons (HFE, HLE). SpliceImpactR is built to take output primarily from the [HIT Index](https://github.com/thepailab/HITindex) and [rMATS](https://github.com/Xinglab/rmats-turbo). 

HIT Index data outputs, such as .exon files, are necessary for much of the process

## Features
Identification of alternative splicing events from RNA-seq data.
Analysis of the potential impact of splicing events on protein structure.
Functional annotation of spliced isoforms to predict their biological impact.
Integration with existing bioinformatics tools and databases for comprehensive analysis.
Holistic analysis of how the use of different RNA processing events differs.

## Installation
You can install SpliceImpactR directly from GitHub using the devtools package. If you do not have devtools installed, you can install it first with:

```r
install.packages("devtools")
```
Then, to install the package:
```r
devtools::install("fiszbein-lab/SpliceImpactR")
```
Or 
``` r
BiocManager::install("fiszbein-lab/SpliceImpactR", version="devel")
```




# Usage
## Access gencode information
__SpliceImpactR__ requires the acceession of various genome annotations, accessed through biomaRt and directly through gencode, here we access the gencode files
```r
## To save time, we add a directory path to save the annotations
translations <- SpliceImpactR::getTranslations("/path/")

transcripts <- SpliceImpactR::getTranscripts("/path/")
```

## Access biomaRt information
Here we obtain further annotations through biomaRt
```r
## Getting access to biomaRt sometimes proves challenging, if errors in completing the ensembl call, run the following commented lines in the same chunk as the call
# new_config <- httr::config(ssl_verifypeer = FALSE)
# httr::set_config(new_config, override = FALSE)

ensembl <- biomaRt::useEnsembl(biomart = "ensembl",
                               dataset = "hsapiens_gene_ensembl")

## To save time, we add a directory path to save the annotations
biomart_data <- SpliceImpactR:::setupBiomart("/path/")

gtf <- SpliceImpactR::getAnnotation(biomart_data)
```

## Toy genome annotation loading
For the sake of this intro, we use toy versions (limited to a handful of genes) of the annotations for speed and ease
```r
# For purposes of these examples, data directory is in extdata, however this should be replaced with the directory path to  data
dataDirectory <- paste0(system.file("extdata", package="SpliceImpactR"), "/")

transcripts_sample <- list(transDF = readr::read_csv(paste0(dataDirectory, "transcripts_limited_transDF.csv"), 
                                                           show_col_types = FALSE),
                       c_trans = readr::read_lines(paste0(dataDirectory, "transcripts_limited_c_trans.csv")))

gtf_sample <- list(gtf = readr::read_csv(paste0(dataDirectory, "gtf_limited.csv"), 
                                                           show_col_types = FALSE),
            transcript_gtf = readr::read_csv(paste0(dataDirectory, "transcript_gtf_limited.csv"), 
                                                           show_col_types = FALSE),
            tgp_biomart = readr::read_csv(paste0(dataDirectory, "tgp_biomart_limited"), 
                                                           show_col_types = FALSE)
            )

translations_sample <- readr::read_lines(paste0(dataDirectory, "translations_limited.csv"))

biomart_data_sample <- list(ip = readr::read_csv(paste0(dataDirectory, "biomart_ip.csv"), 
                                                           show_col_types = FALSE),
                            code_regions = readr::read_csv(paste0(dataDirectory, "biomart_code_regions.csv"), 
                                                           show_col_types = FALSE),
                            pfam_exon_level = readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon_level.csv"), 
                                                           show_col_types = FALSE),
                            fsd_exon_data = readr::read_csv(paste0(dataDirectory, "biomart_data_sample.csv"), 
                                                           show_col_types = FALSE),
                            pfam_data = readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon.csv"), 
                                                           show_col_types = FALSE))

```


## Setting up data_df to identify conditions for comparison
test_group and control_group should be formatted as vectors of the path to the 
files from the HIT Index and/or rMATS -- without the ".{type}PSI" ending 
(eg: if the AFEPSI files is in /files/sample1.AFEPSI, the path given should 
be /files/sample1)
Rmats files should be renamed to be the same as the HIT Index file paths, just 
with the respective alternative splicing event: /files/sample1.SEPSI, 
/files/sample1.RIPSI
```r
# We are feeding the toy data in extdata/rawData here for computational speed
test_group <- paste0(dataDirectory, "rawData/", c("test1","test2", "test3"))
control_group <- paste0(dataDirectory, "rawData/", c("control1", "control2", 
                                                     "control3"))

# the data_df should be organized as such
data_df <- data.frame(
    sample_names = c(control_group, test_group),
    phenotype_names = c(
      rep("control", length(control_group)),
      rep("test", length(test_group))
     ),
   stringsAsFactors = FALSE
  )
```


## Basic Usage
The basic usage is able to be performed in one main call, fullASoutcome. 
This generates output from multiple alternative RNA processing event types 
Note that HIT Index output (eg: .exon file) is required for all the steps 
to properly run in the pipeline, even if just looking at rMATS output 
```r
# fullASoutcome performs the full pipeline on multiple alternative RNA 
# processing events, only performed on AFE, SE, HIT here
twoASfullRun <- fullASoutcome(as_types = c("AFE", "SE", "HIT"),
                              output_directory = NULL,
                              data_directory = dataDirectory,
                              data_df,
                              outlier_handle = "Inf",
                              cutoff = .1,
                              cores = 1,
                              bg_pre = NA,
                              tti_location = NULL,
                              initTTI = NULL,
                              mOverlap = .05,
                              plotAlignments = FALSE,
                              s_gtf = gtf_sample,
                              transcripts = transcripts_sample,
                              translations = translations_sample,
                              biomart_data = biomart_data_sample,
                              max_zero_prop = 1,
                              min_prop_samples = 0,
                              chosen_method = 'qbGLM')

```

## Integrated analysis
In order to perform the holistic analysis, we need multiple types of 
alternative RNA processing events with valid outputs -- in the toy data, 
we just provide valid AFE changes for speed
```r
fg_list <- list("AFE" = twoASfullRun$AFE$fg$proBed)
pfg_list <- list("AFE" = twoASfullRun$AFE$pfg$paired_proBed)
domain_list <- list("AFE" = twoASfullRun$AFE$gD$data)

integrated <- getIntegratedResults(fg_list,
                                   pfg_list,
                                   domain_list)

```

## Individual event analysis
Now we'll walkthrough how to perform more specific analysis.
To get output for a single alternative RNA processing event, we use getfxnlASoutcome
```r
# First we need to add 'utc' column to data_df, which is handled internally in fullASoutcome
data_df$utc <- "control"
data_df$utc[data_df$phenotype_names == unique(data_df$phenotype_names)[2]] <- "test"

# Then we can run the function for AFE
oneASrun <- getfxnlASoutcome(output_location = NULL,
                             test_group,
                             control_group,
                             data_df,
                             exon_type = "AFE",
                             cutoff = .1,
                             outlier_handle = "Inf",
                             cores = 1,
                             tti_location = NULL,
                             full_pipe = FALSE,
                             bg = NA,
                             mOverlap = .05,
                             gtf_sample,
                             plotAlignments = FALSE,
                             transcripts = transcripts_sample,
                             translations = translations_sample,
                             biomart_data = biomart_data_sample,
                             max_zero_prop = 1,
                             min_prop_samples = 0,
                             chosen_method = 'qbGLM',
                             initTTI = NULL)
```

Here we are able to access the same output as from the fullASoutcome run for
the individual event type

## Step-by-step analysis
Now we can look at what is actually happening within getfxnlASoutcome and do 
some step-by-step analysis
First, we are going to perform differential inclusion analysis
```r
diHIT <- differential_inclusion_HITindex(test_names = test_group, 
                                         control_names = control_group, 
                                         et = "AFE",
                                         outlier_threshold = 'Inf', 
                                         minReads = 10,
                                         min_prop_samples = 0, 
                                         chosen_method = 'qbGLM')
```

Then, we match to annotations and extract significant events in the same step 
with getForeground. This leaves us with a collection of vectors and dataframes, 
from different steps of this process. The final output is proBed
```r
fg <- getForeground(input=diHIT,
                    test_names = test_group,
                    control_names = control_group,
                    thresh = 0.1,
                    fdr=.05,
                    mOverlap=0.1,
                    cores=1,
                    exon_type="AFE",
                    output_location=NULL,
                    gtf=gtf_sample,
                    max_zero_prop = 1, 
                    min_prop_samples = 0,
                    translations = translations_sample)
```

After getForeground, we use getPaired to extract 'exon swapping events', where 
conditions strongly prefer different alternative RNA processing events. This 
leaves us primarily with paired_proBed and exon_pairs, two representations 
of the paired events.
```r
pfg <- getPaired(foreground = fg$proBed, 
                 et = "AFE", 
                 nucleotides = transcripts_sample,
                 output_location = NULL, 
                 newGTF = gtf_sample, 
                 saveAlignments = FALSE,
                 exon_data = biomart_data_sample$fsd_exon_data)
```

For the next steps, we need to identify the set of transcripts identified in 
the samples, not restricterd to significant changes. HIT Index output of .exon 
files are necessary for this step. Depending on how many samples with how many 
exons, this step can take some time... 
With the toy data given, it is very quick.
```r
bg <- getBackground(input=c(control_group, test_group),
                    mOverlap = 0.1,
                    cores = 1,
                    exon_type = 'AFE',
                    output_location = NULL,
                    gtf=gtf_sample,
                    translations_sample)
```


Then, pfam domains are matched to foreground and background transcripts...
```r
pfam <- getPfam(foreground = fg,
                background = bg,
                pdir = system.file(package="SpliceImpactR"),
                cores = 1,
                output_location = NULL,
                biomart_data = biomart_data_sample)
```


Using this data info, we can get domain change data along with global 
enrichment of domain changes if any
```r
gD <- getDomainData(fg = fg,
                    bg = bg,
                    pfg=pfg,
                    cores = 1,
                    pfam = pfam,
                    output_location = NULL,
                    fdr_use = .05,
                    min_sample_success = 2,
                    engine = "Pfam",
                    topViz = 15)
```


For the final step of the overall pipeline, we investigate how the changing 
domains may rewire protein-protein interactions using 3did's set of 
domain-domain interactions. We need to initialize this here -- which may be a 
timely function on the full set of exons (given the toy biomart_data, 
this is quick).
getTTI saves plots if desired and outputs the differences between interaction 
networks (along with the networks if desired)
```r
initTTI <- init_ddi(pdir = dataDirectory,
                    output_location = NULL,
                    ppidm_class = "Gold_Standard",
                    cores = 1,
                    removeDups = TRUE,
                    pfam_data = biomart_data_sample$pfam_data)

tti <- getTTI(paired_foreground = pfg$paired_proBed,
              background = bg$proBed,
              steps=1,
              max_vertices_for_viz = 2000,
              fdr = .05,
              ppidm_class = "Gold_Standard",
              write_igraphs = FALSE,
              output_location = NULL,
              tti_location = NULL,
              tgp = gtf_sample$tgp_biomart,
              init_edgelist = initTTI$edgelist)
```


There are also a handful of adjacent methods, useful for overall probing of the 
data, but which aren't involved in the primary flow of the pipeline.
```r
# getOverviewComparison looks at overall number of events and events 
# per gene between the two conditions -- along with the distribution of PSI 
# values between the samples
initial_comparison <- getOverviewComparison(data_df = data_df, 
                                            exon_type = "AFE", 
                                            output_location = NULL)

# getProximalShift is used to probe AFE and ALE for any global shifts in the 
# proximal or distal directions
proxPlot <- getProximalShift(type = "AFE",
                             exon_pairs = pfg$exon_pairs, 
                             ep_supp = pfg$paired_proBed, 
                             output_location = NULL)

# getLengthComparison idenitifes any signficant shifts in leng  th of protein 
# usage due to significantly used alternative RNA processing events
length_comparison <- getLengthComparison(data_df = data_df, 
                                         paired_df = pfg$paired_proBed, 
                                         output_location = NULL)

# getHitCompare uses a wilcox rank sum test to identify significant changes in 
# the HTI Index values for each exon
hitCompare <- getHitCompare(data_df = data_df,
                            output_location = NULL,
                            threshold = .25)
```


## Contributing
Contributions to SpliceImpactR are welcome, including bug reports, feature requests, and pull requests. Please see CONTRIBUTING.md for guidelines on how to contribute.

## Support
If you encounter any problems or have suggestions, please file an issue on the GitHub issue tracker. Or contact zachpw@bu.edu

##Citation
If you use SpliceImpactR in your research, please cite:

```bibtex
Zachary Wakefield
SpliceImpactR: Analyzing the Impact of Alternative Splicing on Protein Structure and Function
2024
https://github.com/fiszbein-lab/SpliceImpactR
```

