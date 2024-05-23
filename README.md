# SpliceImpactR

<p align="center">
<img src="./inst/Screenshot 2024-02-08 at 3.27.19 PM.png" width="180" height="180">
<img src="./inst/Screenshot 2024-02-12 at 1.53.27 PM.png" width="300" height="180">
</p>


by Zachary Wakefield


SpliceImpactR is an R package designed for studying the impact of alternative splicing on protein structure and function. It provides tools for analyzing RNA-seq data to identify differentially included splicing events and predict their consequences on the resulting protein products. SpliceImpactR output involves identifying key changes in proteins at various levels: primary sequence, domain content, and transcript-transcript interactions.

The suite of funcitons is designed to anaylyze the consequences of AFE, ALE, SE, MXE, A5SS, RI, and A3SS, along with hybrid exons (HFE, HLE). SpliceImpactR is built to take output from the [HIT Index](https://github.com/thepailab/HITindex) and [rMATS](https://github.com/Xinglab/rmats-turbo). 

SpliceImpactR first identifies differentially included exons across the input phenotypes. This is performed differently for HIT Index output (AFE/ALE/HFE/HLE) and rMATS output (SE/MXE/A5SS/A3SS).

## Features
Identification of alternative splicing events from RNA-seq data.
Analysis of the potential impact of splicing events on protein structure.
Functional annotation of spliced isoforms to predict their biological impact.
Integration with existing bioinformatics tools and databases for comprehensive analysis.

## Installation
You can install SpliceImpactR directly from GitHub using the devtools package. If you do not have devtools installed, you can install it first with:

```r
install.packages("devtools")
```
Then, to install the package:
```r
devtools::install("zachpwakefield/SpliceImpactR")
```

## Basic Usage

```r
library(SpliceImpactR)
```

SpliceImpactR requires files organized in a particular way:
```
sample_name/
  sample_name.AFEPSI
  sample_name.ALEPSI
  sample_name.exon
  sample_name.SEPSI
  ...
```

In order to accomplish this, there is the function organizeSamples(). This works if your files are currently organized as such:
```
sample_name/
  rmats/
    sample_name.AFEPSI
    sample_name.ALEPSI
    sample_name.exon
  hit/
    SE.MATS.JC.txt
    ...
```
Load paths for your test and control groups
eg: output_location/sample_name/sample_name
```
test_group <- paste0(output_location, 'data/', sample_name, '/', sample_name)
control_group <- paste0(output_location, 'data/', sample_name, '/', sample_name)
```

Download the proper GTF and FASTA files from GENCODE and load them into R
```r
c_trans <- SpliceImpactR::get_c_trans("gencode.v45.pc_translations.fa")

c_nucs <- SpliceImpactR::get_c_nucs("gencode.v45.pc_transcripts.fa")

newGTF <- SpliceImpactR::setup_gtf('gencode.v45.annotation.gtf')
gtf <- newGTF$gtf
```

Run the entire pipeline for all AS events
```r
functionOutcomes <- SpliceImpactR::fullASoutcome(output_directory = "/output_location/",
              data_directory = "/data_directory/",
              test_group = test_group,control_group = control_group)
```

## Contributing
Contributions to SpliceImpactR are welcome, including bug reports, feature requests, and pull requests. Please see CONTRIBUTING.md for guidelines on how to contribute.

## Support
If you encounter any problems or have suggestions, please file an issue on the GitHub issue tracker. Or contact zachpw@bu.edu

## License
SpliceImpactR is available under the 

##Citation
If you use SpliceImpactR in your research, please cite:

```bibtex
Zachary Wakefield
SpliceImpactR: Analyzing the Impact of Alternative Splicing on Protein Structure and Function
2024
https://github.com/zachpwakefield/SpliceImpactR
```
