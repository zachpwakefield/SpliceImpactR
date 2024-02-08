# SpliceImpactR
SpliceImpactR is an R package designed for bioinformatics researchers interested in studying the impact of alternative splicing on protein structure and function. It provides tools for analyzing RNA-seq data to identify splicing events and predict their consequences on the resulting protein products.

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

## Usage
```r
library(SpliceImpactR)

output_location <- "/projectnb/evolution/zwakefield/proteinImpacts/feb6_7/"
system(paste0("mkdir ",  output_location))
pdir <- system.file(package="domainEnrichment")

sample_sheet <- read_csv('/projectnb2/evolution/zwakefield/tcga/metadata/tcga_data/complete_sample_sheet.csv')
tumor <- sample_sheet$File.ID[sample_sheet$Project.ID == "TCGA-BRCA" & sample_sheet$Sample.Type == "Primary Tumor"][1:5]
control <- sample_sheet$File.ID[sample_sheet$Project.ID == "TCGA-BRCA" & sample_sheet$Sample.Type == "Solid Tissue Normal"][1:5]

organizeSamples(paste0('/projectnb2/evolution/zwakefield/tcga/runs/', c(tumor, control), '/'), cores = 2)

test_group <- paste0(output_location, 'data/', tumor, '/', tumor)
control_group <- paste0(output_location, 'data/', control, '/', control)

diSE <- domainEnrichment::differential_inclusion_rMATS(test_names = test_group, control_names = control_group, et = "SE", cores = 8, outlier_threshold = "4/n", min_proportion_samples_per_phenotype = .333)
diHIT <- domainEnrichment::differential_inclusion(test_names = test_group, control_names = control_group, cores = 16, outlier_threshold = "4/n")

c_trans <- domainEnrichment::get_c_trans(pdir)
gtf <- domainEnrichment::get_gtf(pdir)

bg_input <- gsub("[^/]*$", "", c(control_group, test_group))
bg <- getBackground(input=bg_input,
                    mOverlap = .2,
                    cores = 30,
                    nC = length(control_group),
                    nE = length(test_group),
                    exon_type = "SE", 
                    pdir = pdir,
                    output_location = output_location)

fg <- getForeground(input=diSE,
                    test_names = test_group, 
                    control_names = control_group,
                    thresh = .2,
                    fdr=.05,
                    mOverlap=.2,
                    cores=4,
                    nC=length(control_group),
                    nE=length(test_group),
                    exon_type="SE",
                    pdir=pdir,
                    output_location=output_location)

#####
## Extract pairs with +/-

pfg <- getPaired(foreground = fg$proBed)

#####
pfam <- domainEnrichment::getPfam(foreground = fg, background = bg, pdir = pdir, cores = 4, output_location = output_location)


gD <- getData(output_location = output_location, fdr_use = .05, min_sample_success = 3, engine = "Pfam")

iDDI <- domainEnrichment::init_ddi(pdir = pdir, output_location = output_location, ppidm_class = "Gold", removeDups = T)

# ncol.g <- read.table(paste0(output_location,  " tti_igraph_edgelist_ Gold _removeDups"), sep = " ", row.names = NULL)
tti <- getTTI(paired_foreground = pfg$paired_proBed,
       pdir = pdir,
       steps=1,
       max_vertices_for_viz = 5000,
       fdr = .05,
       plot_bool = T,
       ppidm_class = "Gold",
       write_igraphs = T,
       ddi = "Gold",
       ddi_type = "pdm",
       output_location = output_location)
```

## Contributing
Contributions to SpliceImpactR are welcome, including bug reports, feature requests, and pull requests. Please see CONTRIBUTING.md for guidelines on how to contribute.

## Support
If you encounter any problems or have suggestions, please file an issue on the GitHub issue tracker.

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
