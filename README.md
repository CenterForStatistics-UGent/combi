
Manual for the use of the compIntegrate package
===============================================

Install and load packages
-------------------------

This repo contains R-code to fit and plot the mode-based integration models for compositional omics data. The basic usage is demonstrated here.

The package can be installed from BioConductor and loaded using the following commands:

``` r
library(BiocManager)
install("compIntegrate", update = FALSE)
```

``` r
suppressPackageStartupMessages(library(compIntegrate))
cat("compIntegrate package version", as.character(packageVersion("compIntegrate")), "\n")
```

    ## compIntegrate package version 0.1.0

Alternatively, the latest version can be installed directly from this GitHub repo as follows:

``` r
library(devtools)
install_github("CenterForStatistics-UGent/compIntegrate")
```
