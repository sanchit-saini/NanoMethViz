# NanoMethViz

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/Shians/NanoMethViz/branch/master/graph/badge.svg)](https://codecov.io/gh/Shians/NanoMethViz?branch=master)
[![R-CMD-check](https://github.com/Shians/NanoMethViz/workflows/R-CMD-check/badge.svg)](https://github.com/Shians/NanoMethViz/actions)
<!-- badges: end -->

NanoMethViz is a toolkit for visualising methylation data from Oxford Nanopore sequencing.

## Installation

You can install NanoMethViz from Bioconductor with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("NanoMethViz")
```

To install the latest developmental version, use:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version='devel')

BiocManager::install("NanoMethViz")
```

## Usage

This package currently works with data from megalodon, nanopolish and f5c, to import your data please see the following vignette

``` r
vignette("ImportingData", package = "NanoMethViz")
```

An introductory example for plotting can be found in the package vignette:

``` r
vignette("Introduction", package = "NanoMethViz")
```

Other vignettes are provided for various features:

``` r
# how to use dimensionality reduction plots
vignette("DimensionalityReduction", package = "NanoMethViz")

# how to import external annotations
vignette("ExonAnnotations", package = "NanoMethViz")
```

## Examples

### MDS Plot

The MDS plot is used to visualise differences in the methylation profiles of
multiple samples.

![](img/mds.png)

### Feature Aggregation

The feature aggregation plot can average the methylation profiles across a set
of features.

![](img/agg_genes.png)

### Region Methylation Plot

The region methylation plot can visualise the methylation profile of a region
of interest. As well as provide a heatmap of the methylation along individual
reads.

![](img/peg3_gene.png)


## License

This project is licensed under Apache License, Version 2.0.
