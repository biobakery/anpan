
<!-- README.md is generated from README.Rmd. Please edit that file -->

# anpan

<!-- badges: start -->
<!-- badges: end -->

The goal of anpan is to consolidate statistical methods for strain
analysis. This includes testing genetic elements for association with
outcomes, phylogenetic association testing, and ….

## Requirements

anpan requires the following R packages, most of which which can be
installed from CRAN:

    install.packages(c("ape", 
                       "brms", 
                       "data.table", 
                       "dplyr", 
                       "purrr", 
                       "tibble", 
                       "tidyr", 
                       "furrr", 
                       "ggnewscale",
                       "patchwork", 
                       "progressr",
                       "rstanarm")) # Ncpus = 4 to go faster

    install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

If the `cmdstanr` installation doesn’t work you can find more detailed
instructions [at this link](https://mc-stan.org/cmdstanr/).

## Installation

You can install anpan from github with:

``` r
devtools::install("biobakery/anpan")
```

## Example - element testing

You can filter large uniref tables and look for associations with
outcomes (while controlling for covariates) with `anpan_batch`. Right
now this only works for glm-based models.

``` r
library(tidyverse)
library(data.table)
library(furrr); plan(multisession, workers = 2) # Run with two background processes on the local machine
library(progressr); handlers(global=TRUE) # show progress bars for long computations
library(anpan)


anpan_batch(bug_dir = "/path/to/gene_families/",
            meta_file = "/path/to/metadata.tsv",
            out_dir = "/path/to/output",
            annotation_file = "/path/to/annotation.tsv", #optional, used for plots
            filtering_method = "kmeans",
            model_type = "glm",
            covariates = c("age", "gender"),
            plot_ext = "pdf",
            save_filter_stats = TRUE)
```

This will create the output directory with tables of regression
coefficients, top-results plots, and filter statistics.

## Example - phylogenetic generalized linear mixed model

You can run a PGLMM with the `anpan_pglmm` function. This fits a model
where the covariance matrix of outcomes is determined from the
phylogenetic tree. Right now there’s only gaussian outcomes, and no
regularization on the covariance matrix.

``` r
meta_file = "/path/to/metadata.tsv" 
tree_file = "/path/to/file.tre"

anpan_pglmm(meta_file,
            tree_file,
            trim_pattern = "_bowtie2",
            covariates = NULL,
            plot_cov_mat = TRUE,
            outcome = "age",
            verbose = TRUE,
            cores = 4)
```
