
<!-- README.md is generated from README.Rmd. Please edit that file -->

# proteomicsHF

<!-- badges: start -->

[![GitHub
license](https://img.shields.io/github/license/Naereen/StrapDown.js.svg)](https://github.com/pranavdorbala/proteomicsHF/blob/main/LICENSE)
[![ForTheBadge
built-by-developers](http://ForTheBadge.com/images/badges/built-by-developers.svg)](https://github.com/pranavdorbala/proteomicsHF)
<!-- badges: end -->

proteomicsHF is a package with Functions necessary to analyze and
replicate results of biomarker discovery pipeline for SOMALogic 5K Assay
proteins associated with incident heart failure. Results are cached and
derived from the Atherosclerosis Risk in Communities (ARIC) cohort from
visits three and five. Includes functions to set adjustment variables,
impute missing variables, summarize data, conduct uni-protein and
multi-protein analyses and visualize data. Also includes a number of
utility functions for data manipulation and reporting.

## Installation

Install the latest version of proteomicsHF from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pranavdorbala/proteomicsHF")
```

## Set Up

This is is how you set up multicore processing:

``` r
library(proteomicsHF)
#> Loading required package: magrittr
#> Loading required package: rlang
#> 
#> Attaching package: 'rlang'
#> The following object is masked from 'package:magrittr':
#> 
#>     set_names
#> 
#> Attaching package: 'proteomicsHF'
#> The following object is masked from 'package:stats':
#> 
#>     rf

set.seed(101010)

#Set number of cores according to system
doParallel::registerDoParallel(8) # 2 * system cores
future::plan(future::multisession, workers = 4) # number of system cores
```
