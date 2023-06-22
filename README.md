
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bifurcatingr

<!-- badges: start -->
<!-- badges: end -->

Bifurcating autoregressive (BAR) models are commonly used to model
binary tree-structured data that appear in many applications, most
famously cell-lineage applications. The BAR model is an extension of the
autoregressive (AR) model where each line of descent considers an AR
process with the modification that the observations on the two sibling
cells who share the same parent are correlated. In practice, the BAR
model is used to explain the progression of single-cell proliferation.
The goal of the `bifurcatingr` package is to provide a collection of
functions that can be used for analyzing bifurcating autoregressive
data. The package implements the least squares estimation of bifurcating
autoregressive models of any order, p, BAR(p), and allows for executing
several types of bias correction on the least-squares estimators of the
autoregressive parameters including different types of confidence
intervals. Currently, the bias correction methods supported include
bootstrap (single, double, and fast-double) bias correction and
linear-bias-function -based bias correction. The library also contains
functions for generating and plotting bifurcating autoregressive data
from any BAR(p) model.

## Installation

You can install the development version of bifurcatingr like so:

``` r
You can install the released version of `bifurcatingr` from
[CRAN](https://CRAN.R-project.org) with:

    install.packages("bifurcatingr")
```

## Example

This is a basic example which shows you how to use `bifurcatingr` to fit
a bifurcating autoregressive model to the `ecoli` dataset which contains
the lifetimes of lineage E. coli cells.

Loading the `bifurcatingr` library and the `ecoli` dataset:

``` r
library(bifurcatingr)
#> Loading required package: fMultivar
data("ecoli")
## Fitting a BAR(p=1) model to the `ecoli` dataset:

    bfa.ls(ecoli$lifetime, p = 1, conf = TRUE, conf.level = 0.95, 
           p.value = TRUE, cov.matrix = TRUE)
#> $coef
#>      Intercept  X_[t/2]
#> [1,]  17.61658 0.355198
#> 
#> $p.value
#>       Intercept  X_[t/2]
#> [1,] 0.01079535 0.181183
#> 
#> $error.cor
#> [1] 0.5836975
#> 
#> $cov.matrix
#>            Intercept    X_[t/2]
#> Intercept 1480.39874 -56.147524
#> X_[t/2]    -56.14752   2.187566
#> 
#> $asymptotic.ci
#>                 2.5%      97.5%
#> Intercept  4.0722831 31.1608852
#> X_[t/2]   -0.1654543  0.8758503
#> 
#> $bootstarp.ci
#>                2.5%     97.5%
#> Intercept  2.787469 32.445699
#> X_[t/2]   -0.162691  0.873087
#> 
#> $percentile.ci
#>                  2.5%      97.5%
#> Intercept  3.91627094 29.0253668
#> X_[t/2]   -0.05318311  0.8051654
```
