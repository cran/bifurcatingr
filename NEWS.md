# bifurcatingr 2.1.0
The function naming convention was separated with dots  (`.`). It is modified to use underscore instead ( '_'). For example, the function "bfa.ls" becomes "bfa_ls".  in addition, vignette file is developed and published with this version.


# bifurcatingr 2.0.0

## NEW FEATURES
Several corrected and uncorrected confidence intervals for the lease squares estimators of the bifurcating autoregressive model with any order p are added. The following confidence intervals are available in the new version (2.0.0):

 * Uncorrected standard normal bootstrap confidence interval.
 - Standard normal bootstrap confidence interval with single bootstrap bias correction.
 - Standard normal bootstrap confidence interval with fast double bootstrap bias correction.
 - Uncorrected percentile bootstrap confidence interval.
 - Percentile bootstrap confidence interval with single bootstrap bias correction.
 - Percentile bootstrap confidence interval with fast double bootstrap bias correction.
 - Bias-corrected and accelerated bootstrap confidence interval.


## Corrections
* bfa.scatterplot function that allows to create the scatterplots between observations at time t and the lagged observations from the given bifurcating autoregressive tree data fixed. Currently, it gives correct axes' labels.    
