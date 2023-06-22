#' Percentile Bootstrap Confidence Interval Function
#'
#' This function calculates the percentile  bootstrap CI for the least squares estimator of the
#' bifurcating autoregressive model.
#' @param a1.ls.star A numeric vector representing B replicates of the least squares estimator.
#' @param conf.level A numeric value representing the confidence level. Defaults to 0.95.
#' @importFrom stats quantile rnorm
#' @return A numeric vector representing the lower and upper limits of the confidence interval.
#'
#' @export
#' @examples
#' a1.ls.star <- c(rnorm(100,0.7,0.05))
#' bfa.perc.ci(a1.ls.star, conf.level= 0.95)
bfa.perc.ci <- function(a1.ls.star, conf.level=0.95){
  alpha = 1 - conf.level
  return(perc.ci = quantile(a1.ls.star, prob=c(alpha/2, 1-alpha/2)))
}
