#' Percentile Bootstrap Confidence Interval Function
#'
#' This function calculates the percentile  bootstrap CI for the least squares estimator of the
#' bifurcating autoregressive model.
#' @param a1_ls_star A numeric vector representing B replicates of the least squares estimator.
#' @param conf_level A numeric value representing the confidence level. Defaults to 0.95.
#' @importFrom stats quantile rnorm
#' @return A numeric vector representing the lower and upper limits of the confidence interval.
#'
#' @export
#' @examples
#' a1_ls_star <- c(rnorm(100,0.7,0.05))
#' bfa_perc_ci(a1_ls_star, conf_level= 0.95)
bfa_perc_ci <- function(a1_ls_star, conf_level=0.95){
  alpha <- 1 - conf_level
  perc_ci = quantile(a1_ls_star, prob=c(alpha/2, 1-alpha/2))
}
