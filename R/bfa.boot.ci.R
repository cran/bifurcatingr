#' Standard Normal Bootstrap Confidence Interval Function
#'
#' This function calculates the standard normal bootstrap CI for the least squares estimator of the
#' bifurcating autoregressive model.
#'
#' @param a1.ls A numeric value of the least squares estimator of bifurcating autoregressive model
#' @param a1.ls.star A numeric vector representing B replicates of the least squares estimator
#' @param conf.level A numeric value representing the confidence level. Defaults to 0.95.
#' @return A numeric vector representing the lower and upper limits of the confidence interval
#' @importFrom stats qnorm sd rnorm
#' @export
#' @examples
#' a1.ls <- 0.7
#' a1.ls.star <- c(rnorm(100,0.7,0.05))
#' bfa.boot.ci(a1.ls, a1.ls.star, conf.level= 0.95)
bfa.boot.ci <- function(a1.ls, a1.ls.star, conf.level = 0.95){
  alpha = 1 - conf.level
  se = sd(a1.ls.star)
  a1.L = a1.ls - qnorm(1-alpha/2)*se
  a1.U = a1.ls + qnorm(1-alpha/2)*se
  return(boot.ci = c(a1.L,a1.U))
}
