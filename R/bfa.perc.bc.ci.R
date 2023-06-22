#' Bias Correction Percentile Confidence Interval Function
#'
#' This function calculates the Bias-Corrected percentile CI for the least
#' squares estimator of the bifurcating autoregressive model.
#' @param z a numeric vector containing the tree data.
#' @param a1.ls A numeric value of the least squares estimator of bifurcating autoregressive model
#' @param a1.ls.star A numeric vector representing B replicates of the least squares estimator
#' @param conf.level A numeric value representing the confidence level. Defaults to 0.95.
#' @return A numeric vector representing the lower and upper limits of the bias corrected percentile
#' confidence interval for the autoregressive coefficients of BAR model.
#' @importFrom stats qnorm pnorm rnorm
#' @export
#' @examples
#' # Generating Non-contaminated normal BAR(1) tree and calculating the bias
#' # corrected percentile CI for the autoregressive coefficients of the BAR(1) model
#' z <- bfa.tree.gen(31, 1, 1, 1, 0.5, 0.5, 0, 10, c(0.7))
#' a1.ls <- bfa.ls(z,1)$coef[2]
#' a1.ls.star <- (rnorm(100,0.7,0.05))
#' bfa.perc.bc.ci(z, a1.ls, a1.ls.star, conf.level=0.95)
bfa.perc.bc.ci <- function(z, a1.ls, a1.ls.star, conf.level=0.95) {
  # generate the bootstrap samples
  n <- length(z)
  B <- length(a1.ls.star)
  # compute z.hat
  z.hat <- qnorm(sum(a1.ls.star < a1.ls)/B)
  # pick out the BC quantiles
  alpha = 1 - conf.level
  p <- c(alpha/2, 1 - alpha/2)
  w <- 2*z.hat + qnorm(p)
  q <- pnorm(w)
  ind <- ceiling(q*B)
  a1.ls.star <- sort(a1.ls.star)
  return(percbc.ci = a1.ls.star[ind])
}
