#' Bias Correction Percentile Confidence Interval Function
#'
#' This function calculates the Bias-Corrected percentile CI for the least
#' squares estimator of the bifurcating autoregressive model.
#' @param z a numeric vector containing the tree data.
#' @param a1_ls A numeric value of the least squares estimator of bifurcating autoregressive model
#' @param a1_ls_star A numeric vector representing B replicates of the least squares estimator
#' @param conf_level A numeric value representing the confidence level. Defaults to 0.95.
#' @return A numeric vector representing the lower and upper limits of the bias corrected percentile
#' confidence interval for the autoregressive coefficients of BAR model.
#' @importFrom stats qnorm pnorm rnorm
#' @export
#' @examples
#' # Generating Non-contaminated normal BAR(1) tree and calculating the bias
#' # corrected percentile CI for the autoregressive coefficients of the BAR(1) model
#' z <- bfa_tree_gen(31, 1, 1, 1, 0.5, 0.5, 0, 10, c(0.7))
#' a1_ls <- bfa_ls(z,1)$coef[2]
#' a1_ls_star <- (rnorm(100,0.7,0.05))
#' bfa_perc_bc_ci(z, a1_ls, a1_ls_star, conf_level=0.95)
bfa_perc_bc_ci <- function(z, a1_ls, a1_ls_star, conf_level=0.95) {
  # generate the bootstrap samples
  n <- length(z)
  B <- length(a1_ls_star)
  # compute z_hat
  z_hat <- qnorm(sum(a1_ls_star < a1_ls)/B)
  # pick out the BC quantiles
  alpha = 1 - conf_level
  p <- c(alpha/2, 1 - alpha/2)
  w <- 2*z_hat + qnorm(p)
  q <- pnorm(w)
  ind <- ceiling(q*B)
  a1_ls_star <- sort(a1_ls_star)
  percbc_ci <- a1_ls_star[ind]
}
