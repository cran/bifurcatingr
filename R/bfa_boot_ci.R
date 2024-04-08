#' Standard Normal Bootstrap Confidence Interval Function
#'
#' This function calculates the standard normal bootstrap CI for the least squares estimator of the
#' bifurcating autoregressive model.
#'
#' @param a1_ls A numeric value of the least squares estimator of bifurcating autoregressive model
#' @param a1_ls_star A numeric vector representing B replicates of the least squares estimator
#' @param conf_level A numeric value representing the confidence level. Defaults to 0.95.
#' @return A numeric vector representing the lower and upper limits of the confidence interval
#' @importFrom stats qnorm sd rnorm
#' @export
#' @examples
#' a1_ls <- 0.7
#' a1_ls_star <- c(rnorm(100,0.7,0.05))
#' bfa_boot_ci(a1_ls, a1_ls_star, conf_level= 0.95)
bfa_boot_ci <- function(a1_ls, a1_ls_star, conf_level = 0.95){
  alpha <- 1 - conf_level
  se <- sd(a1_ls_star)
  a1_L <- a1_ls - qnorm(1-alpha/2)*se
  a1_U <- a1_ls + qnorm(1-alpha/2)*se
  boot_ci <- c(a1_L,a1_U)
}
