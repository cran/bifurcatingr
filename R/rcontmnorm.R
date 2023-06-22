#' Contaminated Normal Generator
#'
#' This function generates contaminated multivariate normal errors that are used
#' in the generation of the Bifurcating autoregressive tree.
#' @param n sample size
#' @param d dimension. Defaults to 2 for bivariate normal errors.
#' @param mu1 mean vector for first multivariate normal distribution. Defaults
#'   to the zero vector.
#' @param mu2 mean vector for second multivariate normal distribution. Defaults
#'   to the zero vector.
#' @param sigma1 variance-covariance matrix for first multivariate normal
#'   distribution. Defaults to the d by d identity matrix, where d is the
#'   dimension.
#' @param sigma2 variance-covariance matrix for second multivariate normal
#'   distribution. Defaults to the d by d identity matrix, where d is the
#'   dimension.
#' @param g proportion of contamination. Defaults to zero producing
#'   non-contaminated multivariate normal data.
#' @return An \code{n} by \code{d} contaminated multivariate normal matrix.
#' @export
#' @examples
#' #Non-contaminated bivariate normal:
#' rcontmnorm(10, sigma2=2^2*matrix(c(1,0,0,1),nrow=2) , g=0)
#' #Contaminated bivariate normal with 20% contamination:
#' rcontmnorm(10, sigma2=2^2*matrix(c(1,0,0,1),nrow=2) , g=0.2)
rcontmnorm <- function(n, d=2, mu1=rep(0,d), sigma1=diag(d), mu2=rep(0,d), sigma2=diag(d), g){
  bin <- as.numeric(stats::runif(n)<g)
  (1-bin)*rmnorm(n, d, mu1, sigma1) + bin*rmnorm(n, d, mu2, sigma2)
}
