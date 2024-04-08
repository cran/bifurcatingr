#' Multivariate Normal Generator
#'
#' This function generates multivariate normal errors that are used in the
#' generation of the Bifurcating autoregressive tree.
#' @param n sample size
#' @param d dimension. Defaults to 2 for bivariate normal errors.
#' @param mu mean vector. Defaults to the zero vector.
#' @param sigma variance-covariance matrix. Defaults to the \code{d} by \code{d}
#'   identity matrix, where d is the dimension.
#' @return An \code{n} by \code{d} multivariate normal matrix.
#' @export
#' @examples
#' rmnorm(10)
#' rmnorm(10, 3)
rmnorm <- function(n, d=2, mu=rep(0,d), sigma=diag(d)){
  temp <- eigen(sigma)
  lambda <- diag(temp$values)
  A <- temp$vectors%*%(lambda^0.5)%*%t(temp$vectors)
  z <- matrix(stats::rnorm(n*d,0,1),ncol=n,nrow=d)
  x <- A%*%z + mu
  t(x)
}
