#' Scatterplots for Bifurcating Autoregressive Data
#'
#' Draw scatterplots between observations at time \code{t} and the lagged
#' observations from the given bifurcating autoregressive tree data.
#' @param p an integer determining the order of the bifurcating autoregressive
#'   model that is believed to best fit the data
#' @inheritParams bfa.ls
#' @param ... other graphical parameters that can be passed to \code{plot()} or
#'   \code{pairs()} (see \code{\link[graphics]{par}} and
#'   \code{\link[graphics]{pairs}})
#' @return A single scatterplot when \code{p=1} or a matrix of scatterplots when
#'   \code{p>1}.
#' @export
#' @examples
#' z <- bfa.tree.gen(63, 1, 1, 2, 0.5, 0.5, 0.2, 10, c(0.7))
#' bfa.scatterplot(z,1)
#' z<-bfa.tree.gen(63, 2, 1, 2, 0.5, 0.5, 0.2, 10, c(0.7,0.2))
#' bfa.scatterplot(z,2)
#' bfa.scatterplot(z,2,lower.panel=NULL)

bfa.scatterplot <- function(z, p, ...){
  n <- length(z)
  m <- (n-1)/2
  f <- c(1:p)
  y <- z[c((2^p):n)]
  x <- matrix(rep(NA,p*length(y)),ncol=p)
  d <- as.numeric(1:p)
  for (i in 1:p){
    d[i] <- floor(n/(2^i))}
  g <- rev(c(2^(1:p-1)))
  for (k in 1:p){
    x[,k] <- rep(z[g[k]:d[k]],each=2^f[k])}
  dat <- cbind(matrix(x, ncol=p),y)
  colnames(dat) <- c( paste0('X_[t/',2^c(1:p),"]"),'X_t')
  if(p==1) plot(dat, ...)
  if(p>1) graphics::pairs(dat, labels=colnames(dat), ...)
}
