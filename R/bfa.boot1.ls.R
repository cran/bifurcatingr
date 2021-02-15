#' Single Bootstrap of Least Squares Estimators of BAR(p) Models
#'
#' This function performs single bootstrapping of the least squares estimators
#' of the autoregressive coefficients in a bifurcating autoregressive (BAR)
#' model of any order \code{p} as described in Elbayoumi & Mostafa (2020).
#'
#' @param z a numeric vector containing the tree data
#' @param p an integer determining the order of bifurcating autoregressive model
#'   to be fit to the data
#' @param burn number of tree generations to discard before starting the
#'   bootstrap sample (replicate)
#' @param B number of bootstrap samples (replicates)
#' @param boot.est a logical that determines whether the bootstrapped least
#'   squares estimates of the autoregressive coefficients should be returned.
#'   Defaults to TRUE.
#' @param boot.data a logical that determines whether the bootstrap samples
#'   should be returned. Defaults to FALSE.
#' @return \item{boot.est}{a matrix containing the bootstrapped least squares
#'   estimates of the autoregressive coefficients} \item{boot.data}{a matrix
#'   containing the bootstrap samples used}
#' @references Elbayoumi, T. M. & Mostafa, S. A. (2020). On the estimation bias
#'   in bifurcating autoregressive models. \emph{Stat}, 1-16.
#' @export
#' @examples
#' z <- bfa.tree.gen(31, 1, 1, 1, 0.5, 0.5, 0, 10, c(0.7))
#' bfa.boot1.ls(z, p=1, B=999)
bfa.boot1.ls <- function(z, p, burn = 5, B, boot.est=TRUE, boot.data=FALSE){  # burn can be determined by simulation
  n=length(z)
  boot.dat <- matrix(rep(NA, n*B), ncol=n)
  est <- bfa.ls(z,p,resids=TRUE)
  intercept <- est$coef[1,1]
  coef <- t(as.vector(est$coef[1,2:(p+1)]))
  e <- est$resids-mean(est$resids,na.rm=TRUE)
  e.e=matrix(e,ncol=2, byrow=T)
  ee <- e.e[stats::complete.cases(e.e),]
  for (h in 1:B){
    # start-up tree
    n0 = 2^(p + burn) - 1
    term.matrix0 <- matrix(rep(NA,p*n0),ncol=p)
    dat0 <- rep(NA,n0)
    er0 <- matrix(rep(NA,((n0-1)/2)),ncol=2,nrow=((n0-1)/2))
    samp0 <- sample(1:nrow(ee), ((n0-1)/2) ,replace=T)
    for (s in 1: nrow(er0)) {
      l <- samp0[s]
      er0[s,] <- ee[l,]
    }
    error0 <- rep(NA,n0)
    for(i in 1:((n0-1)/2)){
      error0[(i*2)] <- er0[i,1]
      error0[(i*2)+1] <- er0[i,2]
    }
    dat0[1:((2^p)-1)] <- z[1:((2^p)-1)]
    for(j in (2^p) : n0){
      for(k in 1:p){
        term.matrix0[,k] <- coef[1,k]*dat0[floor(j/(2^k))]
      }
      dat0[j] <- intercept+sum(term.matrix0[j,])+error0[j]
    }
  # bootstrap tree
    term.matrix <- matrix(rep(NA,p*n),ncol=p)
    dat <- rep(NA,n)
    er <- matrix(rep(NA,((n-1)/2)),ncol=2,nrow=((n-1)/2))
    samp <- sample(1:nrow(ee), ((n-1)/2) ,replace=T)
    for (s in 1: nrow(er)) {
      l <- samp[s]
      er[s,] <- ee[l,]
    }
    error <- rep(NA,n)
    for(i in 1:((n-1)/2)){
      error[(i*2)] <- er[i,1]
      error[(i*2)+1] <- er[i,2]
    }
    dat[1:((2^p)-1)] <- dat0[bfa.subtree(n0,p)]
    for(j in (2^p) : n){
      for(k in 1:p){
        term.matrix[,k] <- coef[1,k]*dat[floor(j/(2^k))]
      }
      dat[j] <- intercept+sum(term.matrix[j,])+error[j]
    }
  # copying missing values as in the original tree
     for (k in 1:n){
      dat[k] <- ifelse(is.na(z[k])==T, NA, dat[k])
     }
    boot.dat[h,] <- dat
  }
  out <- list()
  if(boot.est)
    out$boot.est=matrix(unlist(apply(boot.dat, 1, bfa.ls, p=p, error.cor=FALSE)),nrow=B, byrow=TRUE)
  if(boot.data)
    out$boot.data <- boot.dat
  return(out)
}
