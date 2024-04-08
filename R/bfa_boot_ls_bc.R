#' Bootstrap of Bias-Correction Least Squares Estimators of BAR(p) Models
#'
#' This function performs linear-bias-function bias-correction (LBC), single
#' bootstrap, double bootstrap, fast-double bootstrap of the bias-correction
#' least squares estimators of the autoregressive coefficients in a bifurcating
#' autoregressive (BAR) model of any order \code{p} as described in Elbayoumi &
#' Mostafa (2020).
#'
#' @param z a numeric vector containing the tree data
#' @param p an integer determining the order of bifurcating autoregressive model
#'   to be fit to the data
#' @param method method of bias correction. Currently, "boot1", "boot2",
#'   "boot2fast" and "LBC" are supported and they implement single bootstrap,
#'   double bootstrap, fast-double bootstrap, and linear-bias-function
#'   bias-correction, respectively.
#' @param burn number of tree generations to discard before starting the
#'   bootstrap sample (replicate)
#' @param B number of bootstrap samples (replicates)
#' @param boot_est a logical that determines whether the bootstrapped least
#'   squares estimates of the autoregressive coefficients should be returned.
#'   Defaults to TRUE.
#' @param boot_data a logical that determines whether the bootstrap samples
#'   should be returned. Defaults to FALSE.
#' @return \item{boot_bcest}{a matrix containing the bootstrapped bias-correction
#'   least squares estimates of the autoregressive coefficients} \item{boot_data}{a matrix
#'   containing the bootstrap samples used}
#' @references Elbayoumi, T. M. & Mostafa, S. A. (2020). On the estimation bias
#'   in bifurcating autoregressive models. \emph{Stat}, 1-16.
#' @export
#' @examples
#' z <- bfa_tree_gen(31, 1, 1, 1, 0.5, 0.5, 0, 10, c(0.7))
#' bfa_boot_ls_bc(z, p=1, method="LBC", B=500)
#' hist(bfa_boot_ls_bc(z, p=1, method="LBC", B=500)$boot_bcest)
bfa_boot_ls_bc <- function(z, p, method="boot1",burn = 5, B, boot_est=TRUE, boot_data=FALSE){
  n <- length(z)
  boot_dat <- matrix(rep(NA, n*B), ncol=n)
  est <- bfa_ls(z,p,resids=TRUE)
  intercept <- est$coef[1,1]
  coef <- t(as.vector(est$coef[1,2:(p+1)]))
  e <- est$resids-mean(est$resids,na.rm=TRUE)
  e_e <- matrix(e,ncol=2, byrow=T)
  ee <- e_e[stats::complete.cases(e_e),]
  for (h in 1:B){
    # start-up tree
    n0 <- 2^(p + burn) - 1
    term_matrix0 <- matrix(rep(NA,p*n0),ncol=p)
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
        term_matrix0[,k] <- coef[1,k]*dat0[floor(j/(2^k))]
      }
      dat0[j] <- intercept+sum(term_matrix0[j,])+error0[j]
    }
    # bootstrap tree
    term_matrix <- matrix(rep(NA,p*n),ncol=p)
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
    dat[1:((2^p)-1)] <- dat0[bfa_subtree(n0,p)]
    for(j in (2^p) : n){
      for(k in 1:p){
        term_matrix[,k] <- coef[1,k]*dat[floor(j/(2^k))]
      }
      dat[j] <- intercept+sum(term_matrix[j,])+error[j]
    }
    # copying missing values as in the original tree
    for (k in 1:n){
      dat[k] <- ifelse(is.na(z[k])==T, NA, dat[k])
    }
    boot_dat[h,] <- dat
  }
  out <- list()
  if(boot_est)
    out$boot_bcest=matrix(unlist(apply(boot_dat, 1, bfa_ls_bc, p=p, method=method)),nrow=B, byrow=TRUE)
  if(boot_data)
    out$boot_data <- boot_dat
  out
}
