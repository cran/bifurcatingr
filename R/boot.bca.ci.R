#' Bias-Corrected and Accelerated bootstrap Confidence Interval (BCa) Function.
#'
#' This function calculates the Bias-Corrected and Accelerated bootstrap (BCa)
#' CI for the least squares estimator of the bifurcating autoregressive model.
#'
#' @param z a numeric vector containing the tree data.
#' @param p an integer determining the order of bifurcating autoregressive model to be fit to the data.
#' @param conf.level A numeric value representing the confidence level. Defaults to 0.95.
#' @param B number of bootstrap samples (replicates). Defaults to 99.
#' @param J an integer determining the number of vectors of bootstrap replicates. Defaults to 2.
#' @return \item{coef}{a matrix containing the least squares estimates of the autoregressive
#' coefficients}\item{limits}{A numeric vector representing the lower and upper limits of the BCa
#' confidence interval for the autoregressive coefficients of BAR model}
#' @importFrom utils tail
#' @importFrom stats pnorm qnorm quantile
#' @export
#' @examples
#' # Generating Non-contaminated normal BAR(1) tree and calculating the BCa CI for
#' # the autoregressive coefficients of the BAR(1) model
#' z <- bfa.tree.gen(31, 1, 1, 1, 0.5, 0.5, 0, 10, c(0.7))
#' boot.bca.ci(z,p=1,B=99,J=2,conf.level=0.95)
#' # Generating Non-contaminated normal BAR(2) tree and calculating the BCa CI for
#' # the autoregressive coefficients of the BAR(2) model
#' z <- bfa.tree.gen(127, 2, 1, 1, 0.5, 0.5, 0, 10, c(0.5, 0.3))
#' boot.bca.ci(z,p=2,B=99,J=2,conf.level=0.95)
boot.bca.ci <- function(z,p,B=99,J=2,conf.level=0.95){
  n=length(z)
  nn=floor((n-1)/J)
  Rs=c(J:((2*J)-1))
  M=length(Rs)
  coef.ls=bfa.ls(z,p)$coef
  dat.matrx = matrix(NA, nrow=J, ncol=nn)
  limits = matrix(NA, nrow=p, ncol=2)
  dat=list()
  for (s in 1:p) {
    # Bootstarp LS estimates
    phis <- matrix(unlist(apply(matrix(unlist(bfa.boot1.ls(z,p,B=B, boot.est=FALSE,boot.data = TRUE)), ncol=n, nrow = B),1,p=p, bfa.ls)),ncol=p+2,byrow=TRUE)[,s+1]
    alpha <- (1 + c(-conf.level, conf.level))/2
    zalpha <- qnorm(alpha)
    # the bias correction factor
    z0 <- qnorm(sum(as.vector(phis) < coef.ls[1,s+1]) / B)
    # the acceleration factor (jackknife est.)
    for(R in Rs){
      indx=NULL
      for (i in 0:10){
        ind=((2^i)*R):((2^i)*R+(2^i-1))
        indx=c(indx,if(tail(ind, n=1)<=n){ind})
      }
      dat=append(dat, list(indx))
    }
    for (f in 1:M){
      dat.matrx[f,] = z[dat[[f]]]
    }
    phi.jack <- matrix(unlist(apply(dat.matrx,1,p=p, bfa.ls)),ncol=p+2,byrow=TRUE)[,s+1]
    L <- mean(phi.jack) - phi.jack
    a <- sum(L^3)/(6 * sum(L^2)^1.5)
    # BCa conf. limits
    adj.alpha <- pnorm(z0 + (z0+zalpha)/(1-a*(z0+zalpha)))
    limits[s,1] <- quantile(phis, adj.alpha, names=F, type=8,na.rm = TRUE)[1]
    limits[s,2] <- quantile(phis, adj.alpha, names=F, type=8,na.rm = TRUE)[2]
  }
  return(list("est"=coef.ls, "BCa"=limits))

}
