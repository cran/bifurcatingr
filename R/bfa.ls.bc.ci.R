#' Bias-Corrected Confidence intervals of Least Squares Estimators for Bifurcating
#' Autoregressive Models
#'
#' This function performs bias correction confidence intervals on the least squares
#' estimators of the autoregressive coefficients in a BAR(p) model using single,
#' fast-double, and the Bias-corrected and accelerated bootstrapping as described
#' in Elbayoumi and Mostafa (2023).
#'
#' @param z a numeric vector containing the tree data
#' @param p an integer determining the order of bifurcating autoregressive model
#'   to be fit to the data
#' @param method method of bias correction. Currently, "boot1", "boot2fast", and
#'  "BCa" are supported and they implement single bootstrap, fast-double bootstrap,
#'   and bias-corrected and accelerated bootstrap, respectively. Defaults to "BCa".
#' @param conf.int type of the confidence interval. Currently, "standard.normal.bc",
#'  "percentile", and "percentile.bc" are supported and they implement corrected
#'  standard normal bootstrap CI, uncorrected percentile bootstrap CI, and corrected
#'  percentile bootstrap CI, respectively. If "boot1" method is selected, the
#'  "standard.normal.bc", "percentile", "percentile.bc" confidence intervals can
#'  be obtained. If "boot2fast" method is selected, the "standard.normal.bc" and
#'  "percentile.bc" confidence intervals can be obtained. No effect for conf.int,
#'  the "BCa" method is selected. Defaults to standard.normal.bc".
#' @param conf.level confidence level to be used in computing confidence intervals
#'   for model coefficients. Defaults to \code{0.95}.
#' @param B number of bootstrap samples (replicates).
#' @param burn number of tree generations to discard before starting the
#'  bootstrap sample (replicate). Defaults to 5.
#' @return \item{Bias.corrected.coef}{a matrix containing the
#'  bias-correction least squares estimates of the autoregressive coefficients}
#'  \item{BCa.ci}{a matrix containing the lower and upper limits of corrected BCa
#'  confidence intervals,if \code{method="BCa"}}
#'  \item{standard.normal.bc.ci}{a matrix containing the lower and upper limits of
#'  corrected confidence intervals, if \code{method="boot1"} and \code{conf.int="standard.normal.bc"}
#'  or \code{conf.int="All"}}
#'  \item{percentile.ci}{a matrix containing the lower and upper limits of uncorrected
#'  percentile confidence intervals, if \code{method="boot1"} and \code{conf.int="percentile"}
#'  or \code{conf.int="All"}}
#'  \item{percentile.bc.ci}{a matrix containing the lower and upper limits of corrected
#'  percentile confidence intervals, if \code{method="boot1"} and \code{conf.int="percentile.bc"}
#'  or \code{conf.int="All"}}
#'  \item{standard.normal.bc.ci}{a matrix containing the lower and upper limits of corrected
#'  confidence intervals, if \code{method="boot2fast"} and \code{conf.int="standard.normal.bc"}
#'  or \code{conf.int="All"}}
#'  \item{percentile.bc.ci}{a matrix containing the lower and upper limits of corrected percentile
#'  confidence intervals, if \code{method="boot2fast"} and \code{conf.int="percentile.bc"} or \code{conf.int="All"}}
#' @references Elbayoumi, T. M. & Mostafa, S. A. (2023). Impact of Bias Correction
#'  of the Least Squares Estimation on Bootstrap Confidence Intervals for Bifurcating
#'  Autoregressive Models. \emph{Journal of Data Science}, 1-20, doi.org/10.6339/23-JDS1092.
#' @export
#' @examples
#' # Generating Non-contaminated normal BAR(1) tree and calculating the bias corrected
#' # standard normal CI for the autoregressive coefficients of the BAR(1) model
#' # Note that in this example (B=2) for speeding up the calculations.
#' # B must be set to 499 or more for calculation accuracy.
#' z <- bfa.tree.gen(15, 1, 1, 1, -0.9, -0.9, 0, 10, c(-0.5))
#' bfa.ls.bc.ci(z, p=1, method="boot1", B=2)
bfa.ls.bc.ci <- function(z, p, method="BCa", conf.int="standard.normal.bc", conf.level=0.95, B=5, burn=5){
  if (method=="boot1"){
    coef.ls <- bfa.ls(z, p, error.cor = FALSE)
    boot.ci=matrix(NA, nrow=p,ncol=2)
    perc.boot1.ci=matrix(NA, nrow=p,ncol=2)
    percbc.ci=matrix(NA, nrow=p,ncol=2)
    a1.ls.boot1.star <- bfa.boot.ls.bc(z, p, method=method, burn = burn, B=B, boot.est=TRUE)$boot.bcest
    coef.ls.bc <- 2*coef.ls$coef[-1] - colMeans(a1.ls.boot1.star)
    if (conf.int=="All"){
      a1.ls.star = bfa.boot1.ls(z, p, burn = burn, B=B, boot.est = TRUE, boot.data = FALSE)
    for (i in 1:p) {
      boot.ci[i,] <- bfa.boot.ci(coef.ls.bc[i], a1.ls.boot1.star[,i], conf.level)
      perc.boot1.ci[i,] <- bfa.perc.ci(a1.ls.boot1.star[,i], conf.level)
      percbc.ci[i,] <- bfa.perc.bc.ci(z, coef.ls$coef[1,(i+1)], a1.ls.star$boot.est[,(i+1)], conf.level)
    }
    out <- list(Bias.corrected.coef=as.matrix(coef.ls.bc), standard.normal.bc.ci=boot.ci,percentile.ci=perc.boot1.ci, percentile.bc.ci=percbc.ci)
    colnames(out$Bias.corrected.coef) <- c(paste0("Estimates"))
    rownames(out$Bias.corrected.coef) <- c(paste0('X_[t/',2^c(1:p),"]"))
    colnames(out$standard.normal.bc.ci) <- c(paste0(((1-conf.level)/2)*100,"%"), paste0((1-(1-conf.level)/2)*100,"%"))
    rownames(out$standard.normal.bc.ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    colnames(out$percentile.ci) <- c(paste0(((1-conf.level)/2)*100,"%"), paste0((1-(1-conf.level)/2)*100,"%"))
    rownames(out$percentile.ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    colnames(out$percentile.bc.ci) <- c(paste0(((1-conf.level)/2)*100,"%"), paste0((1-(1-conf.level)/2)*100,"%"))
    rownames(out$percentile.bc.ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    }
    if (conf.int=="standard.normal.bc"){
      for (i in 1:p) {
        boot.ci[i,] <- bfa.boot.ci(coef.ls.bc[i], a1.ls.boot1.star[,i], conf.level)
      }
      out <- list(Bias.corrected.coef=as.matrix(coef.ls.bc),standard.normal.ci=boot.ci)
      colnames(out$standard.normal.ci) <- c(paste0(((1-conf.level)/2)*100,"%"), paste0((1-(1-conf.level)/2)*100,"%"))
      rownames(out$standard.normal.ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    }
    if (conf.int=="percentile"){
      for (i in 1:p) {
        perc.boot1.ci[i,] <- bfa.perc.ci(a1.ls.boot1.star[,i], conf.level)
      }
      out <- list(Bias.corrected.coef=as.matrix(coef.ls.bc),percentile.ci=perc.boot1.ci)
      colnames(out$Bias.corrected.coef) <- c(paste0("Estimates"))
      rownames(out$Bias.corrected.coef) <- c(paste0('X_[t/',2^c(1:p),"]"))
      colnames(out$percentile.ci) <- c(paste0(((1-conf.level)/2)*100,"%"), paste0((1-(1-conf.level)/2)*100,"%"))
      rownames(out$percentile.ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    }
    if (conf.int=="percentile.bc"){
      a1.ls.star = bfa.boot1.ls(z, p, burn = burn, B=B, boot.est = TRUE, boot.data = FALSE)
      for (i in 1:p) {
        percbc.ci[i,] <- bfa.perc.bc.ci(z, coef.ls$coef[1,(i+1)], a1.ls.star$boot.est[,(i+1)], conf.level)
      }
      out <- list(Bias.corrected.coef=as.matrix(coef.ls.bc), percentile.bc.ci=percbc.ci)
      colnames(out$Bias.corrected.coef) <- c(paste0("Estimates"))
      rownames(out$Bias.corrected.coef) <- c(paste0('X_[t/',2^c(1:p),"]"))
      colnames(out$percentile.bc.ci) <- c(paste0(((1-conf.level)/2)*100,"%"), paste0((1-(1-conf.level)/2)*100,"%"))
      rownames(out$percentile.bc.ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    }
  }
  if (method=="boot2fast"){
    coef.ls <- bfa.ls(z, p, error.cor = FALSE)
    boot2fast.ci=matrix(NA, nrow=p,ncol=2)
    perc.boot2fast.ci=matrix(NA, nrow=p,ncol=2)
    a1.ls.boot2fast.star <- bfa.boot.ls.bc(z, p, method=method, burn = burn, B=B, boot.est=TRUE)$boot.bcest
    bt2 <- bfa.boot2fast.ls(z, p, burn = burn, B=B)
    bias1 <- colMeans(bt2[[1]])[-1]-coef.ls$coef[-1]
    bt2.1 <- matrix(unlist(lapply(bt2[[2]], FUN=unlist)), nrow=B, byrow=TRUE)
    bias2 <- bias1 - colMeans(bt2.1)[-1] + colMeans(bt2[[1]])[-1]
    coef.ls.bc <- coef.ls$coef[-1] - bias1 - bias2
    if (conf.int=="All"){
      a1.ls.star = bfa.boot2.ls(z, p, burn = burn, B1=B, B2=1)
      for (i in 1:p) {
        boot2fast.ci[i,] <- bfa.boot.ci(coef.ls.bc[i], a1.ls.boot2fast.star[,i], conf.level)
        perc.boot2fast.ci[i,] <- bfa.perc.ci(a1.ls.boot2fast.star[,i], conf.level)
      }
      out <- list(Bias.corrected.coef=as.matrix(coef.ls.bc),standard.normal.bc.ci=boot2fast.ci,percentile.bc.ci=perc.boot2fast.ci)
      colnames(out$Bias.corrected.coef) <- c(paste0("Estimates"))
      rownames(out$Bias.corrected.coef) <- c(paste0('X_[t/',2^c(1:p),"]"))
      colnames(out$standard.normal.bc.ci) <- c(paste0(((1-conf.level)/2)*100,"%"), paste0((1-(1-conf.level)/2)*100,"%"))
      rownames(out$standard.normal.bc.ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
      colnames(out$percentile.bc.ci) <- c(paste0(((1-conf.level)/2)*100,"%"), paste0((1-(1-conf.level)/2)*100,"%"))
      rownames(out$percentile.bc.ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    }
    if (conf.int=="standard.normal.bc"){
      for (i in 1:p) {
        boot2fast.ci[i,] <- bfa.boot.ci(coef.ls.bc[i], a1.ls.boot2fast.star[,i], conf.level)
      }
      out <- list(Bias.corrected.coef=as.matrix(coef.ls.bc), standard.normal.bc.ci=boot2fast.ci)
      colnames(out$Bias.corrected.coef) <- c(paste0("Estimates"))
      rownames(out$Bias.corrected.coef) <- c(paste0('X_[t/',2^c(1:p),"]"))
      colnames(out$standard.normal.bc.ci) <- c(paste0(((1-conf.level)/2)*100,"%"), paste0((1-(1-conf.level)/2)*100,"%"))
      rownames(out$standard.normal.bc.ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    }
    if (conf.int=="percentile.bc"){
      for (i in 1:p) {
        perc.boot2fast.ci[i,] <- bfa.perc.ci(a1.ls.boot2fast.star[,i], conf.level)
      }
      out <- list(Bias.corrected.coef=as.matrix(coef.ls.bc),percentile.bc.ci=perc.boot2fast.ci)
      colnames(out$Bias.corrected.coef) <- c(paste0("Estimates"))
      rownames(out$Bias.corrected.coef) <- c(paste0('X_[t/',2^c(1:p),"]"))
      colnames(out$percentile.bc.ci) <- c(paste0(((1-conf.level)/2)*100,"%"), paste0((1-(1-conf.level)/2)*100,"%"))
      rownames(out$percentile.bc.ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    }
  }
  if (method=="BCa"){
    coef.ls <- bfa.ls(z, p, error.cor = FALSE)
    a1.ls.boot1.star <- bfa.boot.ls.bc(z, p, burn = burn, B=B, boot.est=TRUE)$boot.bcest
    coef.ls.bc <- 2*coef.ls$coef[-1] - colMeans(a1.ls.boot1.star)
    bca.ci <- boot.bca.ci(z,p,B,J=4,conf.level)$BCa
    out <- list(Bias.corrected.coef=as.matrix(coef.ls.bc),BCa.ci=bca.ci)
    colnames(out$Bias.corrected.coef) <- c(paste0("Estimates"))
    rownames(out$Bias.corrected.coef) <- c(paste0('X_[t/',2^c(1:p),"]"))
    colnames(out$BCa.ci) <- c(paste0(((1-conf.level)/2)*100,"%"), paste0((1-(1-conf.level)/2)*100,"%"))
    rownames(out$BCa.ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
  }
  return(out)
}
