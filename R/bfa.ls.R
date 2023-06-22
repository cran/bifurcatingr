#' Least Squares Estimation of Bifurcating Autoregressive Models
#'
#' This function performs Least Squares estimation of bifurcating autoregressive
#' (BFA) models of any order as described in Zhou and Basawa (2005).
#'
#' @param z a numeric vector containing the tree data
#' @param p an integer determining the order of bifurcating autoregressive model
#'   to be fit to the data
#' @param x.data a logical that determines whether the x data used in fitting
#'   the model should be returned. Defaults to FALSE.
#' @param y.data a logical that determines whether the y data used in fitting
#'   the model should be returned. Defaults to FALSE.
#' @param resids a logical that determines whether the model residuals should be
#'   returned. Defaults to FALSE.
#' @param error.cor a logical that determines whether the estimated correlation
#'   between pairs of model errors \eqn{(e_{2t}, e_{2t+1})} should be returned.
#'   Defaults to TRUE.
#' @param error.var a logical that determines whether the estimated variance of
#'   the model errors should be returned. Defaults to FALSE.
#' @param cov.matrix a logical that determines whether the estimated
#'   variance-covariance matrix of the least squares estimates should be
#'   returned. Defaults to FALSE.
#' @param conf a logical that determines whether confidence intervals for model
#'   coefficients should be returned. Defaults to FALSE. If TRUE, asymptotic normal
#'   confidence intervals for the intercept and the slops are calculated using
#'  \code{cov.matrix}. In addition, normal bootstrap confidence interval, and
#'   percentile confidence interval for the slop are calculated.
#' @param conf.level confidence level to be used in computing the normal
#'   confidence intervals for model coefficients when \code{conf=TRUE}. Defaults
#'   to \code{0.95}.
#' @param B number of bootstrap samples (replicates)
#' @param p.value a logical that determines whether p-values for model
#'   coefficients should be returned. Defaults to FALSE. If TRUE, p-values are
#'   computed from normal distribution using estimated coefficients and
#'   \code{cov.matrix}.
#' @return \item{coef}{a matrix containing the least squares estimates of the
#'   autoregressive coefficients} \item{error.cor}{the least squares estimate of
#'   the correlation between pairs of model errors \eqn{(e_{2t}, e_{2t+1})}.
#'   Only returned if \code{error.cor=TRUE}} \item{x}{a matrix containing the x
#'   data used in fitting the model. Only returned if \code{x.data=TRUE}}
#'   \item{y}{a vector containing the y data used in fitting the model. Only
#'   returned if \code{y.data=TRUE}} \item{resids}{the model residuals. Only
#'   returned if \code{resids=TRUE}} \item{error.var}{the estimated variance of
#'   the model errors. Only returned if \code{error.var=TRUE}}
#'   \item{cov.matrix}{the estimated variance-covariance matrix of the least
#'   squares coefficients. Only returned if \code{cov.matrix=TRUE}}\item{conf}{a
#'   matrix of normal confidence intervals for model coefficients. Only returned
#'   if \code{conf=TRUE}}\item{p.value}{a matrix of two-sided p-values for
#'   testing the significance of model coefficients. Computed from normal
#'   distribution and using the estimated covariance matrix \code{cov.matrix}.
#'   Only returned if \code{p.value=TRUE}}
#' @references Zhou, J. & Basawa, I. V. (2005). Least squares estimation for
#'   bifurcating autoregressive processes. \emph{Statistics & Probability
#'   Letters}, 74(1):77-88.
#' @export
#' @examples
#' z <- bfa.tree.gen(127, 1, 1, 1, -0.9, -0.9, 0, 10, c(0.7))
#' bfa.ls(z, p=1)
#' bfa.ls(z,p=1,conf=TRUE,cov.matrix = TRUE,conf.level = 0.9,p.value=TRUE)
bfa.ls <- function(z, p, x.data=FALSE, y.data=FALSE, resids=FALSE, error.cor=TRUE, error.var=FALSE, cov.matrix=FALSE, conf=FALSE, conf.level=0.95, B = 49, p.value=FALSE){
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
  model <- stats::lm(y~x,na.action="na.exclude")
  coef <- t(as.vector(model$coefficients))
  colnames(coef) <- c('Intercept', paste0('X_[t/',2^c(1:p),"]"))
  r1 <- numeric(length(stats::resid(model))/2)
  for (i in 1:length(r1)) r1[i] <- stats::resid(model)[i*2-1]
  r2 <- numeric(length(stats::resid(model))/2)
  for (i in 1:length(r2)) r2[i] <- stats::resid(model)[(i*2)]
  sigma2 <- sum((model$residuals)^2)/(n-(2^p)-p)
  out <- list(coef=coef)
  if(p.value){
    x.matrix <- cbind(rep(1,length(y)),x)
    a <- (1/(length(y)))*t(x.matrix)%*%x.matrix
    coef.cov <- sigma2*(1+error.cor)*MASS::ginv(a)
    coef.pval <- 2*stats::pnorm(abs(coef/sqrt(diag(coef.cov)/n)), lower.tail = FALSE)
    out$p.value <- coef.pval
  }
  if (x.data)
    out$x <- x
  if(y.data)
    out$y <- y
  if(resids)
    out$resids <- stats::resid(model)
  if(error.cor)
    out$error.cor <- sum(r1*r2,na.rm=TRUE)/(sigma2*m)
  if(error.var)
    out$error.var <- sigma2
  if(cov.matrix){
    x.matrix <- cbind(rep(1,length(y)),x)
    a <- (1/(length(y)))*t(x.matrix)%*%x.matrix
    coef.cov <- sigma2*(1+error.cor)*MASS::ginv(a)
    out$cov.matrix <- coef.cov
    colnames(out$cov.matrix) <- c('Intercept', paste0('X_[t/',2^c(1:p),"]"))
    rownames(out$cov.matrix) <- c('Intercept', paste0('X_[t/',2^c(1:p),"]"))
  }
  if(conf){
    x.matrix <- cbind(rep(1,length(y)),x)
    a <- (1/(length(y)))*t(x.matrix)%*%x.matrix
    coef.cov <- sigma2*(1+error.cor)*MASS::ginv(a)
    coef.conf.ll <- coef + stats::qnorm((1-conf.level)/2)*sqrt(diag(coef.cov)/n)
    coef.conf.ul <- coef - stats::qnorm((1-conf.level)/2)*sqrt(diag(coef.cov)/n)
    coef.conf <- t(rbind(coef.conf.ll, coef.conf.ul))
    colnames(coef.conf) <- c(paste0(((1-conf.level)/2)*100,"%"), paste0((1-(1-conf.level)/2)*100,"%"))
    out$asymptotic.ci <- coef.conf

    asyb.matrs=matrix(rep(NA,(p+1)*2),nrow=p+1)
    perc.matrs=matrix(rep(NA,(p+1)*2),nrow=p+1)
    a1.ls.star = bfa.boot1.ls(z, p, burn = 5, B , boot.est = TRUE, boot.data = FALSE)
    for (i in 1:(p+1)) {
      asyb.ci = bfa.boot.ci(coef[1,i], a1.ls.star$boot.est[,i], conf.level)
      asyb.matrs[i,1]=asyb.ci[1]
      asyb.matrs[i,2]=asyb.ci[2]
      perc.ci = bfa.perc.ci(a1.ls.star$boot.est[,i], conf.level)
      perc.matrs[i,1]=perc.ci[1]
      perc.matrs[i,2]=perc.ci[2]

    }

    out$bootstarp.ci <- asyb.matrs
    colnames(out$bootstarp.ci) <- c(paste0(((1-conf.level)/2)*100,"%"), paste0((1-(1-conf.level)/2)*100,"%"))
    rownames(out$bootstarp.ci) <- c('Intercept', paste0('X_[t/',2^c(1:p),"]"))
    out$percentile.ci <- perc.matrs
    colnames(out$percentile.ci) <- c(paste0(((1-conf.level)/2)*100,"%"), paste0((1-(1-conf.level)/2)*100,"%"))
    rownames(out$percentile.ci) <- c('Intercept', paste0('X_[t/',2^c(1:p),"]"))
    }
  return(out)
}
