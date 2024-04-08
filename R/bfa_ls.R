#' Least Squares Estimation of Bifurcating Autoregressive Models
#'
#' This function performs Least Squares estimation of bifurcating autoregressive
#' (BFA) models of any order as described in Zhou and Basawa (2005).
#'
#' @param z a numeric vector containing the tree data
#' @param p an integer determining the order of bifurcating autoregressive model
#'   to be fit to the data
#' @param x_data a logical that determines whether the x data used in fitting
#'   the model should be returned. Defaults to FALSE.
#' @param y_data a logical that determines whether the y data used in fitting
#'   the model should be returned. Defaults to FALSE.
#' @param resids a logical that determines whether the model residuals should be
#'   returned. Defaults to FALSE.
#' @param error_cor a logical that determines whether the estimated correlation
#'   between pairs of model errors \eqn{(e_{2t}, e_{2t+1})} should be returned.
#'   Defaults to TRUE.
#' @param error_var a logical that determines whether the estimated variance of
#'   the model errors should be returned. Defaults to FALSE.
#' @param cov_matrix a logical that determines whether the estimated
#'   variance-covariance matrix of the least squares estimates should be
#'   returned. Defaults to FALSE.
#' @param conf a logical that determines whether confidence intervals for model
#'   coefficients should be returned. Defaults to FALSE. If TRUE, asymptotic normal
#'   confidence intervals for the intercept and the slops are calculated using
#'  \code{cov_matrix}. In addition, normal bootstrap confidence interval, and
#'   percentile confidence interval for the slop are calculated.
#' @param conf_level confidence level to be used in computing the normal
#'   confidence intervals for model coefficients when \code{conf=TRUE}. Defaults
#'   to \code{0.95}.
#' @param B number of bootstrap samples (replicates)
#' @param p_value a logical that determines whether p-values for model
#'   coefficients should be returned. Defaults to FALSE. If TRUE, p-values are
#'   computed from normal distribution using estimated coefficients and
#'   \code{cov_matrix}.
#' @return \item{coef}{a matrix containing the least squares estimates of the
#'   autoregressive coefficients} \item{error_cor}{the least squares estimate of
#'   the correlation between pairs of model errors \eqn{(e_{2t}, e_{2t+1})}.
#'   Only returned if \code{error_cor=TRUE}} \item{x}{a matrix containing the x
#'   data used in fitting the model. Only returned if \code{x_data=TRUE}}
#'   \item{y}{a vector containing the y data used in fitting the model. Only
#'   returned if \code{y_data=TRUE}} \item{resids}{the model residuals. Only
#'   returned if \code{resids=TRUE}} \item{error_var}{the estimated variance of
#'   the model errors. Only returned if \code{error_var=TRUE}}
#'   \item{cov_matrix}{the estimated variance-covariance matrix of the least
#'   squares coefficients. Only returned if \code{cov_matrix=TRUE}}\item{conf}{a
#'   matrix of normal confidence intervals for model coefficients. Only returned
#'   if \code{conf=TRUE}}\item{p_value}{a matrix of two-sided p-values for
#'   testing the significance of model coefficients. Computed from normal
#'   distribution and using the estimated covariance matrix \code{cov_matrix}.
#'   Only returned if \code{p_value=TRUE}}
#' @references Zhou, J. & Basawa, I. V. (2005). Least squares estimation for
#'   bifurcating autoregressive processes. \emph{Statistics & Probability
#'   Letters}, 74(1):77-88.
#' @export
#' @examples
#' z <- bfa_tree_gen(127, 1, 1, 1, -0.9, -0.9, 0, 10, c(0.7))
#' bfa_ls(z, p=1)
#' bfa_ls(z,p=1,conf=TRUE,cov_matrix = TRUE,conf_level = 0.9,p_value=TRUE)
bfa_ls <- function(z, p, x_data=FALSE, y_data=FALSE, resids=FALSE, error_cor=TRUE, error_var=FALSE, cov_matrix=FALSE, conf=FALSE, conf_level=0.95, B = 49, p_value=FALSE){
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
  if(p_value){
    x_matrix <- cbind(rep(1,length(y)),x)
    a <- (1/(length(y)))*t(x_matrix)%*%x_matrix
    coef_cov <- sigma2*(1+error_cor)*MASS::ginv(a)
    coef_pval <- 2*stats::pnorm(abs(coef/sqrt(diag(coef_cov)/n)), lower.tail = FALSE)
    out$p_value <- coef_pval
  }
  if (x_data)
    out$x <- x
  if(y_data)
    out$y <- y
  if(resids)
    out$resids <- stats::resid(model)
  if(error_cor)
    out$error_cor <- sum(r1*r2,na.rm=TRUE)/(sigma2*m)
  if(error_var)
    out$error_var <- sigma2
  if(cov_matrix){
    x_matrix <- cbind(rep(1,length(y)),x)
    a <- (1/(length(y)))*t(x_matrix)%*%x_matrix
    coef_cov <- sigma2*(1+error_cor)*MASS::ginv(a)
    out$cov_matrix <- coef_cov
    colnames(out$cov_matrix) <- c('Intercept', paste0('X_[t/',2^c(1:p),"]"))
    rownames(out$cov_matrix) <- c('Intercept', paste0('X_[t/',2^c(1:p),"]"))
  }
  if(conf){
    x_matrix <- cbind(rep(1,length(y)),x)
    a <- (1/(length(y)))*t(x_matrix)%*%x_matrix
    coef_cov <- sigma2*(1+error_cor)*MASS::ginv(a)
    coef_conf_ll <- coef + stats::qnorm((1-conf_level)/2)*sqrt(diag(coef_cov)/n)
    coef_conf_ul <- coef - stats::qnorm((1-conf_level)/2)*sqrt(diag(coef_cov)/n)
    coef_conf <- t(rbind(coef_conf_ll, coef_conf_ul))
    colnames(coef_conf) <- c(paste0(((1-conf_level)/2)*100,"%"), paste0((1-(1-conf_level)/2)*100,"%"))
    out$asymptotic_ci <- coef_conf

    asyb_matrs <- matrix(rep(NA,(p+1)*2),nrow=p+1)
    perc_matrs <- matrix(rep(NA,(p+1)*2),nrow=p+1)
    a1_ls_star <- bfa_boot1_ls(z, p, burn = 5, B , boot_est = TRUE, boot_data = FALSE)
    for (i in 1:(p+1)) {
      asyb_ci <- bfa_boot_ci(coef[1,i], a1_ls_star$boot_est[,i], conf_level)
      asyb_matrs[i,1] <- asyb_ci[1]
      asyb_matrs[i,2] <- asyb_ci[2]
      perc_ci <- bfa_perc_ci(a1_ls_star$boot_est[,i], conf_level)
      perc_matrs[i,1] <- perc_ci[1]
      perc_matrs[i,2] <- perc_ci[2]

    }

    out$bootstarp_ci <- asyb_matrs
    colnames(out$bootstarp_ci) <- c(paste0(((1-conf_level)/2)*100,"%"), paste0((1-(1-conf_level)/2)*100,"%"))
    rownames(out$bootstarp_ci) <- c('Intercept', paste0('X_[t/',2^c(1:p),"]"))
    out$percentile_ci <- perc_matrs
    colnames(out$percentile_ci) <- c(paste0(((1-conf_level)/2)*100,"%"), paste0((1-(1-conf_level)/2)*100,"%"))
    rownames(out$percentile_ci) <- c('Intercept', paste0('X_[t/',2^c(1:p),"]"))
    }
  out
}
